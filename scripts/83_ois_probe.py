#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 83_ois_probe.py
#
# Generation-readiness probe: dynamic range, motif insertion, ISM. GPU STEP, ~1.7 h.
#
# CONSUMES: results/head_ois.pt, pred_*.npy, /workspace/genomes/hg38/
# PRODUCES: results/{results_probe.json,motif_insertion.csv,ism_importance.npz,ism_per_region.csv,locus_checks.csv,PROBE_LOG.txt}
# ---------------------------------------------------------------------------
"""
STEP 3 -- GENERATION-READINESS PROBE (the decision gate).

Runs ONLY if the shuffled-label control passed. Classification skill and design skill are
different questions: a head can rank existing elements from trunk features without the model
being able to say what to BUILD.

COORDINATE SYSTEM: this probe uses design_anchor_hg38 -- the ATAC summit, a median 245 bp
from the H3K27ac summit. H3K27ac marks the nucleosomes FLANKING an element and is depleted
over the nucleosome-free region where TFs bind, so a 200 bp cassette centred on a K27ac
summit would sit in the acetylation trough, on a nucleosome. Training and caching used
region_center_hg38; this probe does not.

Prediction for a modified sequence = head applied to the trunk embedding mean-pooled over
the central 1 kb (POOL_W) of the predicted span -- the same width the head was trained on.
"""
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import torch
from scipy import stats

import grelu.resources
import grelu.sequence.format
from grelu.model.heads import ConvHead

torch.backends.cudnn.allow_tf32 = False
torch.backends.cuda.matmul.allow_tf32 = False

DEV = "cuda"
CACHE = Path("/workspace/ois_cache")
OUT = Path("/workspace/ois_results")
TRAINSET = "/workspace/output/ois_enhancer_trainset/ois_enhancer_trainset.h5ad"

SEQ_LEN = 524_288
PRED_SPAN = 196_608
N_BINS = 6144
BIN = 32
SPAN_OFF = (SEQ_LEN - PRED_SPAN) // 2
POOL_W = 1000                 # head was trained on 1 kb pooled embeddings
CASSETTE = 200                # the design target length
ISM_W = 500                   # ISM scan width (PWMs find matches at this width; see report)
N_ISM_REGIONS = 20
ISM_BUDGET_S = 6000

JASPAR = {"NFKB_RELA": "MA0107.1", "CEBPB": "MA0466.2",
          "AP1_FOSJUN": "MA0099.3", "ETS1": "MA0098.3"}
MOTIF_SEQ = {"NFKB_RELA": "GGGACTTTCC", "CEBPB": "TTGCGCAA",
             "AP1_FOSJUN": "TGACTCA", "ETS1": "AGGAAGT"}
K_VALUES = [0, 1, 2, 3, 4, 6, 8]
N_SEEDS = 5
BG = {"A": .295, "C": .205, "G": .205, "T": .295}
REL_THRESH = 0.80

probe = {"coordinate_system": "design_anchor_hg38 (ATAC summit), NOT region_center_hg38",
         "pool_width": POOL_W, "cassette_len": CASSETTE, "ism_width": ISM_W,
         "replicates": "n=2 vs 2 in both studies"}
LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


tr = json.load(open(OUT / "results_train.json"))
if not tr.get("control_gate_passed"):
    R("CONTROL GATE NOT PASSED -- probe must not run. Exiting.")
    sys.exit(1)
R("control gate passed; running probe")

ad = anndata.read_h5ad(TRAINSET)
var = ad.var
TASKS = list(ad.obs.index)
IMR_PRO = np.asarray(ad.X, np.float32)[TASKS.index("IMR90_SEN_vs_PRO")]
GM21 = np.asarray(ad.X, np.float32)[TASKS.index("GM21_RAS_vs_EV")]

head = ConvHead(n_tasks=1, in_channels=1920, act_func=None, pool_func="avg",
                norm=False).to(DEV)
head.load_state_dict(torch.load(OUT / "head_ois.pt", map_location=DEV))
head.eval()
model = grelu.resources.load_model(repo_id="Genentech/borzoi-model",
                                   filename="human_rep0.ckpt")
model.eval()
core = model.model.to("cuda")
R("head + trunk loaded")

ALLOWED = set(grelu.sequence.format.ALLOWED_BASES)
_SAN = {c: "N" for c in map(chr, range(256)) if c.upper() not in ALLOWED}
BASES = ["A", "C", "G", "T"]
CB0 = (PRED_SPAN // 2 - POOL_W // 2) // BIN
CB1 = -((-(PRED_SPAN // 2 + POOL_W // 2)) // BIN)


def onehot(s):
    if set(s) - ALLOWED:
        s = s.translate(str.maketrans(_SAN))
    x = grelu.sequence.format.convert_input_type([s], output_type="one_hot")
    if not torch.is_tensor(x):
        x = torch.as_tensor(np.asarray(x))
    x = x.float()
    return x if x.dim() == 3 else x.unsqueeze(0)


def pick_batch():
    for b in (8, 6, 4, 3, 2, 1):
        try:
            x = torch.zeros(b, 4, SEQ_LEN, device=DEV)
            with torch.no_grad():
                core.embedding(x)
            del x; torch.cuda.empty_cache(); return b
        except RuntimeError:
            torch.cuda.empty_cache()
    return 1


BATCH = pick_batch()
R(f"trunk batch size {BATCH} (peak {torch.cuda.max_memory_allocated()/1e9:.2f} GB)")
probe["batch"] = BATCH


def score(seqs):
    """List[str] -> predicted H3K27ac response, pooled over the central POOL_W."""
    out = []
    for b0 in range(0, len(seqs), BATCH):
        x = torch.cat([onehot(s) for s in seqs[b0:b0 + BATCH]]).to(DEV)
        with torch.no_grad():
            t = core.embedding(x)[:, :, CB0:CB1].mean(dim=-1, keepdim=True)
            out.append(head(t).squeeze(-1).squeeze(-1).cpu().numpy())
        del x
    return np.concatenate(out)


def fetch(chrom, start, end):
    iv = pd.DataFrame({"chrom": [chrom], "start": [int(start)], "end": [int(end)],
                       "strand": ["+"]})
    s = grelu.sequence.format.convert_input_type(iv, output_type="strings", genome="hg38")
    s = s[0] if isinstance(s, list) else s
    return s.translate(str.maketrans(_SAN)) if (set(s) - ALLOWED) else s


sizes = {}
for line in open("/workspace/genomes/hg38/hg38.fa.sizes"):
    a, b = line.split(); sizes[a] = int(b)


def window_on(anchor_chrom, anchor_pos):
    ws = int(anchor_pos) - SEQ_LEN // 2
    ws = max(0, min(ws, sizes[anchor_chrom] - SEQ_LEN))
    return fetch(anchor_chrom, ws, ws + SEQ_LEN), ws


# =====================================================================================
# (c) DYNAMIC RANGE + SIGN AGREEMENT  (cheap, from the cached predictions)
# =====================================================================================
R("\n" + "=" * 78 + "\n(c) DYNAMIC RANGE AND SIGN AGREEMENT\n" + "=" * 78)
pr = np.load(OUT / "pred_rest.npy"); ir = np.load(OUT / "idx_rest.npy")
pt = np.load(OUT / "pred_test.npy"); it = np.load(OUT / "idx_test.npy")
allp = np.concatenate([pr, pt]); alli = np.concatenate([ir, it])
meas = IMR_PRO[alli]
qs = [0, 1, 5, 25, 50, 75, 95, 99, 100]
pcs = {str(q): float(np.percentile(allp, q)) for q in qs}
R(f"  predicted: sd {allp.std():.4f}  range [{allp.min():+.4f}, {allp.max():+.4f}]  "
  f"IQR {np.subtract(*np.percentile(allp,[75,25])):.4f}")
R(f"  measured IMR90 (n=2v2): sd {meas.std():.4f}  range [{meas.min():+.4f}, {meas.max():+.4f}]")
R(f"  predicted sd / measured sd = {allp.std()/meas.std():.4f}")
R(f"  predicted percentiles: " + ", ".join(f"{k}%={v:+.3f}" for k, v in pcs.items()))
need_95 = float(np.percentile(allp, 95) - np.percentile(allp, 50))
R(f"  delta required to move median -> 95th percentile: {need_95:+.4f}")
top = np.argsort(-np.abs(meas))[:50]
sign_ok = int((np.sign(allp[top]) == np.sign(meas[top])).sum())
R(f"  sign agreement on the 50 most responsive held-out regions: {sign_ok}/50 "
  f"({100*sign_ok/50:.0f}%)")
sign20 = int((np.sign(allp[np.argsort(-np.abs(meas))[:20]]) ==
              np.sign(meas[np.argsort(-np.abs(meas))[:20]])).sum())
R(f"  sign agreement on the top 20: {sign20}/20 ({100*sign20/20:.0f}%)")
probe["dynamic_range"] = dict(pred_sd=float(allp.std()), pred_min=float(allp.min()),
                              pred_max=float(allp.max()), percentiles=pcs,
                              measured_sd=float(meas.std()),
                              pred_over_measured_sd=float(allp.std() / meas.std()),
                              delta_median_to_p95=need_95,
                              sign_agreement_top50=sign_ok,
                              sign_agreement_top20=sign20, n_regions=int(len(allp)))

# =====================================================================================
# (a) MOTIF-INSERTION DOSE RESPONSE
# =====================================================================================
R("\n" + "=" * 78 + "\n(a) MOTIF-INSERTION DOSE RESPONSE (design_anchor coordinates)\n" + "=" * 78)
cand = np.where(np.isfinite(meas) & (np.abs(meas) < 0.1))[0]
neutral_i = alli[cand[np.argmin(np.abs(meas[cand]))]] if len(cand) else alli[0]
nrow = var.iloc[neutral_i]
anchor = int(nrow["design_anchor_hg38"]) if nrow["design_anchor_hg38"] > 0 else int(nrow["region_center_hg38"])
R(f"  neutral background context: {nrow['chrom']}:{anchor:,} "
  f"(design_anchor; measured IMR90 {IMR_PRO[neutral_i]:+.4f})")
bgseq, ws = window_on(nrow["chrom"], anchor)
lo = SEQ_LEN // 2 - CASSETTE // 2
rng = np.random.default_rng(7)


def scramble(m, s):
    a = list(m); np.random.default_rng(s).shuffle(a); return "".join(a)


rows = []
for mname, mseq in MOTIF_SEQ.items():
    for variant in ("real", "scrambled"):
        for seed in range(N_SEEDS):
            bg = "".join(rng.choice(list("ACGT"), size=CASSETTE,
                                    p=[BG["A"], BG["C"], BG["G"], BG["T"]]))
            seqs = []
            for k in K_VALUES:
                cas = list(bg)
                use = mseq if variant == "real" else scramble(mseq, seed * 31 + k)
                L = len(use)
                if k:
                    span = CASSETTE // k
                    for j in range(k):
                        st = max(0, min(CASSETTE - L, j * span + (span - L) // 2))
                        cas[st:st + L] = list(use)
                seqs.append(bgseq[:lo] + "".join(cas) + bgseq[lo + CASSETTE:])
            sc = score(seqs)
            for k, v in zip(K_VALUES, sc):
                rows.append(dict(motif=mname, variant=variant, seed=seed, k=int(k),
                                 pred=float(v)))
ins = pd.DataFrame(rows)
ins.to_csv(OUT / "motif_insertion.csv", index=False)
summ = {}
R(f"\n  {'motif':12s} {'variant':10s} {'rho(pred,k)':>12s} {'slope/copy':>12s} {'delta k0->8':>12s}")
for mname in MOTIF_SEQ:
    for variant in ("real", "scrambled"):
        d = ins[(ins.motif == mname) & (ins.variant == variant)]
        rho = stats.spearmanr(d["k"], d["pred"]).statistic
        sl = float(np.polyfit(d["k"], d["pred"], 1)[0])
        d0 = d[d.k == 0]["pred"].mean(); d8 = d[d.k == max(K_VALUES)]["pred"].mean()
        summ[f"{mname}_{variant}"] = dict(rho=float(rho), slope=sl, k0=float(d0),
                                          kmax=float(d8), delta=float(d8 - d0))
        R(f"  {mname:12s} {variant:10s} {rho:12.3f} {sl:12.5f} {d8-d0:12.5f}")
best = max((v["delta"] for k, v in summ.items() if k.endswith("_real")))
R(f"\n  THE NUMBER THAT MATTERS")
R(f"    largest real-motif delta at saturation (k=8): {best:+.5f}")
R(f"    prediction sd across held-out regions        : {allp.std():.5f}")
R(f"    delta required, median -> 95th percentile    : {need_95:+.5f}")
R(f"    lever length = achieved / required           : {best/need_95:.3f}x")
R(f"    (expression model for comparison: +0.132 achieved vs ~+0.65 required = 0.20x)")
probe["motif_insertion"] = dict(summary=summ, best_real_delta=float(best),
                                required_delta=need_95,
                                lever_ratio=float(best / need_95),
                                neutral_context=f"{nrow['chrom']}:{anchor}")

# =====================================================================================
# (b) IN-SILICO SATURATION MUTAGENESIS + PWM ENRICHMENT
# =====================================================================================
R("\n" + "=" * 78 + f"\n(b) ISM ({ISM_W} bp around design_anchor) + JASPAR PWM enrichment\n"
  + "=" * 78)


def fetch_pfm(mid):
    for ep in ["https://jaspar.elixir.no/api/v1/matrix/{}/?format=json",
               "https://jaspar.genereg.net/api/v1/matrix/{}/?format=json"]:
        try:
            r = subprocess.run(["curl", "-sS", "-f", "--max-time", "40", ep.format(mid)],
                               capture_output=True, text=True, timeout=60)
            if r.returncode == 0 and r.stdout.strip().startswith("{"):
                j = json.loads(r.stdout)
                if j.get("pfm"):
                    return j
        except Exception:
            continue
    return None


PWM, MINS, MAXS, srcs = {}, {}, {}, {}
for nm, mid in JASPAR.items():
    j = fetch_pfm(mid)
    if j is None:
        R(f"  JASPAR fetch FAILED for {nm} {mid}")
        continue
    pfm = np.array([j["pfm"][b] for b in "ACGT"], float) + 0.25
    pfm = pfm / pfm.sum(0, keepdims=True)
    w = np.log2(pfm / np.array([[BG[b]] for b in "ACGT"]))
    PWM[nm] = w; MINS[nm] = w.min(0).sum(); MAXS[nm] = w.max(0).sum()
    srcs[nm] = dict(matrix_id=j.get("matrix_id", mid), name=j.get("name"),
                    length=int(w.shape[1]))
    R(f"  {nm}: {j.get('matrix_id', mid)} ({j.get('name')}) length {w.shape[1]}")
probe["motif_sources"] = srcs
IDXB = {b: i for i, b in enumerate("ACGT")}


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def scan_mask(seq, key):
    w = PWM[key]; L = w.shape[1]
    cut = MINS[key] + REL_THRESH * (MAXS[key] - MINS[key])
    mask = np.zeros(len(seq), bool)
    for s, flip in ((seq, False), (revcomp(seq), True)):
        for i in range(len(s) - L + 1):
            sub = s[i:i + L]
            if "N" in sub:
                continue
            if sum(w[IDXB[b], j] for j, b in enumerate(sub)) >= cut:
                a, b_ = (len(seq) - i - L, len(seq) - i) if flip else (i, i + L)
                mask[a:b_] = True
    return mask


order = np.argsort(-np.abs(meas))
sel = []
for j in order:
    gi = alli[j]
    r = var.iloc[gi]
    if r["design_anchor_hg38"] > 0:
        sel.append(gi)
    if len(sel) >= N_ISM_REGIONS:
        break
R(f"  regions: top {len(sel)} held-out by |measured IMR90 response| with a design anchor")

per_region, imp_store = [], {}
t0 = time.time()
for gi in sel:
    if time.time() - t0 > ISM_BUDGET_S:
        R(f"  [time budget reached after {len(per_region)} regions]")
        break
    r = var.iloc[gi]
    anchor = int(r["design_anchor_hg38"])
    ref, ws = window_on(r["chrom"], anchor)
    ctr = anchor - ws
    lo_, hi_ = ctr - ISM_W // 2, ctr + ISM_W // 2
    win = ref[lo_:hi_]
    muts, meta = [], []
    for p in range(ISM_W):
        wt = win[p]
        for b in BASES:
            if b == wt:
                continue
            muts.append(ref[:lo_ + p] + b + ref[lo_ + p + 1:])
            meta.append(p)
    sc = score([ref] + muts)
    delta = sc[1:] - sc[0]
    imp = np.zeros(ISM_W)
    for p, d in zip(meta, delta):
        imp[p] += abs(d) / 3.0
    masks = {k: scan_mask(win, k) for k in PWM}
    anym = np.any(np.stack(list(masks.values())), axis=0) if masks else np.zeros(ISM_W, bool)
    top10 = int(np.ceil(0.10 * ISM_W))
    conc = float(np.sort(imp)[::-1][:top10].sum() / imp.sum()) if imp.sum() else np.nan
    enr = float(imp[anym].mean() / imp[~anym].mean()) if (anym.any() and (~anym).any()) else np.nan
    per_region.append(dict(region=var.index[gi], chrom=r["chrom"], anchor=anchor,
                           measured_IMR90=float(IMR_PRO[gi]), wt_pred=float(sc[0]),
                           max_abs_delta=float(np.abs(delta).max()),
                           mean_abs_delta=float(np.abs(delta).mean()),
                           top10_concentration=conc, motif_bases=int(anym.sum()),
                           motif_enrichment=enr))
    imp_store[str(var.index[gi])] = imp
    R(f"    {var.index[gi]:22s} wt {sc[0]:+.4f} max|d| {np.abs(delta).max():.5f} "
      f"conc {conc:.3f} motif_bases {int(anym.sum()):3d} enr "
      f"{enr:.3f}" if enr == enr else
      f"    {var.index[gi]:22s} wt {sc[0]:+.4f} (no motif match in window)")

pr_df = pd.DataFrame(per_region)
pr_df.to_csv(OUT / "ism_per_region.csv", index=False)
np.savez_compressed(OUT / "ism_importance.npz", **imp_store)
ev = pr_df["motif_enrichment"].dropna()
t_, p_ = stats.ttest_1samp(ev, 1.0) if len(ev) > 1 else (np.nan, np.nan)
R(f"\n  regions run: {len(pr_df)}")
R(f"  windows with >=1 PWM match: {int((pr_df['motif_bases']>0).sum())}/{len(pr_df)}")
R(f"  top-10% attribution concentration: mean {pr_df['top10_concentration'].mean():.3f} "
  f"(0.10 == perfectly diffuse; expression model gave 0.291)")
R(f"  PWM enrichment (in-motif / out-motif importance): mean "
  f"{ev.mean():.3f} +/- {ev.std(ddof=1) if len(ev)>1 else 0:.3f} over n={len(ev)}, "
  f"p vs 1.0 = {p_:.4f}")
R(f"  mean |delta| per single mutation: {pr_df['mean_abs_delta'].mean():.6f} "
  f"({100*pr_df['mean_abs_delta'].mean()/allp.std():.2f}% of prediction sd)")
probe["ism"] = dict(n_regions=len(pr_df), ism_width=ISM_W,
                    windows_with_match=int((pr_df["motif_bases"] > 0).sum()),
                    concentration_mean=float(pr_df["top10_concentration"].mean()),
                    concentration_null=0.10,
                    enrichment_mean=float(ev.mean()) if len(ev) else None,
                    enrichment_sd=float(ev.std(ddof=1)) if len(ev) > 1 else None,
                    enrichment_n=int(len(ev)),
                    enrichment_p=float(p_) if p_ == p_ else None,
                    mean_abs_delta=float(pr_df["mean_abs_delta"].mean()),
                    seconds=float(time.time() - t0))

json.dump(probe, open(OUT / "results_probe.json", "w"), indent=2, default=str)
(OUT / "PROBE_LOG.txt").write_text("\n".join(LOG))
R("\nPROBE_DONE")
