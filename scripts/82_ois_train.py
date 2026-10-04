#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 82_ois_train.py
#
# Train the ConvHead on cached embeddings; 200-permutation shuffled-label control. GPU STEP.
#
# CONSUMES: ois_cache/region_{fwd,rc}.npy, ois_enhancer_trainset.h5ad
# PRODUCES: results/{head_ois.pt,results_train.json,pred_*.npy,idx_*.npy,TRAIN_LOG.txt}
# ---------------------------------------------------------------------------
"""
STEP 2 -- train the OIS enhancer head with a cross-study holdout, plus controls.

Train on GM21 (skin, ER:RAS OIS, empty-vector controls, n=2 vs 2).
Evaluate on IMR90 (lung, HRasV12 OIS, n=2 vs 2) -- a different lab, tissue and induction
construct.

DELIBERATE DEVIATION, stated: the brief says "early stopping on validation loss" with IMR90
as validation. Early-stopping on IMR90 would select the checkpoint using the held-out study,
leaking it into training and inflating the primary metric. Early stopping therefore uses a
GM21-INTERNAL chromosome holdout (chr8/chr16/chr20); IMR90 is never touched until final
evaluation. chr9+chr6 remain held out from everything.

COORDINATE SYSTEM: region_center_hg38 throughout -- the H3K27ac measurement coordinate.

n=2 PER ARM IN BOTH STUDIES. Stated with every number below.
"""
import json
import sys
import time
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import torch
from scipy import stats
from sklearn.metrics import roc_auc_score

from grelu.model.heads import ConvHead

torch.backends.cudnn.allow_tf32 = False
torch.backends.cuda.matmul.allow_tf32 = False

DEV = "cuda"
CACHE = Path("/workspace/ois_cache")
OUT = Path("/workspace/ois_results")
OUT.mkdir(parents=True, exist_ok=True)
TRAINSET = "/workspace/output/ois_enhancer_trainset/ois_enhancer_trainset.h5ad"

HP = dict(lr=1e-4, batch_size=1024, max_epochs=200, patience=10, optimizer="adam",
          loss="mse", seed=0, rc_augment=True, eval_rc_average=True,
          early_stop_on="GM21-internal chr8/chr16/chr20 (NOT IMR90)")
LOG2FC_THRESH = 1.0
N_PERM = 200
CEILING = 0.637
ES_CHROMS = ["chr8", "chr16", "chr20"]

res = {"hyperparameters": HP, "log2fc_threshold": LOG2FC_THRESH, "n_permutations": N_PERM,
       "reproducibility_ceiling": CEILING,
       "replicates": "n=2 vs 2 in BOTH studies (GM21 and IMR90)"}
LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


# ---------------------------------------------------------------- data
ad = anndata.read_h5ad(TRAINSET)
var = ad.var
Y = np.asarray(ad.X, np.float32)
TASKS = list(ad.obs.index)
R(f"tasks: {TASKS}")
GM21 = Y[TASKS.index("GM21_RAS_vs_EV")]
IMR_PRO = Y[TASKS.index("IMR90_SEN_vs_PRO")]
IMR_QUI = Y[TASKS.index("IMR90_SEN_vs_QUI")]

Xf = torch.tensor(np.load(CACHE / "region_fwd.npy"), device=DEV)
Xr = torch.tensor(np.load(CACHE / "region_rc.npy"), device=DEV)
n_reg = Xf.shape[0]
assert n_reg == len(var), f"cache/trainset mismatch {n_reg} vs {len(var)}"
cached = (Xf.abs().sum(1) > 0).cpu().numpy()
R(f"cache {tuple(Xf.shape)}; regions with a non-empty embedding: {int(cached.sum()):,}")

chrom = var["chrom"].to_numpy()
is_test = np.isin(chrom, ["chr9", "chr6"])
is_es = np.isin(chrom, ES_CHROMS)
tr_idx = np.where(cached & ~is_test & ~is_es)[0]
es_idx = np.where(cached & ~is_test & is_es)[0]
te_idx = np.where(cached & is_test)[0]
rest_idx = np.where(cached & ~is_test)[0]
R(f"GM21 train regions {len(tr_idx):,} | GM21 early-stop {len(es_idx):,} | "
  f"chr9+chr6 held out {len(te_idx):,}")
res["region_counts"] = dict(train=len(tr_idx), early_stop=len(es_idx),
                            test_chr9_chr6=len(te_idx), non_test=len(rest_idx))


def make_head(seed):
    torch.manual_seed(seed)
    return ConvHead(n_tasks=1, in_channels=1920, act_func=None, pool_func="avg",
                    norm=False).to(DEV)


h = make_head(999)
with torch.no_grad():
    p = h(Xf[:512].unsqueeze(-1)).squeeze(-1)
n_neg = int((p < 0).sum())
R(f"\nhead signed-output check: range [{p.min():.4f}, {p.max():.4f}] negatives "
  f"{n_neg}/{p.numel()} ({100*n_neg/p.numel():.1f}%) -> {'PASS' if n_neg else 'FAIL'}")
res["signed_output_check"] = dict(pred_min=float(p.min()), pred_max=float(p.max()),
                                  n_negative=n_neg, verdict="PASS" if n_neg else "FAIL")
if not n_neg:
    json.dump(res, open(OUT / "results_train.json", "w"), indent=2); sys.exit(1)


def predict(head, idx):
    head.eval()
    gi = torch.as_tensor(idx, device=DEV)
    with torch.no_grad():
        a = head(Xf[gi].unsqueeze(-1)).squeeze(-1).squeeze(-1)
        if HP["eval_rc_average"]:
            b = head(Xr[gi].unsqueeze(-1)).squeeze(-1).squeeze(-1)
            return ((a + b) / 2).cpu().numpy()
        return a.cpu().numpy()


def train_head(target, seed=0, verbose=False, max_epochs=None):
    tgt = torch.tensor(target, device=DEV)
    head = make_head(seed)
    opt = torch.optim.Adam(head.parameters(), lr=HP["lr"])
    g = torch.Generator(); g.manual_seed(seed)
    best, best_state, bad, n_ep = np.inf, None, 0, 0
    me = max_epochs or HP["max_epochs"]
    for ep in range(me):
        head.train()
        perm = torch.randperm(len(tr_idx), generator=g).numpy()
        for b in range(0, len(tr_idx), HP["batch_size"]):
            bi = tr_idx[perm[b:b + HP["batch_size"]]]
            gi = torch.as_tensor(bi, device=DEV)
            if HP["rc_augment"]:
                m = torch.rand(len(bi), device=DEV) < 0.5
                emb = torch.where(m.unsqueeze(1), Xf[gi], Xr[gi])
            else:
                emb = Xf[gi]
            pred = head(emb.unsqueeze(-1)).squeeze(-1).squeeze(-1)
            loss = torch.nn.functional.mse_loss(pred, tgt[gi])
            opt.zero_grad(); loss.backward(); opt.step()
        head.eval()
        with torch.no_grad():
            gi = torch.as_tensor(es_idx, device=DEV)
            pv = (head(Xf[gi].unsqueeze(-1)).squeeze(-1).squeeze(-1) +
                  head(Xr[gi].unsqueeze(-1)).squeeze(-1).squeeze(-1)) / 2
            vl = torch.nn.functional.mse_loss(pv, tgt[gi]).item()
        n_ep += 1
        if vl < best - 1e-7:
            best, bad = vl, 0
            best_state = {k: v.detach().clone() for k, v in head.state_dict().items()}
        else:
            bad += 1
            if bad >= HP["patience"]:
                break
        if verbose and ep % 10 == 0:
            R(f"      ep {ep:3d} GM21-internal val MSE {vl:.5f}")
    head.load_state_dict(best_state)
    return head, dict(best_val_mse=float(best), epochs_run=n_ep)


def auc_block(score, measured, tag):
    pos = measured >= LOG2FC_THRESH
    neg = measured <= -LOG2FC_THRESH
    keep = pos | neg
    n_tot, n_keep = len(measured), int(keep.sum())
    out = dict(tag=tag, threshold=LOG2FC_THRESH, n_regions_total=n_tot, n_evaluable=n_keep,
               n_ambiguous=int(n_tot - n_keep),
               frac_evaluable=float(n_keep / n_tot) if n_tot else 0.0,
               n_positive=int(pos.sum()), n_negative=int(neg.sum()))
    if n_keep >= 10 and pos.sum() and neg.sum():
        out["auc"] = float(roc_auc_score(pos[keep].astype(int), score[keep]))
    else:
        out["auc"] = None
    return out


def evaluate(head, idx, tag):
    pred = predict(head, idx)
    block = {"tag": tag, "n": len(idx)}
    for nm, meas in [("IMR90_SEN_vs_PRO", IMR_PRO[idx]), ("IMR90_SEN_vs_QUI", IMR_QUI[idx]),
                     ("GM21_RAS_vs_EV", GM21[idx])]:
        sp = stats.spearmanr(pred, meas).statistic
        pe = stats.pearsonr(pred, meas).statistic
        a = auc_block(pred, meas, f"{tag}:{nm}")
        block[nm] = dict(spearman=float(sp), pearson=float(pe),
                         frac_of_ceiling=float(sp / CEILING), **a)
    block["pred_sd"] = float(pred.std())
    block["measured_sd_IMR90"] = float(IMR_PRO[idx].std())
    block["collapse_ratio"] = float(pred.std() / IMR_PRO[idx].std())
    return block, pred


# =====================================================================================
R("\n" + "=" * 78 + "\nTRAINING (GM21 target; IMR90 untouched)\n" + "=" * 78)
t0 = time.time()
head, tlog = train_head(GM21, seed=HP["seed"], verbose=True)
R(f"  trained: GM21-internal val MSE {tlog['best_val_mse']:.5f}, "
  f"{tlog['epochs_run']} epochs, {time.time()-t0:.0f}s")
res["training"] = tlog
torch.save(head.state_dict(), OUT / "head_ois.pt")

R("\n--- PRIMARY: held-out IMR90 (n=2 vs 2), all non-chr9/chr6 regions ---")
blk_rest, pred_rest = evaluate(head, rest_idx, "IMR90_rest")
res["eval_rest"] = blk_rest
for nm in ["IMR90_SEN_vs_PRO", "IMR90_SEN_vs_QUI"]:
    b = blk_rest[nm]
    R(f"  {nm:18s} Spearman {b['spearman']:+.4f} ({100*b['frac_of_ceiling']:.1f}% of the "
      f"{CEILING:.3f} ceiling)  Pearson {b['pearson']:+.4f}")
    R(f"  {'':18s} AUC {b['auc']:.4f}  evaluable {b['n_evaluable']:,}/{b['n_regions_total']:,} "
      f"({100*b['frac_evaluable']:.1f}%)  ambiguous {b['n_ambiguous']:,}  "
      f"pos {b['n_positive']:,} neg {b['n_negative']:,}")
R(f"  in-sample check GM21 Spearman {blk_rest['GM21_RAS_vs_EV']['spearman']:+.4f}")

R("\n--- chr9 + chr6 (held out from everything) ---")
blk_test, pred_test = evaluate(head, te_idx, "IMR90_chr9_chr6")
res["eval_chr9_chr6"] = blk_test
for nm in ["IMR90_SEN_vs_PRO", "IMR90_SEN_vs_QUI"]:
    b = blk_test[nm]
    R(f"  {nm:18s} Spearman {b['spearman']:+.4f} ({100*b['frac_of_ceiling']:.1f}% of ceiling)"
      f"  Pearson {b['pearson']:+.4f}  AUC {b['auc']:.4f}  evaluable {b['n_evaluable']:,}"
      f"/{b['n_regions_total']:,}  ambiguous {b['n_ambiguous']:,}")

R("\n--- COLLAPSE CHECK ---")
R(f"  prediction sd {blk_rest['pred_sd']:.4f} vs measured IMR90 sd "
  f"{blk_rest['measured_sd_IMR90']:.4f}  ratio {blk_rest['collapse_ratio']:.3f}")
R(f"  (a collapsed head would give ratio ~0; the target is signed so sd is the right check)")

# ---- CDKN2A / CDKN1A distal ----
R("\n--- CDKN2A / CDKN1A, promoter vs distal (chr9/chr6, fully held out) ---")
LOCI = {"CDKN2A": ("chr9", 21_967_752, 21_995_301), "CDKN1A": ("chr6", 36_676_461, 36_687_339)}
allpred = np.full(n_reg, np.nan)
allpred[rest_idx] = pred_rest
allpred[te_idx] = pred_test
loci_rows = []
for g, (c, s, e) in LOCI.items():
    m = (chrom == c) & (var["region_center_hg38"].to_numpy() >= s - 200_000) & \
        (var["region_center_hg38"].to_numpy() <= e + 200_000) & cached
    d = var[m]
    dist = np.where(d["region_center_hg38"] < s, s - d["region_center_hg38"],
                    np.where(d["region_center_hg38"] > e, d["region_center_hg38"] - e, 0))
    prox = dist <= 2000
    R(f"\n  {g}: {int(m.sum())} regions (promoter {int(prox.sum())}, distal {int((~prox).sum())})")
    for lab, mm in [("promoter", prox), ("distal", ~prox)]:
        if mm.sum() == 0:
            R(f"    {lab:9s} none"); continue
        idxs = np.where(m)[0][mm]
        row = dict(locus=g, region_class=lab, n=int(mm.sum()),
                   predicted=float(np.nanmean(allpred[idxs])),
                   measured_IMR90=float(np.nanmean(IMR_PRO[idxs])),
                   measured_GM21=float(np.nanmean(GM21[idxs])))
        row["direction_correct"] = bool(np.sign(row["predicted"]) == np.sign(row["measured_IMR90"]))
        loci_rows.append(row)
        R(f"    {lab:9s} predicted {row['predicted']:+.4f} | measured IMR90 "
          f"{row['measured_IMR90']:+.4f} | measured GM21 {row['measured_GM21']:+.4f} -> "
          f"{'CORRECT' if row['direction_correct'] else 'WRONG'}")
res["loci"] = loci_rows

# =====================================================================================
R("\n" + "=" * 78 + "\nCONTROLS\n" + "=" * 78)

R(f"\n--- CONTROL 1: shuffled-label permutations (n={N_PERM}) ---")
real_sp = blk_rest["IMR90_SEN_vs_PRO"]["spearman"]
real_auc = blk_rest["IMR90_SEN_vs_PRO"]["auc"]
rng = np.random.default_rng(20260826)
t0 = time.time()
sps, aucs = [], []
for k in range(N_PERM):
    sh = GM21.copy()
    sh[tr_idx] = sh[rng.permutation(tr_idx)]
    sh[es_idx] = sh[rng.permutation(es_idx)]
    h2, _ = train_head(sh, seed=HP["seed"], max_epochs=40)
    pr = predict(h2, rest_idx)
    sps.append(stats.spearmanr(pr, IMR_PRO[rest_idx]).statistic)
    a = auc_block(pr, IMR_PRO[rest_idx], f"sh{k}")
    if a["auc"] is not None:
        aucs.append(a["auc"])
    del h2
    if (k + 1) % 25 == 0:
        R(f"      {k+1}/{N_PERM}  ({time.time()-t0:.0f}s elapsed)")
sps, aucs = np.array(sps), np.array(aucs)
nb_sp = int((sps >= real_sp).sum()); p_sp = (nb_sp + 1) / (len(sps) + 1)
nb_au = int((aucs >= real_auc).sum()); p_au = (nb_au + 1) / (len(aucs) + 1)
res["control1_shuffled_labels"] = dict(
    n_permutations=len(sps),
    spearman=dict(real=float(real_sp), mean=float(sps.mean()), sd=float(sps.std(ddof=1)),
                  min=float(sps.min()), max=float(sps.max()), n_beating=nb_sp,
                  empirical_p=float(p_sp), passes=bool(p_sp < 0.05)),
    auc=dict(real=float(real_auc), mean=float(aucs.mean()), sd=float(aucs.std(ddof=1)),
             min=float(aucs.min()), max=float(aucs.max()), n_beating=nb_au,
             empirical_p=float(p_au), passes=bool(p_au < 0.05)),
    resolution=f"min attainable p = 1/{len(sps)+1} = {1/(len(sps)+1):.5f}",
    all_spearman=[float(x) for x in sps])
R(f"    Spearman real {real_sp:+.4f} | shuffled {sps.mean():+.4f} +/- {sps.std(ddof=1):.4f} "
  f"range [{sps.min():+.4f}, {sps.max():+.4f}] | beat {nb_sp}/{len(sps)} | p={p_sp:.5f} "
  f"{'PASS' if p_sp<0.05 else 'FAIL'}")
R(f"    AUC      real {real_auc:.4f} | shuffled {aucs.mean():.4f} +/- {aucs.std(ddof=1):.4f} "
  f"range [{aucs.min():.4f}, {aucs.max():.4f}] | beat {nb_au}/{len(aucs)} | p={p_au:.5f} "
  f"{'PASS' if p_au<0.05 else 'FAIL'}")
R(f"    resolution: min attainable p = 1/{len(sps)+1} = {1/(len(sps)+1):.5f}")

R("\n--- CONTROL 2: untrained randomly-initialised head ---")
un = make_head(321)
pu = predict(un, rest_idx)
b_un = dict(spearman=float(stats.spearmanr(pu, IMR_PRO[rest_idx]).statistic),
            **auc_block(pu, IMR_PRO[rest_idx], "untrained"))
res["control2_untrained"] = b_un
R(f"    Spearman {b_un['spearman']:+.4f}  AUC {b_un['auc']:.4f}  "
  f"evaluable {b_un['n_evaluable']:,}")

R("\n--- CONTROL 3: shuffled region-to-embedding assignment (does it use sequence?) ---")
ep = rng.permutation(n_reg)
Xf_s, Xr_s = Xf[ep].clone(), Xr[ep].clone()
Xf_o, Xr_o = Xf, Xr
Xf, Xr = Xf_s, Xr_s
h3, _ = train_head(GM21, seed=HP["seed"])
p3 = predict(h3, rest_idx)
Xf, Xr = Xf_o, Xr_o
b_sh = dict(spearman=float(stats.spearmanr(p3, IMR_PRO[rest_idx]).statistic),
            **auc_block(p3, IMR_PRO[rest_idx], "shuffled_embeddings"))
res["control3_shuffled_embeddings"] = b_sh
R(f"    Spearman {b_sh['spearman']:+.4f}  AUC {b_sh['auc']:.4f}")
R(f"    real {real_sp:+.4f} vs shuffled-embedding {b_sh['spearman']:+.4f} -- a large gap "
  f"means the model is using sequence, not per-region marginals")
del Xf_s, Xr_s

gate = bool(p_sp < 0.05 and p_au < 0.05)
res["control_gate_passed"] = gate
R("\n" + "=" * 78)
R(f"CONTROL GATE: {'PASSED -- generation probe may run' if gate else 'FAILED -- STOP'}")
R("=" * 78)

np.save(OUT / "pred_rest.npy", pred_rest)
np.save(OUT / "pred_test.npy", pred_test)
np.save(OUT / "idx_rest.npy", rest_idx)
np.save(OUT / "idx_test.npy", te_idx)
pd.DataFrame(loci_rows).to_csv(OUT / "locus_checks.csv", index=False)
json.dump(res, open(OUT / "results_train.json", "w"), indent=2)
(OUT / "TRAIN_LOG.txt").write_text("\n".join(LOG))
R("TRAIN_DONE")
