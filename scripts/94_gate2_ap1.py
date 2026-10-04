#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 94_gate2_ap1.py
#
# Gate 2 H1/H2: ISM enrichment, native-motif ablation, multi-context insertion, variant panel. GPU STEP.
#
# CONSUMES: gate2_motif_calls.csv, ism_importance.npz, via 90_gate_common.py
# PRODUCES: results/gate2_h1a_ism_enrichment.csv, gate2_h1b_ablation.csv, gate2_h1c_{contexts,h2_raw}.csv, gate2_h2_variants.csv
# ---------------------------------------------------------------------------
"""
GATE 2 -- RESOLVE THE AP-1 ANOMALY.

Prior probe: inserting the FOS::JUN consensus into random 200 bp background LOWERS
predicted senescence H3K27ac (delta -0.464, monotonic, rho -0.375, against a flat
scrambled control of -0.016). Three hypotheses are separated here.

H1  CONTEXT DEPENDENCE -- AP-1 in real regions carries positive model attribution;
    the negative insertion result is a property of the random background.
      H1a  the already-saved ISM importance arrays: is attribution elevated at AP-1
           positions in real high-response regions? NOTE: those arrays store
           mean |delta| over the three alternative bases, so they are UNSIGNED and
           can show only that the model ATTENDS to a position, never the direction.
      H1b  the direct signed test, and the one that actually decides H1: ablate the
           native AP-1 match in real regions and read the SIGNED change. If AP-1
           contributes positively in context, ablating it must LOWER the prediction.
           Ablation = mononucleotide shuffle of the matched bases, so length and
           composition are preserved. Two controls: an equally-sized shuffle at a
           motif-free position in the same window, and NF-kB ablation in regions
           carrying a strong NF-kB site (which should lower the prediction).
      H1c  is the negative insertion result specific to the one background used
           before? The dose response is repeated in several independent neutral
           genomic contexts.

H2  INSERTION ARTIFACT -- the effect comes from the particular 7-mer, not from AP-1
    recognition. Tested with a variant panel: the TRE heptamer, the CRE octamer,
    and all 21 single-base mutants of TGACTCA. If the response tracks AP-1 PWM
    affinity it is AP-1 recognition; if it jumps around independently of affinity
    it is an artifact. Cassettes are also rescanned for unintended sites of the
    other three motifs created by the insertion.

H3 is handled separately in gate2_h3.py (no GPU needed).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/workspace/ois_gates/scripts")
from gate_common import (BG, CASSETTE, MOTIF_SEQ, OUT, RUN, SEQ_LEN, build_pwms,
                         load_model, load_sizes, make_scorer, pick_batch,
                         scan_pwm_fast, scramble, window_on)

ISM_W = 500
CORE = 300                       # the scan window used by gate2_scan.py
N_SEEDS = 10
K_MAIN = 8                       # the dose the prior probe reported deltas at
K_DOSE = [2, 4, 6, 8]
AP1 = "TGACTCA"
REQUIRED_DELTA = 0.8612759113311768

N_AP1_ABLATE = 24                # per response stratum
N_NFKB_ABLATE = 24
N_CTX = 4

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    t0 = time.time()
    outdir = Path(OUT) / "results"
    calls = pd.read_csv(outdir / "gate2_motif_calls.csv")
    pwms = build_pwms()
    head, core = load_model()
    batch = pick_batch(core)
    score = make_scorer(head, core, batch)
    sizes = load_sizes()
    R(f"loaded: {len(calls)} regions, trunk batch {batch}")

    res = {}

    # =================================================================== H1a
    R("\n" + "=" * 100)
    R("H1a -- SAVED ISM IMPORTANCE AT AP-1 POSITIONS (unsigned; direction comes from H1b)")
    R("=" * 100)
    z = np.load(f"{RUN}/results/ism_importance.npz")
    ismtab = pd.read_csv(f"{RUN}/results/ism_per_region.csv")
    anchors = dict(zip(ismtab["region"], zip(ismtab["chrom"], ismtab["anchor"])))
    h1a = []
    for key in z.files:
        imp = z[key]
        chrom, anchor = anchors[key]
        ref, ws, off = window_on(chrom, int(anchor), sizes)
        ctr = int(anchor) - ws
        win = ref[ctr - ISM_W // 2: ctr + ISM_W // 2]
        assert len(win) == ISM_W == len(imp)
        row = dict(region=key, measured=float(
            ismtab.loc[ismtab.region == key, "measured_IMR90"].iloc[0]))
        for nm in pwms:
            mask = np.zeros(ISM_W, bool)
            for st, rel, strand in scan_pwm_fast(win, pwms[nm], 0.85):
                mask[st:st + pwms[nm]["length"]] = True
            if mask.any() and (~mask).any():
                row[f"{nm}_enr"] = float(imp[mask].mean() / imp[~mask].mean())
                row[f"{nm}_bases"] = int(mask.sum())
            else:
                row[f"{nm}_enr"] = np.nan
                row[f"{nm}_bases"] = int(mask.sum())
        h1a.append(row)
    h1a = pd.DataFrame(h1a)
    h1a.to_csv(outdir / "gate2_h1a_ism_enrichment.csv", index=False)
    R(f"  n = {len(h1a)} high-response held-out regions (rel>=0.85 PWM matches)")
    R(f"  {'motif':12s} {'n_regions':>10s} {'mean enr':>9s} {'sd':>7s} {'p vs 1.0':>10s}")
    h1a_sum = {}
    for nm in pwms:
        v = h1a[f"{nm}_enr"].dropna()
        if len(v) > 1:
            t = stats.ttest_1samp(v, 1.0)
            R(f"  {nm:12s} {len(v):10d} {v.mean():9.3f} {v.std(ddof=1):7.3f} "
              f"{t.pvalue:10.4f}")
            h1a_sum[nm] = dict(n=int(len(v)), mean_enrichment=float(v.mean()),
                               sd=float(v.std(ddof=1)), p_vs_1=float(t.pvalue))
    res["H1a"] = dict(note="unsigned importance; cannot give direction of effect",
                      n_regions=int(len(h1a)), per_motif=h1a_sum)

    # =================================================================== H1b
    R("\n" + "=" * 100)
    R("H1b -- SIGNED ABLATION OF NATIVE MOTIFS IN REAL REGIONS")
    R("  effect = pred(ablated) - pred(wild type).  NEGATIVE => the motif contributes")
    R("  POSITIVELY to predicted senescence H3K27ac in that real context.")
    R("=" * 100)

    def pick(motif, n_each, strata=True, rel=0.95, near=100):
        d = calls[(calls[f"{motif}_best_rel"] >= rel)
                  & (np.abs(calls[f"{motif}_best_pos"] + pwms[motif]["length"] / 2
                            - CORE / 2) <= near)].copy()
        d = d[np.isfinite(d["GM21_RAS_vs_EV"])]
        d = d.sort_values("GM21_RAS_vs_EV")
        if not strata:
            return d.tail(n_each).assign(stratum="high")
        assert len(d) > 2 * n_each, f"{motif}: only {len(d)} candidates for 2x{n_each}"
        return pd.concat([d.head(n_each).assign(stratum="low"),
                          d.tail(n_each).assign(stratum="high")])

    sel_ap1 = pick("AP1_FOSJUN", N_AP1_ABLATE)
    sel_nfkb = pick("NFKB_RELA", N_NFKB_ABLATE, strata=False)
    R(f"  AP-1 ablation set : n={len(sel_ap1)} "
      f"({(sel_ap1.stratum=='high').sum()} high / {(sel_ap1.stratum=='low').sum()} low "
      f"GM21 response), best_rel >= 0.95 within +-100 bp of anchor")
    R(f"  NF-kB control set : n={len(sel_nfkb)} (high GM21 response)")

    jobs, meta = [], []
    for motif, sel in (("AP1_FOSJUN", sel_ap1), ("NFKB_RELA", sel_nfkb)):
        L = pwms[motif]["length"]
        for _, r in sel.iterrows():
            ref, ws, off = window_on(r["chrom"], int(r["anchor"]), sizes)
            ctr = int(r["anchor"]) - ws
            c0 = ctr - CORE // 2
            mst = c0 + int(r[f"{motif}_best_pos"])
            win = ref[c0:c0 + CORE]
            rng = np.random.default_rng(int(r["region_index"]))
            # motif-free position control, same length, nearest to centre
            occupied = np.zeros(CORE, bool)
            for nm in pwms:
                for st, _, _ in scan_pwm_fast(win, pwms[nm], 0.80):
                    occupied[st:st + pwms[nm]["length"]] = True
            free = [p for p in range(CORE - L)
                    if not occupied[p:p + L].any()]
            base = dict(region_index=int(r["region_index"]), motif=motif,
                        chrom=r["chrom"], anchor=int(r["anchor"]),
                        stratum=r["stratum"], best_rel=float(r[f"{motif}_best_rel"]),
                        gm21=float(r["GM21_RAS_vs_EV"]),
                        imr90=float(r["IMR90_SEN_vs_PRO"]),
                        atac=float(r["GM21_ATAC_RAS_vs_EV"]))
            jobs.append(ref); meta.append({**base, "arm": "wt"})
            abl = scramble(ref[mst:mst + L], int(r["region_index"]))
            jobs.append(ref[:mst] + abl + ref[mst + L:])
            meta.append({**base, "arm": "ablated"})
            if free:
                p = c0 + min(free, key=lambda q: abs(q + L / 2 - CORE / 2))
                ctl = scramble(ref[p:p + L], int(r["region_index"]) + 7)
                jobs.append(ref[:p] + ctl + ref[p + L:])
                meta.append({**base, "arm": "position_control"})
    R(f"  {len(jobs)} sequences for H1b")

    preds = []
    for i in range(0, len(jobs), 50):
        preds.append(score(jobs[i:i + 50]))
        R(f"    H1b scored {min(i+50,len(jobs))}/{len(jobs)} ({time.time()-t0:.0f}s)")
    ab = pd.DataFrame(meta); ab["pred"] = np.concatenate(preds)
    wt = ab[ab.arm == "wt"].set_index(["region_index", "motif"])["pred"]
    ab["effect"] = ab["pred"] - ab.set_index(["region_index", "motif"]).index.map(wt)
    ab.to_csv(outdir / "gate2_h1b_ablation.csv", index=False)

    R(f"\n  {'motif':12s} {'arm':18s} {'stratum':8s} {'n':>4s} {'mean effect':>12s} "
      f"{'+-sem':>8s} {'p vs 0':>9s}")
    h1b_sum = []
    for (motif, arm, st), g in ab[ab.arm != "wt"].groupby(["motif", "arm", "stratum"]):
        v = g["effect"].to_numpy()
        t = stats.ttest_1samp(v, 0.0)
        R(f"  {motif:12s} {arm:18s} {st:8s} {len(v):4d} {v.mean():+12.4f} "
          f"{v.std(ddof=1)/np.sqrt(len(v)):8.4f} {t.pvalue:9.4f}")
        h1b_sum.append(dict(motif=motif, arm=arm, stratum=st, n=int(len(v)),
                            mean_effect=float(v.mean()),
                            sem=float(v.std(ddof=1) / np.sqrt(len(v))),
                            p=float(t.pvalue)))
    # ablated vs position control, paired within region
    for motif in ("AP1_FOSJUN", "NFKB_RELA"):
        a = ab[(ab.motif == motif) & (ab.arm == "ablated")].set_index("region_index")["effect"]
        c = ab[(ab.motif == motif) & (ab.arm == "position_control")].set_index("region_index")["effect"]
        common = a.index.intersection(c.index)
        if len(common) > 1:
            t = stats.ttest_rel(a[common], c[common])
            R(f"  {motif}: ablated vs position control, paired n={len(common)}: "
              f"{a[common].mean():+.4f} vs {c[common].mean():+.4f}, "
              f"diff {(a[common]-c[common]).mean():+.4f}, p={t.pvalue:.4f}")
            h1b_sum.append(dict(motif=motif, arm="ablated_vs_position_control",
                                stratum="all", n=int(len(common)),
                                mean_effect=float((a[common] - c[common]).mean()),
                                sem=float((a[common] - c[common]).std(ddof=1) / np.sqrt(len(common))),
                                p=float(t.pvalue)))
    res["H1b"] = h1b_sum

    # =================================================================== H1c + H2
    R("\n" + "=" * 100)
    R("H1c -- AP-1 INSERTION DOSE RESPONSE ACROSS INDEPENDENT NEUTRAL CONTEXTS")
    R("=" * 100)
    neutral = calls[(np.abs(calls["GM21_RAS_vs_EV"]) < 0.05)
                    & (np.abs(calls["IMR90_SEN_vs_PRO"]) < 0.05)
                    & (calls["chrom"] != "chr1")].copy()
    neutral = neutral.drop_duplicates("chrom")
    ctxs = [("chr1", 265888)]
    for _, r in neutral.iterrows():
        if len(ctxs) >= N_CTX:
            break
        ch, pos = r["chrom"], int(r["anchor"])
        if ch not in sizes or pos - SEQ_LEN // 2 < 0 or pos + SEQ_LEN // 2 > sizes[ch]:
            continue            # window would be clipped; anchor must stay centred
        ctxs.append((ch, pos))
    R(f"  contexts ({len(ctxs)}, each on a different chromosome, window unclipped): "
      f"{', '.join(f'{c}:{p:,}' for c, p in ctxs)}")

    rng = np.random.default_rng(20261004)
    bgs = ["".join(rng.choice(list("ACGT"), size=CASSETTE,
                              p=[BG["A"], BG["C"], BG["G"], BG["T"]]))
           for _ in range(N_SEEDS)]

    def place(bg, motif, k, variant, seed):
        cas = list(bg)
        if k == 0:
            return "".join(cas)
        span = CASSETTE // k
        for j in range(k):
            use = motif if variant == "real" else scramble(motif, seed * 1009 + j * 31)
            L = len(use)
            st = max(0, min(CASSETTE - L, j * span + (span - L) // 2))
            cas[st:st + L] = list(use)
        return "".join(cas)

    jobs, meta = [], []
    ctx_seq = {}
    for ci, (ch, pos) in enumerate(ctxs):
        s, ws, off = window_on(ch, pos, sizes)
        assert off == SEQ_LEN // 2, f"context {ch}:{pos} window clipped"
        ctx_seq[ci] = s
        lo = SEQ_LEN // 2 - CASSETTE // 2
        for seed in range(N_SEEDS):
            bg = bgs[seed]
            jobs.append(s[:lo] + bg + s[lo + CASSETTE:])
            meta.append(dict(part="H1c", ctx=f"{ch}:{pos}", ctx_i=ci, seed=seed,
                             variant="real", k=0, name="empty"))
            for k in (4, K_MAIN):
                for variant in ("real", "scrambled"):
                    jobs.append(s[:lo] + place(bg, AP1, k, variant, seed)
                                + s[lo + CASSETTE:])
                    meta.append(dict(part="H1c", ctx=f"{ch}:{pos}", ctx_i=ci,
                                     seed=seed, variant=variant, k=k, name="TRE7"))

    R("\n" + "=" * 100)
    R("H2 -- AP-1 VARIANT PANEL (context 0, the background used by the prior probe)")
    R("=" * 100)
    variants = {"TRE7": AP1, "CRE8": "TGACGTCA"}
    for p in range(len(AP1)):
        for b in "ACGT":
            if b != AP1[p]:
                variants[f"m{p+1}{AP1[p]}{b}"] = AP1[:p] + b + AP1[p + 1:]
    R(f"  {len(variants)} variants (TRE heptamer, CRE octamer, "
      f"{len(variants)-2} single-base mutants) + scrambled control")
    s0 = ctx_seq[0]
    lo = SEQ_LEN // 2 - CASSETTE // 2
    for seed in range(N_SEEDS):
        bg = bgs[seed]
        jobs.append(s0[:lo] + bg + s0[lo + CASSETTE:])
        meta.append(dict(part="H2", ctx=f"{ctxs[0][0]}:{ctxs[0][1]}", ctx_i=0,
                         seed=seed, variant="real", k=0, name="empty"))
        for nm, v in variants.items():
            jobs.append(s0[:lo] + place(bg, v, K_MAIN, "real", seed) + s0[lo + CASSETTE:])
            meta.append(dict(part="H2", ctx=f"{ctxs[0][0]}:{ctxs[0][1]}", ctx_i=0,
                             seed=seed, variant="real", k=K_MAIN, name=nm))
        jobs.append(s0[:lo] + place(bg, AP1, K_MAIN, "scrambled", seed) + s0[lo + CASSETTE:])
        meta.append(dict(part="H2", ctx=f"{ctxs[0][0]}:{ctxs[0][1]}", ctx_i=0, seed=seed,
                         variant="scrambled", k=K_MAIN, name="TRE7"))
        for k in (2, 6):
            for nm in ("TRE7", "CRE8"):
                jobs.append(s0[:lo] + place(bg, variants[nm], k, "real", seed)
                            + s0[lo + CASSETTE:])
                meta.append(dict(part="H2", ctx=f"{ctxs[0][0]}:{ctxs[0][1]}", ctx_i=0,
                                 seed=seed, variant="real", k=k, name=nm))

    R(f"\n  {len(jobs)} sequences for H1c + H2")
    preds = []
    for i in range(0, len(jobs), 50):
        preds.append(score(jobs[i:i + 50]))
        if (i // 50) % 4 == 0:
            el = time.time() - t0
            R(f"    scored {min(i+50,len(jobs))}/{len(jobs)} ({el:.0f}s)")
    ins = pd.DataFrame(meta); ins["pred"] = np.concatenate(preds)
    basemap = ins[ins["name"] == "empty"].set_index(["part", "ctx_i", "seed"])["pred"]
    ins["delta"] = ins["pred"] - ins.set_index(["part", "ctx_i", "seed"]).index.map(basemap)
    ins.to_csv(outdir / "gate2_h1c_h2_raw.csv", index=False)

    R(f"\n  H1c: AP-1 (TRE7) dose response by context, adjusted = real - scrambled")
    R(f"  {'context':20s} {'k':>3s} {'raw':>9s} {'scram':>9s} {'adj':>9s} {'+-sem':>7s} {'p':>9s}")
    h1c = []
    for (ctx, k), g in ins[(ins.part == "H1c") & (ins.k > 0)].groupby(["ctx", "k"]):
        a = g[g.variant == "real"].sort_values("seed")["delta"].to_numpy()
        b = g[g.variant == "scrambled"].sort_values("seed")["delta"].to_numpy()
        adj = a - b
        t = stats.ttest_1samp(adj, 0.0)
        R(f"  {ctx:20s} {k:3d} {a.mean():+9.4f} {b.mean():+9.4f} {adj.mean():+9.4f} "
          f"{adj.std(ddof=1)/np.sqrt(len(adj)):7.4f} {t.pvalue:9.4f}")
        h1c.append(dict(context=ctx, k=int(k), n_seeds=len(adj), raw=float(a.mean()),
                        scrambled=float(b.mean()), adjusted=float(adj.mean()),
                        sem=float(adj.std(ddof=1) / np.sqrt(len(adj))),
                        p=float(t.pvalue)))
    pd.DataFrame(h1c).to_csv(outdir / "gate2_h1c_contexts.csv", index=False)
    res["H1c"] = h1c

    # ---- H2 summary: delta vs realised PWM affinity ---------------------------
    scr = ins[(ins.part == "H2") & (ins.variant == "scrambled")].set_index("seed")["delta"]
    h2 = []
    for nm, v in variants.items():
        g = ins[(ins.part == "H2") & (ins["name"] == nm) & (ins.k == K_MAIN)
                & (ins.variant == "real")].sort_values("seed")
        adj = g["delta"].to_numpy() - scr.loc[g["seed"]].to_numpy()
        # realised AP-1 affinity of the built cassette, averaged over seeds
        rels, extra = [], {n: 0 for n in pwms if n != "AP1_FOSJUN"}
        for seed in range(N_SEEDS):
            cas = place(bgs[seed], v, K_MAIN, "real", seed)
            hits = scan_pwm_fast(cas, pwms["AP1_FOSJUN"], 0.0)
            rels.append(max(h[1] for h in hits) if hits else np.nan)
            for n in extra:
                extra[n] += len(scan_pwm_fast(cas, pwms[n], 0.90))
        t = stats.ttest_1samp(adj, 0.0)
        h2.append(dict(variant=nm, seq=v, n_seeds=len(adj),
                       adjusted_delta=float(adj.mean()),
                       sem=float(adj.std(ddof=1) / np.sqrt(len(adj))),
                       p=float(t.pvalue), ap1_best_rel=float(np.nanmean(rels)),
                       **{f"created_{n}": extra[n] / N_SEEDS for n in extra}))
    h2 = pd.DataFrame(h2).sort_values("ap1_best_rel", ascending=False)
    h2.to_csv(outdir / "gate2_h2_variants.csv", index=False)
    R(f"\n  H2 variant panel at k={K_MAIN} (adjusted vs scrambled TRE7), "
      f"n={N_SEEDS} seeds each")
    R(f"  {'variant':10s} {'seq':10s} {'AP1 rel':>8s} {'adj delta':>10s} {'+-sem':>7s} "
      f"{'p':>8s} {'created NFKB/CEBPB/ETS':>24s}")
    for _, r in h2.iterrows():
        R(f"  {r.variant:10s} {r.seq:10s} {r.ap1_best_rel:8.3f} {r.adjusted_delta:+10.4f} "
          f"{r.sem:7.4f} {r.p:8.4f}   "
          f"{r.created_NFKB_RELA:.1f}/{r.created_CEBPB:.1f}/{r.created_ETS1:.1f}")
    ok = h2.dropna(subset=["ap1_best_rel"])
    sp = stats.spearmanr(ok["ap1_best_rel"], ok["adjusted_delta"])
    pe = stats.pearsonr(ok["ap1_best_rel"], ok["adjusted_delta"])
    R(f"\n  affinity vs effect across {len(ok)} variants: "
      f"Spearman rho = {sp.statistic:+.3f} (p={sp.pvalue:.4f}), "
      f"Pearson r = {pe.statistic:+.3f} (p={pe.pvalue:.4f})")
    R("  (a strong NEGATIVE correlation means the effect tracks AP-1 affinity, i.e. the")
    R("   model is responding to AP-1 recognition rather than to an incidental 7-mer)")
    res["H2"] = dict(n_variants=int(len(h2)), n_seeds=N_SEEDS,
                     affinity_vs_effect_spearman=float(sp.statistic),
                     affinity_vs_effect_spearman_p=float(sp.pvalue),
                     affinity_vs_effect_pearson=float(pe.statistic),
                     affinity_vs_effect_pearson_p=float(pe.pvalue),
                     variants=h2.to_dict("records"))

    json.dump(res, open(outdir / "gate2_ap1_results.json", "w"), indent=2, default=str)
    Path(outdir, "GATE2_AP1_LOG.txt").write_text("\n".join(LOG))
    R(f"\nGATE2_AP1_DONE in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
