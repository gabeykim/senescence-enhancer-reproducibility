#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 95_gate2_finish.py
#
# Recovers the H2 summary and adds the H1b interaction test (94 died in a print statement). CPU.
#
# CONSUMES: results/gate2_h1b_ablation.csv, gate2_h2_variants.csv, gate2_h1c_contexts.csv
# PRODUCES: results/gate2_ap1_results.json, GATE2_AP1_FINISH_LOG.txt
# ---------------------------------------------------------------------------
"""
GATE 2 -- recover the H2 summary (the run died in a print statement, after every
CSV had been written: `r.sem` resolved to the pandas Series.sem METHOD rather than
the column of that name), and add the H1b stratum-interaction test, which is the
statistic that actually decides H1.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/workspace/ois_gates/scripts")
from gate_common import OUT

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


outdir = Path(OUT) / "results"
res = json.load(open(outdir / "gate2_ap1_partial.json")) if (
    outdir / "gate2_ap1_partial.json").exists() else {}

# ------------------------------------------------------------------ H1b interaction
R("=" * 104)
R("H1b -- DOES THE SIGN OF THE AP-1 ABLATION EFFECT DEPEND ON THE CONTEXT'S RESPONSE?")
R("  This is the test that decides H1. Per region the AP-1-specific effect is")
R("  (ablate AP-1) - (ablate a motif-free stretch of the same length), which removes")
R("  any generic 'shuffling 10 bp changes the prediction' component. Those per-region")
R("  values are then compared between high- and low-response regions.")
R("=" * 104)
ab = pd.read_csv(outdir / "gate2_h1b_ablation.csv")
h1b = {}
for motif in ab["motif"].unique():
    m = ab[ab["motif"] == motif]
    piv = m.pivot_table(index=["region_index", "stratum", "gm21", "imr90", "atac"],
                        columns="arm", values="effect").reset_index()
    if "position_control" not in piv:
        continue
    piv["specific"] = piv["ablated"] - piv["position_control"]
    R(f"\n  {motif}")
    R(f"    {'stratum':8s} {'n':>4s} {'ablated':>10s} {'poscontrol':>11s} "
      f"{'specific':>10s} {'+-sem':>8s} {'p vs 0':>9s}")
    byst = {}
    for st, g in piv.groupby("stratum"):
        v = g["specific"].to_numpy()
        t = stats.ttest_1samp(v, 0.0)
        R(f"    {st:8s} {len(v):4d} {g['ablated'].mean():+10.4f} "
          f"{g['position_control'].mean():+11.4f} {v.mean():+10.4f} "
          f"{v.std(ddof=1)/np.sqrt(len(v)):8.4f} {t.pvalue:9.4f}")
        byst[st] = dict(n=int(len(v)), ablated=float(g["ablated"].mean()),
                        position_control=float(g["position_control"].mean()),
                        specific=float(v.mean()),
                        sem=float(v.std(ddof=1) / np.sqrt(len(v))),
                        p=float(t.pvalue))
    h1b[motif] = dict(by_stratum=byst)
    if {"high", "low"} <= set(piv["stratum"]):
        hi = piv[piv.stratum == "high"]["specific"].to_numpy()
        lo = piv[piv.stratum == "low"]["specific"].to_numpy()
        t = stats.ttest_ind(hi, lo, equal_var=False)
        u = stats.mannwhitneyu(hi, lo)
        R(f"    INTERACTION high vs low: {hi.mean():+.4f} vs {lo.mean():+.4f}, "
          f"difference {hi.mean()-lo.mean():+.4f}, Welch t={t.statistic:+.3f} "
          f"p={t.pvalue:.5f}, Mann-Whitney p={u.pvalue:.5f}")
        h1b[motif]["interaction"] = dict(high=float(hi.mean()), low=float(lo.mean()),
                                         difference=float(hi.mean() - lo.mean()),
                                         welch_t=float(t.statistic),
                                         welch_p=float(t.pvalue),
                                         mannwhitney_p=float(u.pvalue),
                                         n_high=int(len(hi)), n_low=int(len(lo)))
        # does the effect track the region's measured response continuously?
        sp = stats.spearmanr(piv["gm21"], piv["specific"])
        R(f"    continuous: Spearman(region GM21 H3K27ac response, AP-1-specific "
          f"ablation effect) = {sp.statistic:+.4f}, p={sp.pvalue:.5f}, n={len(piv)}")
        h1b[motif]["continuous_vs_gm21"] = dict(rho=float(sp.statistic),
                                                p=float(sp.pvalue), n=int(len(piv)))
res["H1b_interaction"] = h1b

# ------------------------------------------------------------------ H2 panel
R("\n" + "=" * 104)
R("H2 -- AP-1 VARIANT PANEL at k=8 (adjusted against the scrambled TRE7 control)")
R("=" * 104)
h2 = pd.read_csv(outdir / "gate2_h2_variants.csv").sort_values("ap1_best_rel",
                                                               ascending=False)
R(f"  {'variant':10s} {'seq':10s} {'AP1 rel':>8s} {'adj delta':>10s} {'+-sem':>7s} "
  f"{'p':>8s}  {'created NFKB/CEBPB/ETS':>22s}")
for _, r in h2.iterrows():
    R(f"  {r['variant']:10s} {r['seq']:10s} {r['ap1_best_rel']:8.3f} "
      f"{r['adjusted_delta']:+10.4f} {r['sem']:7.4f} {r['p']:8.4f}  "
      f"{r['created_NFKB_RELA']:6.1f}/{r['created_CEBPB']:.1f}/{r['created_ETS1']:.1f}")

ok = h2.dropna(subset=["ap1_best_rel"])
sp = stats.spearmanr(ok["ap1_best_rel"], ok["adjusted_delta"])
pe = stats.pearsonr(ok["ap1_best_rel"], ok["adjusted_delta"])
R(f"\n  affinity vs effect across n={len(ok)} variants:")
R(f"    Spearman rho = {sp.statistic:+.4f} (p={sp.pvalue:.5f})")
R(f"    Pearson    r = {pe.statistic:+.4f} (p={pe.pvalue:.5f})")
R("  A strong NEGATIVE correlation means a better AP-1 match gives a more negative")
R("  effect, i.e. the model is reading AP-1 affinity rather than reacting to an")
R("  incidental 7-mer. A near-zero or erratic relation would indicate an artifact.")

# single-base mutants only: the sharpest form of the H2 test
mut = ok[ok["variant"].str.startswith("m")]
spm = stats.spearmanr(mut["ap1_best_rel"], mut["adjusted_delta"])
R(f"\n  restricted to the {len(mut)} single-base mutants of TGACTCA: "
  f"Spearman rho = {spm.statistic:+.4f} (p={spm.pvalue:.5f})")
tre = ok[ok["variant"] == "TRE7"].iloc[0]
cre = ok[ok["variant"] == "CRE8"].iloc[0]
R(f"  TRE heptamer TGACTCA : adjusted {tre['adjusted_delta']:+.4f} "
  f"+- {tre['sem']:.4f} (AP-1 rel {tre['ap1_best_rel']:.3f})")
R(f"  CRE octamer TGACGTCA : adjusted {cre['adjusted_delta']:+.4f} "
  f"+- {cre['sem']:.4f} (AP-1 rel {cre['ap1_best_rel']:.3f})")
created = h2[["created_NFKB_RELA", "created_CEBPB", "created_ETS1"]]
R(f"\n  unintended sites created by the insertion (mean per cassette, rel>=0.90): "
  f"NF-kB {created['created_NFKB_RELA'].mean():.2f}, "
  f"C/EBPb {created['created_CEBPB'].mean():.2f}, "
  f"ETS {created['created_ETS1'].mean():.2f}; "
  f"range across variants "
  f"{created.min().min():.1f}-{created.max().max():.1f}")
res["H2"] = dict(n_variants=int(len(h2)),
                 affinity_vs_effect_spearman=float(sp.statistic),
                 affinity_vs_effect_spearman_p=float(sp.pvalue),
                 affinity_vs_effect_pearson=float(pe.statistic),
                 affinity_vs_effect_pearson_p=float(pe.pvalue),
                 mutants_only_spearman=float(spm.statistic),
                 mutants_only_spearman_p=float(spm.pvalue),
                 tre7=float(tre["adjusted_delta"]), cre8=float(cre["adjusted_delta"]),
                 variants=h2.to_dict("records"))

res["H1c"] = pd.read_csv(outdir / "gate2_h1c_contexts.csv").to_dict("records")
res["H1a"] = pd.read_csv(outdir / "gate2_h1a_ism_enrichment.csv").to_dict("records")
json.dump(res, open(outdir / "gate2_ap1_results.json", "w"), indent=2, default=str)
Path(outdir, "GATE2_AP1_FINISH_LOG.txt").write_text("\n".join(LOG))
R("\nGATE2_FINISH_DONE")
