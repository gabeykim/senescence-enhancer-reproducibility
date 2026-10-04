#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 44_verify_expanded_trainset.py
#
# Verify the expanded trainset; cross-study ON-OFF agreement (+0.531/+0.142/+0.093).
#
# CONSUMES: output/trainset_expanded/senescence_trainset_multistudy.h5ad
# PRODUCES: output/trainset_expanded/VERIFY_LOG.txt
# ---------------------------------------------------------------------------
"""
Verification + honest replication accounting for the expanded training set.

The headline number a task count gives you (24 ON tasks) is NOT the number that governs
whether the shuffled-label control can be beaten. What governs that is the number of
INDEPENDENT BIOLOGICAL UNITS -- distinct cultures that could have been sampled
differently. Serial timepoints from one culture, and replicate wells of one experiment,
are not independent draws. This script counts both and reports the gap.
"""
import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
OUT = ROOT / "output" / "trainset_expanded"
LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def biological_unit(row):
    """Map a task to the independent culture it came from.

    GSE175533 : one WI-38 serial-passage lineage. The PDL timepoints are sequential
                samples of THAT SAME culture as it ages, so they are not independent
                draws; the replicate letter (A/B/C) marks parallel culture replicates
                and IS treated as the unit. This is the generous reading -- even these
                descend from a single WI-38 stock at a single starting PDL.
    GSE205692 : 'biol rep N' is explicit in the sample titles; the many day-timepoints
                within a replicate are serial samples of that one culture.
    GSE74324  : the trailing _1/_2/_3 are the biological replicates.
    GSE99028  : two replicate wells (ON) / two population-doubling samples (OFF).
    """
    s, d, n = row["study"], str(row["detail"]), str(row["sample"])
    if s == "GSE175533":
        rep = n.split("_")[-1]
        return f"GSE175533|WI38_lineage|rep{rep}"
    if s == "GSE205692":
        import re
        m = re.search(r"biol_rep=(\d)", d)
        return f"GSE205692|GM21|rep{m.group(1) if m else '?'}"
    if s == "GSE74324":
        return f"GSE74324|IMR90|rep{n.split('_')[-1]}"
    if s == "GSE99028":
        return f"GSE99028|IMR90|{n}"
    return f"{s}|{n}"


def main():
    p = OUT / "senescence_trainset_multistudy.h5ad"
    ad = anndata.read_h5ad(p)
    R("=" * 86)
    R("VERIFICATION -- expanded multi-study senescence training set")
    R("=" * 86)
    R(f"\n  file      : {p}")
    R(f"  shape     : {ad.shape}  (n_tasks, n_intervals)")
    R(f"  X dtype   : {ad.X.dtype}")
    R(f"  layers    : {sorted(k for k in ad.layers.keys() if k)}")
    R(f"  finite    : X {np.isfinite(ad.X).all()}  "
      f"log_native {np.isfinite(ad.layers['log_native']).all()}")

    # ---------- contract ----------
    R("\n[A] CONTRACT CHECKS\n" + "-" * 86)
    w = (ad.var["end"] - ad.var["start"]).unique()
    R(f"  interval widths distinct values : {w}  (expect [524288])")
    assert list(w) == [524288]
    centred = (ad.var["tss"] - ad.var["start"]).unique()
    R(f"  TSS offset within window        : {centred}  (expect [262144])")
    assert list(centred) == [262144]
    R(f"  var columns                     : {list(ad.var.columns)}")
    R(f"  negative coordinates            : {int((ad.var['start'] < 0).sum())} (expect 0)")

    # ---------- leakage ----------
    R("\n[B] LEAKAGE CHECKS\n" + "-" * 86)
    for a, b in [("train", "val"), ("train", "test"), ("val", "test")]:
        ca = set(ad.var.loc[ad.var["split"] == a, "chrom"])
        cb = set(ad.var.loc[ad.var["split"] == b, "chrom"])
        R(f"  {a:5s} vs {b:5s} chromosome overlap: {sorted(ca & cb) or 'none'}")
        assert not (ca & cb)
    R(f"  duplicate gene names            : {int(ad.var.index.duplicated().sum())}")
    R(f"  duplicate task names            : {int(ad.obs.index.duplicated().sum())}")
    R(f"  CDKN2A on {ad.var.loc['CDKN2A','chrom'] if 'CDKN2A' in ad.var.index else 'ABSENT':6s}"
      f" -> split {ad.var.loc['CDKN2A','split'] if 'CDKN2A' in ad.var.index else '-'}")
    R(f"  CDKN1A on {ad.var.loc['CDKN1A','chrom'] if 'CDKN1A' in ad.var.index else 'ABSENT':6s}"
      f" -> split {ad.var.loc['CDKN1A','split'] if 'CDKN1A' in ad.var.index else '-'}")

    # ---------- does the corrected target still carry senescence signal? ----------
    R("\n[C] SIGNAL SURVIVES CORRECTION?\n" + "-" * 86)
    R("  Batch correction can remove the biology along with the batch. Checked directly:")
    obs = ad.obs
    trn = obs["role"] == "train_eligible"
    on = (trn & obs["cls"].eq("senescent")).to_numpy()
    off = (trn & obs["cls"].eq("proliferating")).to_numpy()
    gi = {g: i for i, g in enumerate(ad.var.index)}
    R("\n   marker    corrected ON-OFF (SD units)   native log2FC")
    for g in ["CDKN1A", "IL6", "CXCL8", "LMNB1", "MKI67", "CDKN2A"]:
        if g not in gi:
            continue
        j = gi[g]
        dc = ad.X[on, j].mean() - ad.X[off, j].mean()
        dn = (ad.layers["log_native"][on, j].mean()
              - ad.layers["log_native"][off, j].mean()) / np.log(2)
        R(f"   {g:8s} {dc:+27.3f}   {dn:+13.3f}")

    d = ad.X[on].mean(axis=0) - ad.X[off].mean(axis=0)
    R(f"\n  genome-wide corrected ON-OFF difference: mean {d.mean():+.4f}, "
      f"sd {d.std():.4f}, |d|>1 SD for {int((np.abs(d) > 1).sum()):,}/{len(d):,} genes")

    # per-study agreement of the ON-OFF response -- do the studies agree on direction?
    R("\n  Cross-study agreement of the ON-OFF response (Spearman of per-gene deltas):")
    deltas = {}
    for s in sorted(obs.loc[trn, "study"].unique()):
        m_on = (trn & obs["study"].eq(s) & obs["cls"].eq("senescent")).to_numpy()
        m_off = (trn & obs["study"].eq(s) & obs["cls"].eq("proliferating")).to_numpy()
        if m_on.sum() and m_off.sum():
            deltas[s] = ad.X[m_on].mean(axis=0) - ad.X[m_off].mean(axis=0)
    ks = sorted(deltas)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            r = stats.spearmanr(deltas[ks[i]], deltas[ks[j]]).statistic
            R(f"    {ks[i]:11s} vs {ks[j]:11s} rho = {r:+.3f}")
    R("  A positive rho means the studies agree on which genes move; this is the")
    R("  cross-study consistency the single-study build could not measure at all.")

    # ---------- honest replication accounting ----------
    R("\n[D] INDEPENDENT BIOLOGICAL UNITS (the number that actually matters)\n" + "-" * 86)
    obs = obs.copy()
    obs["bio_unit"] = obs.apply(biological_unit, axis=1)
    trn_obs = obs[obs["role"] == "train_eligible"]
    on_obs = trn_obs[trn_obs["cls"] == "senescent"]
    R("\n  ON arm, tasks vs independent units:")
    tab = on_obs.groupby("study").agg(tasks=("cls", "size"),
                                      independent_units=("bio_unit", "nunique"))
    R(tab.to_string())
    R(f"\n  TOTAL ON tasks              : {len(on_obs)}")
    R(f"  TOTAL independent ON units  : {on_obs['bio_unit'].nunique()}")
    R("\n  Why the two numbers differ:")
    R("    GSE205692 contributes 12 ON tasks but only 3 independent cultures -- the other")
    R("    9 are additional day-timepoints sampled from those same 3 cultures.")
    R("    GSE175533 contributes 9 ON tasks from one WI-38 serial-passage lineage; the 3")
    R("    counted units are parallel replicate cultures of that single lineage.")
    R("    Serial timepoints are NOT independent draws, so they inflate the task count")
    R("    without adding the biological replication the shuffled-label control tests.")

    off_obs = trn_obs[trn_obs["cls"] == "proliferating"]
    R(f"\n  OFF arm: {len(off_obs)} tasks, {off_obs['bio_unit'].nunique()} independent units")
    R("\n  Comparison with the build that failed:")
    R("    BEFORE : ON = 9 tasks, 1 cell line (WI-38), 1 lineage, 1 mechanism")
    R("             -> effectively n~3 independent units, all one lineage")
    R(f"    AFTER  : ON = {len(on_obs)} tasks, {on_obs['cell_line'].nunique()} cell lines, "
      f"{on_obs['bio_unit'].nunique()} independent units, "
      f"{on_obs['mechanism_class'].nunique()} distinct mechanisms")

    # persist bio_unit into the object itself -- it is the quantity that governs whether a
    # shuffled-label control can be beaten, so it belongs with the data, not only in a CSV
    ad.obs["bio_unit"] = obs["bio_unit"]
    ad.write_h5ad(p)
    R(f"\n  wrote bio_unit back into {p.name}; obs columns now: {list(ad.obs.columns)}")

    obs.to_csv(OUT / "multistudy_task_table.csv")
    (OUT / "VERIFY_LOG.txt").write_text("\n".join(LOG))
    R(f"  wrote {OUT/'VERIFY_LOG.txt'}")


if __name__ == "__main__":
    main()
