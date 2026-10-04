#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 61_replicative_ceiling.py
#
# Lift anchors hg19->hg18 and quantify every arm in its own build.
#
# CONSUMES: output/enhancer_test/anchors.csv, data/*/bw/
# PRODUCES: output/replicative_ceiling_test/anchors_with_hg18.csv, signal_matrix_replicative.csv, quantification_stats.json, STEP1_LOG.txt
# ---------------------------------------------------------------------------
"""
WITHIN-REPLICATIVE CEILING TEST -- the control missing from the prior cross-mechanism run.

The prior test established an OIS ceiling of Spearman +0.653 (GM21 vs IMR90 H3K27ac,
different labs/cell lines/protocols) but had NO within-replicative ceiling, because the
only second replicative dataset (GSE210285) has no peak files. Without it, "replicative
and OIS enhancer responses do not correlate" could not be separated from "the WI-38 ATAC
measurement cannot detect anything."

DESIGN: the prior run's 237,824 anchors are REUSED verbatim (1 kb, hg19, >=2-sample
consensus from GM21 ATAC+H3K27ac). That makes every number here directly comparable to
+0.653 -- same regions, same widths, same normalization, same thresholds. The anchors are
defined by peak PRESENCE in fibroblasts (either arm), not by OIS responsiveness, so they
are a general fibroblast regulatory atlas rather than an OIS-biased selection; that is
stated as a caveat rather than assumed harmless.

BUILDS: anchors are hg19. GSE106146 is hg19 (queried directly). GSE146585 is hg18 (all 223
samples) -- anchors are lifted hg19->hg18 to query those bigWigs, so signal is always read
in the file's own build and liftover never touches a coverage value.

n IS SMALL EVERYWHERE. n=2/arm for IMR90 replicative; n=1/arm for BJ replicative, IMR90
OIS, and GSE106146 replicative. Stated alongside every correlation, not in a footnote.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
DATA = ROOT / "data"
OUT = ROOT / "output" / "replicative_ceiling_test"
OUT.mkdir(parents=True, exist_ok=True)
PRIOR = ROOT / "output" / "enhancer_test" / "region_responses.csv"

S = "_STARdef_unique_tbp1.ucsc.bigWig"
BW = {
    # label -> (path, build, arm, dataset)
    "IMR90rep_Young_R1":  (DATA/"GSE146585/bw"/f"GSM4395479_TSD311_IMR90_Young_anti-H3K27ac_Rep1{S}", "hg18", "PRO", "GSE146567_IMR90_replicative"),
    "IMR90rep_Young_R2":  (DATA/"GSE146585/bw"/f"GSM4395480_TSD312_IMR90_Young_anti-H3K27ac_Rep2{S}", "hg18", "PRO", "GSE146567_IMR90_replicative"),
    "IMR90rep_Sen_R1":    (DATA/"GSE146585/bw"/f"GSM4395481_TSD313_IMR90_Senescent_anti-H3K27ac_Rep1{S}", "hg18", "SEN", "GSE146567_IMR90_replicative"),
    "IMR90rep_Sen_R2":    (DATA/"GSE146585/bw"/f"GSM4395482_TSD314_IMR90_Senescent_anti-H3K27ac_Rep2{S}", "hg18", "SEN", "GSE146567_IMR90_replicative"),
    "BJrep_Young":        (DATA/"GSE146585/bw"/f"GSM4395086_Tom161_BJ_Young_anti-H3K27ac{S}", "hg18", "PRO", "GSE146559_BJ_replicative"),
    "BJrep_Sen":          (DATA/"GSE146585/bw"/f"GSM4395091_Tom164_BJ_Senescent_anti-H3K27ac{S}", "hg18", "SEN", "GSE146559_BJ_replicative"),
    "IMR90ois_GFP":       (DATA/"GSE146585/bw"/f"GSM4395451_TSC801_IMR90_GFP_Control_anti-H3K27ac{S}", "hg18", "PRO", "GSE146563_IMR90_OIS"),
    "IMR90ois_RAS":       (DATA/"GSE146585/bw"/f"GSM4395452_TSC802_IMR90_RAS_OIS_anti-H3K27ac{S}", "hg18", "SEN", "GSE146563_IMR90_OIS"),
    "Sen2019_Pro":        (DATA/"GSE106146/bw"/"GSE106146_P2-H3K27ac.bigWig", "hg19", "PRO", "GSE106146_replicative"),
    "Sen2019_Sen":        (DATA/"GSE106146/bw"/"GSE106146_S1-H3K27ac.bigWig", "hg19", "SEN", "GSE106146_replicative"),
    "Sen2019_IR_R1":      (DATA/"GSE106146/bw"/"GSM2830447_IR-H3K27ac.Rep1.bigWig", "hg19", "IR", "GSE106146_IR"),
    "Sen2019_IR_R2":      (DATA/"GSE106146/bw"/"GSM2830448_IR-H3K27ac.Rep2.bigWig", "hg19", "IR", "GSE106146_IR"),
}

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    stats = {}
    R("=" * 84)
    R("WITHIN-REPLICATIVE CEILING TEST -- quantification")
    R("=" * 84)

    prior = pd.read_csv(PRIOR)
    R(f"\nreusing prior anchor set: {len(prior):,} anchors (1 kb, hg19)")
    adf = prior[["chrom_hg19", "start_hg19", "end_hg19", "mid_hg19",
                 "in_gm21_atac", "in_gm21_k27ac", "wi38_annot", "wi38_gene",
                 "gm21_atac", "gm21_k27ac", "imr90_k27ac", "wi38_atac"]].copy()

    # ---- lift anchors hg19 -> hg18 so hg18 bigWigs are read in their own build ----
    R("\n[1] LIFT ANCHORS hg19 -> hg18 (GSE146585 is hg18)\n" + "-" * 84)
    from pyliftover import LiftOver
    lo = LiftOver(str(DATA / "annotation" / "hg19ToHg18.over.chain.gz"))
    m18, ok18 = [], []
    for c, m in zip(adf["chrom_hg19"], adf["mid_hg19"]):
        r = lo.convert_coordinate(c, int(m))
        if r and r[0][0] == c:
            m18.append(r[0][1]); ok18.append(True)
        else:
            m18.append(-1); ok18.append(False)
    adf["mid_hg18"] = m18
    adf["lifted_hg18"] = ok18
    nf = int((~adf["lifted_hg18"]).sum())
    R(f"  midpoints lifted hg19->hg18: {len(adf)-nf:,}/{len(adf):,} "
      f"({100*(len(adf)-nf)/len(adf):.2f}%)  FAILED {nf:,} ({100*nf/len(adf):.2f}%)")
    stats["liftover_hg19_to_hg18_pct"] = float(100 * (len(adf) - nf) / len(adf))
    stats["liftover_hg19_to_hg18_failed"] = nf

    rng = np.random.default_rng(0)
    samp = rng.choice(len(adf), size=min(5000, len(adf)), replace=False)
    both, span_ok = 0, 0
    for i in samp:
        c = adf["chrom_hg19"].iloc[i]
        rs = lo.convert_coordinate(c, int(adf["start_hg19"].iloc[i]))
        re_ = lo.convert_coordinate(c, int(adf["end_hg19"].iloc[i]))
        if rs and re_ and rs[0][0] == c and re_[0][0] == c:
            both += 1
            if abs(abs(re_[0][1] - rs[0][1]) - 1000) <= 50:
                span_ok += 1
    R(f"  span integrity on {len(samp):,} sampled anchors: both ends lift {both:,} "
      f"({100*both/len(samp):.2f}%); span within 50 bp of 1 kb: {span_ok:,} "
      f"({100*span_ok/max(both,1):.2f}%)")
    stats["hg18_span_integrity_pct"] = float(100 * span_ok / max(both, 1))

    # ---- quantify ----
    R("\n[2] BIGWIG QUANTIFICATION (each file read in its OWN build)\n" + "-" * 84)
    import pyBigWig
    c19 = adf["chrom_hg19"].to_numpy(); s19 = adf["start_hg19"].to_numpy()
    e19 = adf["end_hg19"].to_numpy()
    m18a = adf["mid_hg18"].to_numpy(); ok18a = adf["lifted_hg18"].to_numpy()
    sig = {}
    for label, (path, build, arm, ds) in BW.items():
        if not path.exists():
            R(f"  MISSING {label}: {path.name}")
            continue
        bw = pyBigWig.open(str(path))
        avail = set(bw.chroms().keys())
        vals = np.full(len(adf), np.nan)
        for i in range(len(adf)):
            c = c19[i]
            cc = c if c in avail else (c[3:] if c[3:] in avail else None)
            if cc is None:
                continue
            if build == "hg18":
                if not ok18a[i]:
                    continue
                st, en = int(m18a[i]) - 500, int(m18a[i]) + 500
            else:
                st, en = int(s19[i]), int(e19[i])
            if st < 0:
                continue
            try:
                v = bw.stats(cc, st, en, type="mean")[0]
            except Exception:
                v = None
            if v is not None:
                vals[i] = v
        bw.close()
        sig[label] = vals
        R(f"  {label:20s} {build}  {ds:30s} covered "
          f"{int(np.isfinite(vals).sum()):,}/{len(adf):,}  mean {np.nanmean(vals):.4f}")

    sm = pd.DataFrame(sig)
    adf.to_csv(OUT / "anchors_with_hg18.csv", index=False)
    sm.to_csv(OUT / "signal_matrix_replicative.csv", index=False)
    json.dump(stats, open(OUT / "quantification_stats.json", "w"), indent=2)
    (OUT / "STEP1_LOG.txt").write_text("\n".join(LOG))
    R(f"\n  wrote signal_matrix_replicative.csv {sm.shape}")
    R("QUANT_DONE")


if __name__ == "__main__":
    main()
