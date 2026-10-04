#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 70_build_ois_regions.py
#
# Build the 112,099 K27ac-summit consensus regions; attach ATAC design anchors; lift to hg38.
#
# CONSUMES: data/GSE205898/peaks/, data/GSE206402/peaks/
# PRODUCES: output/ois_enhancer_trainset/regions_hg19_hg38.csv, region_stats.json, STEP1_LOG.txt
# ---------------------------------------------------------------------------
"""
OIS-scoped enhancer training set -- STEP 1: regions + signal quantification.

SCOPE: OIS only. Replicative senescence was excluded by a prior kill test (within-
replicative H3K27ac reproducibility +0.054 same-lab / +0.144 cross-lab / +0.268, all inside
the -0.164..+0.164 cross-mechanism range). OIS reaches +0.653 across labs and tissues.

PRIMARY TARGET IS H3K27ac, accessibility is a filter -- 35% of differential-accessibility
regions show no H3K27ac change, and where both move they agree in direction 93% of the time.

REGION DEFINITION
  Source: GSE205898 H3K27ac narrowPeaks (the only deposited H3K27ac peak calls; GSE74238
  has bigWigs only). Consensus = summit supported by >= 2 of the 4 samples within
  SUMMIT_CLUSTER bp. Region = fixed REGION_W bp centred on the cluster's median summit.

  WHY SUMMIT-CENTRED AND WHY THIS WIDTH: H3K27ac marks the nucleosomes FLANKING a regulatory
  element and is depleted over the nucleosome-free region where transcription factors
  actually bind. So the H3K27ac summit is where the acetylation signal is maximal, which is
  the right place to MEASURE, but it is NOT where a 200 bp designed element should sit. The
  1 kb measurement window captures the flanking acetylated nucleosomes; the eventual 200 bp
  design coordinate is the nucleosome-free region, which is located by the ATAC summit. This
  script therefore records, per region, the nearest GM21 ATAC summit as a separate
  design-anchor column. Conflating the two would place designed sequence on a nucleosome.

BUILDS: every bigWig here is hg19, so quantification happens in hg19 and no coverage value
is ever read through a lifted coordinate. Region midpoints are lifted hg19 -> hg38 only to
emit Borzoi-compatible intervals (Borzoi is hg38).

REPLICATE COUNTS (stated with every downstream number): GM21 EV n=2, GM21 RAS_D18 n=2,
IMR90 proliferating n=2, IMR90 senescent n=2, IMR90 quiescent n=2, GM21 ATAC 4 EV / 4 RAS.
"""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
DATA = ROOT / "data"
OUT = ROOT / "output" / "ois_enhancer_trainset"
OUT.mkdir(parents=True, exist_ok=True)

REGION_W = 1000          # measurement window; see docstring
SUMMIT_CLUSTER = 250     # summits within this distance are the same element
MIN_SUPPORT = 2          # of 4 GSE205898 K27ac samples
STD = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

K27_PEAKS = {  # GSE205898, hg19
    "EV_R1": "GSM6235098_EV_K27ac_Rep1_peaks.narrowPeak.gz",
    "EV_R2": "GSM6235099_EV_K27ac_Rep2_peaks.narrowPeak.gz",
    "RAS_R1": "GSM6235110_RAS_D18_K27ac_Rep1_peaks.narrowPeak.gz",
    "RAS_R2": "GSM6235111_RAS_D18_K27ac_Rep2_peaks.narrowPeak.gz",
}
BW = {
    # GM21 H3K27ac (GSE205898) -- TRAINING study
    "GM21_K27_EV_R1":  (DATA/"GSE205898/bw/GSM6235098_EV_K27ac_1_dedup_blacklisted_merge.bw", "GM21_K27", "EV"),
    "GM21_K27_EV_R2":  (DATA/"GSE205898/bw/GSM6235099_EV_K27ac_2_dedup_blacklisted_merge.bw", "GM21_K27", "EV"),
    "GM21_K27_RAS_R1": (DATA/"GSE205898/bw/GSM6235110_RAS_D18_K27ac_1_dedup_blacklisted_merge.bw", "GM21_K27", "RAS"),
    "GM21_K27_RAS_R2": (DATA/"GSE205898/bw/GSM6235111_RAS_D18_K27ac_2_dedup_blacklisted_merge.bw", "GM21_K27", "RAS"),
    # IMR90 H3K27ac (GSE74238) -- VALIDATION study, includes the quiescent arm
    "IMR90_K27_PRO_R1": (DATA/"GSE74238/bw/GSM1915113_H3K27ac_proliferating.bigWig", "IMR90_K27", "PRO"),
    "IMR90_K27_PRO_R2": (DATA/"GSE74238/bw/GSM2098176_H3K27ac_proliferating_replicate.bigWig", "IMR90_K27", "PRO"),
    "IMR90_K27_SEN_R1": (DATA/"GSE74238/bw/GSM1915115_H3K27ac_senescent..bigWig", "IMR90_K27", "SEN"),
    "IMR90_K27_SEN_R2": (DATA/"GSE74238/bw/GSM2098178_H3K27ac_senescent_replicate.bigWig", "IMR90_K27", "SEN"),
    "IMR90_K27_QUI_R1": (DATA/"GSE74238/bw/GSM1915114_H3K27ac_quiescent..bigWig", "IMR90_K27", "QUI"),
    "IMR90_K27_QUI_R2": (DATA/"GSE74238/bw/GSM2098177_H3K27ac_quiescent_replicate.bigWig", "IMR90_K27", "QUI"),
    # GM21 ATAC (GSE206402) -- accessibility filter + design anchor
    "GM21_ATAC_EV_R1":  (DATA/"GSE206402/bw/GSM6253126_pBABE_D8_REP1.bw", "GM21_ATAC", "EV"),
    "GM21_ATAC_EV_R2":  (DATA/"GSE206402/bw/GSM6253127_pBABE_D8_REP2.bw", "GM21_ATAC", "EV"),
    "GM21_ATAC_EV_R3":  (DATA/"GSE206402/bw/GSM6253122_pBABE_D32_REP1.bw", "GM21_ATAC", "EV"),
    "GM21_ATAC_EV_R4":  (DATA/"GSE206402/bw/GSM6253123_pBABE_D32_REP2.bw", "GM21_ATAC", "EV"),
    "GM21_ATAC_RAS_R1": (DATA/"GSE206402/bw/GSM6253132_RAS_D18_REP1.bw", "GM21_ATAC", "RAS"),
    "GM21_ATAC_RAS_R2": (DATA/"GSE206402/bw/GSM6253133_RAS_D18_REP2.bw", "GM21_ATAC", "RAS"),
    "GM21_ATAC_RAS_R3": (DATA/"GSE206402/bw/GSM6253134_RAS_D23_REP1.bw", "GM21_ATAC", "RAS"),
    "GM21_ATAC_RAS_R4": (DATA/"GSE206402/bw/GSM6253135_RAS_D23_REP2.bw", "GM21_ATAC", "RAS"),
}

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def read_summits(path):
    """narrowPeak -> [(chrom, summit_abs)]; col10 is the summit offset from start."""
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[0] not in STD:
                continue
            out.append((f[0], int(f[1]) + int(f[9])))
    return out


def main():
    stats = {}
    R("=" * 84)
    R("OIS ENHANCER TRAINING SET -- STEP 1: regions + quantification")
    R("=" * 84)
    R(f"region width {REGION_W} bp (summit-centred) | summit cluster {SUMMIT_CLUSTER} bp | "
      f"min support {MIN_SUPPORT}/4")

    # ---------------- consensus H3K27ac regions ----------------
    R("\n[1] CONSENSUS H3K27ac REGIONS from GSE205898 narrowPeaks (hg19)\n" + "-" * 84)
    per_sample = {}
    for lab, fn in K27_PEAKS.items():
        p = DATA / "GSE205898" / "peaks" / fn
        s = read_summits(p)
        per_sample[lab] = s
        R(f"  {lab:8s} {len(s):7,} peak summits")
    by_chrom = {}
    for lab, s in per_sample.items():
        for c, pos in s:
            by_chrom.setdefault(c, []).append((pos, lab))
    regions = []
    for c, lst in by_chrom.items():
        lst.sort()
        i = 0
        while i < len(lst):
            j = i
            members = {lst[i][1]}
            while j + 1 < len(lst) and lst[j + 1][0] - lst[j][0] <= SUMMIT_CLUSTER:
                j += 1
                members.add(lst[j][1])
            if len(members) >= MIN_SUPPORT:
                pos = int(np.median([lst[k][0] for k in range(i, j + 1)]))
                regions.append((c, pos, len(members)))
            i = j + 1
    regions.sort()
    R(f"\n  consensus regions (summit supported by >= {MIN_SUPPORT} of 4 samples): "
      f"{len(regions):,}")
    stats["n_consensus_regions"] = len(regions)
    stats["peaks_per_sample"] = {k: len(v) for k, v in per_sample.items()}

    rdf = pd.DataFrame(regions, columns=["chrom", "summit_hg19", "n_samples_supporting"])
    rdf["start_hg19"] = rdf["summit_hg19"] - REGION_W // 2
    rdf["end_hg19"] = rdf["summit_hg19"] + REGION_W // 2
    rdf = rdf[rdf["start_hg19"] >= 0].reset_index(drop=True)
    R(f"  after dropping regions with negative start: {len(rdf):,}")

    # ---------------- ATAC summits as design anchors ----------------
    R("\n[2] GM21 ATAC SUMMITS as 200 bp DESIGN ANCHORS\n" + "-" * 84)
    R("  H3K27ac marks the nucleosomes flanking an element and is depleted over the")
    R("  nucleosome-free region where TFs bind. The K27ac summit is the right place to")
    R("  MEASURE but the wrong place to PLACE 200 bp of designed sequence. The nearest ATAC")
    R("  summit is recorded per region as the design coordinate.")
    atac_sum = []
    pk_dir = DATA / "GSE206402" / "peaks"
    n_files = 0
    for p in sorted(pk_dir.glob("*.narrowPeak.gz")):
        atac_sum.extend(read_summits(p))
        n_files += 1
    R(f"  read {n_files} GM21 ATAC narrowPeak files -> {len(atac_sum):,} summits")
    ab = {}
    for c, pos in atac_sum:
        ab.setdefault(c, []).append(pos)
    for c in ab:
        ab[c] = np.array(sorted(ab[c]))
    nearest, dist = np.full(len(rdf), -1, np.int64), np.full(len(rdf), np.iinfo(np.int32).max, np.int64)
    for i, (c, s) in enumerate(zip(rdf["chrom"], rdf["summit_hg19"])):
        arr = ab.get(c)
        if arr is None or not len(arr):
            continue
        j = np.searchsorted(arr, s)
        cand = [arr[k] for k in (j - 1, j) if 0 <= k < len(arr)]
        if cand:
            best = min(cand, key=lambda x: abs(x - s))
            nearest[i], dist[i] = best, abs(best - s)
    rdf["atac_summit_hg19"] = nearest
    rdf["atac_summit_dist"] = dist
    within = int((rdf["atac_summit_dist"] <= REGION_W // 2).sum())
    R(f"  regions with an ATAC summit inside the {REGION_W} bp window: {within:,}/{len(rdf):,} "
      f"({100*within/len(rdf):.1f}%)")
    R(f"  median distance K27ac summit -> nearest ATAC summit: "
      f"{int(np.median(rdf['atac_summit_dist'][rdf['atac_summit_dist']<1e9])):,} bp")
    stats["regions_with_atac_summit_inside"] = within

    # ---------------- liftover hg19 -> hg38 ----------------
    R("\n[3] LIFTOVER hg19 -> hg38 (for Borzoi intervals; signal stays in hg19)\n" + "-" * 84)
    from pyliftover import LiftOver
    lo = LiftOver(str(DATA / "annotation" / "hg19ToHg38.over.chain.gz"))
    m38, ok = [], []
    for c, s in zip(rdf["chrom"], rdf["summit_hg19"]):
        r = lo.convert_coordinate(c, int(s))
        if r and r[0][0] == c:
            m38.append(r[0][1]); ok.append(True)
        else:
            m38.append(-1); ok.append(False)
    rdf["summit_hg38"] = m38
    rdf["lifted"] = ok
    nf = int((~rdf["lifted"]).sum())
    R(f"  summits lifted: {len(rdf)-nf:,}/{len(rdf):,} ({100*(len(rdf)-nf)/len(rdf):.2f}%)  "
      f"FAILED {nf:,} ({100*nf/len(rdf):.2f}%)")
    stats["liftover_pct"] = float(100 * (len(rdf) - nf) / len(rdf))
    stats["liftover_failed"] = nf
    # ATAC design anchor in hg38 too
    a38 = []
    for c, s, d in zip(rdf["chrom"], rdf["atac_summit_hg19"], rdf["atac_summit_dist"]):
        if s < 0:
            a38.append(-1); continue
        r = lo.convert_coordinate(c, int(s))
        a38.append(r[0][1] if (r and r[0][0] == c) else -1)
    rdf["atac_summit_hg38"] = a38

    # ---------------- quantify ----------------
    R("\n[4] SIGNAL QUANTIFICATION (all bigWigs are hg19; read in hg19)\n" + "-" * 84)
    import pyBigWig
    chrom = rdf["chrom"].to_numpy()
    st = rdf["start_hg19"].to_numpy()
    en = rdf["end_hg19"].to_numpy()
    sig = {}
    for lab, (path, assay, arm) in BW.items():
        if not path.exists():
            R(f"  MISSING {lab}: {path.name}")
            continue
        bw = pyBigWig.open(str(path))
        avail = set(bw.chroms().keys())
        vals = np.full(len(rdf), np.nan)
        for i in range(len(rdf)):
            c = chrom[i]
            cc = c if c in avail else (c[3:] if c[3:] in avail else None)
            if cc is None:
                continue
            try:
                v = bw.stats(cc, int(st[i]), int(en[i]), type="mean")[0]
            except Exception:
                v = None
            if v is not None:
                vals[i] = v
        bw.close()
        sig[lab] = vals
        R(f"  {lab:20s} {assay:10s} {arm:4s} covered {int(np.isfinite(vals).sum()):,}"
          f"/{len(rdf):,}  mean {np.nanmean(vals):.4f}")

    sm = pd.DataFrame(sig)
    rdf.to_csv(OUT / "regions_hg19_hg38.csv", index=False)
    sm.to_csv(OUT / "signal_matrix_ois.csv", index=False)
    json.dump(stats, open(OUT / "region_stats.json", "w"), indent=2)
    (OUT / "STEP1_LOG.txt").write_text("\n".join(LOG))
    R(f"\n  wrote regions_hg19_hg38.csv {rdf.shape} and signal_matrix_ois.csv {sm.shape}")
    R("STEP1_DONE")


if __name__ == "__main__":
    main()
