#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 51_build_enhancer_regions.py
#
# Build the 237,824-anchor consensus set and lift hg19->hg38.
#
# CONSUMES: data/GSE205898/peaks/, GSE206402/peaks/
# PRODUCES: output/enhancer_test/anchors.csv, harmonization_stats.json, STEP1_LOG.txt
# ---------------------------------------------------------------------------
"""
Cross-mechanism enhancer kill test -- STEP 1: harmonize regions and quantify signal.

DESIGN DECISIONS, stated (see REPORT for full rationale):

 * ANCHOR SPACE IS hg19. Three of the four studies (GSE206402, GSE205898, GSE74238) are
   hg19 and EVERY bigWig is hg19; only GSE175533 is hg38. Anchoring in hg19 means the
   signal quantification never touches a lifted coordinate -- liftover is used ONLY to
   carry anchors into hg38 to look up GSE175533's deposited response. Anchoring in hg38
   would have required lifting hg38->hg19 to query the bigWigs, putting a liftover in the
   path of every measurement.

 * FIXED-WIDTH RECENTERING (1 kb around the peak midpoint). Peak widths differ by caller
   and pipeline (MACS2 settings differ between these studies), so raw peak intervals are
   not comparable units. A fixed width makes coverage per-base comparable and removes
   peak-width as a confounder. It also means only the MIDPOINT needs lifting, which
   matters: the commonly-cited ~99.9% liftover success rate is for point positions, not
   wide intervals. Span integrity of the lift is checked separately by lifting both ends.

 * CONSENSUS = called in >= MIN_SUPPORT samples of that assay, merged across samples.
   Requiring reproducibility across samples rather than taking a union of single-sample
   calls, which would be dominated by one-off calls.

Outputs: output/enhancer_test/anchors.csv, signal_matrix.csv, harmonization_stats.json
"""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
DATA = ROOT / "data"
OUT = ROOT / "output" / "enhancer_test"
OUT.mkdir(parents=True, exist_ok=True)

WIDTH = 1000            # fixed-width anchor
MIN_SUPPORT = 2         # peak must be called in >= this many samples of that assay
DEDUP_DIST = 500        # anchors whose midpoints are closer than this are merged
STD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


# ------------------------------------------------------------------ sample tables
GM21_ATAC = {  # GSE206402, hg19. RAS D32+ excluded: this study's subject is ESCAPE from
               # OIS and late timepoints revert (established in the expression build).
    "EV":  [("GSM6253126", "pBABE_D8_Rep1"), ("GSM6253127", "pBABE_D8_Rep2"),
            ("GSM6253122", "pBABE_D32_Rep1"), ("GSM6253123", "pBABE_D32_Rep2")],
    "SEN": [("GSM6253132", "RAS_D18_Rep1"), ("GSM6253133", "RAS_D18_Rep2"),
            ("GSM6253134", "RAS_D23_Rep1"), ("GSM6253135", "RAS_D23_Rep2")],
}
GM21_ATAC_BW = {"GSM6253126": "GSM6253126_pBABE_D8_REP1.bw",
                "GSM6253127": "GSM6253127_pBABE_D8_REP2.bw",
                "GSM6253122": "GSM6253122_pBABE_D32_REP1.bw",
                "GSM6253123": "GSM6253123_pBABE_D32_REP2.bw",
                "GSM6253132": "GSM6253132_RAS_D18_REP1.bw",
                "GSM6253133": "GSM6253133_RAS_D18_REP2.bw",
                "GSM6253134": "GSM6253134_RAS_D23_REP1.bw",
                "GSM6253135": "GSM6253135_RAS_D23_REP2.bw"}

GM21_K27 = {   # GSE205898, hg19
    "EV":  [("GSM6235098", "EV_K27ac_Rep1"), ("GSM6235099", "EV_K27ac_Rep2")],
    "SEN": [("GSM6235110", "RAS_D18_K27ac_Rep1"), ("GSM6235111", "RAS_D18_K27ac_Rep2")],
}
GM21_K27_BW = {"GSM6235098": "GSM6235098_EV_K27ac_1_dedup_blacklisted_merge.bw",
               "GSM6235099": "GSM6235099_EV_K27ac_2_dedup_blacklisted_merge.bw",
               "GSM6235110": "GSM6235110_RAS_D18_K27ac_1_dedup_blacklisted_merge.bw",
               "GSM6235111": "GSM6235111_RAS_D18_K27ac_2_dedup_blacklisted_merge.bw"}

IMR90_K27_BW = {  # GSE74238, hg19. NO peaks were ever deposited for this series.
    "PRO": ["GSM1915113_H3K27ac_proliferating.bigWig",
            "GSM2098176_H3K27ac_proliferating_replicate.bigWig"],
    "SEN": ["GSM1915115_H3K27ac_senescent..bigWig",
            "GSM2098178_H3K27ac_senescent_replicate.bigWig"],
}


def read_narrowpeak(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[0] not in STD_CHROMS:
                continue
            rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def consensus(sample_intervals, min_support):
    """Sweep-line: keep genomic segments covered by >= min_support samples.
    Each sample contributes its own merged intervals so one sample cannot self-support."""
    per_chrom = {}
    for iv in sample_intervals:
        # merge within sample first
        by_c = {}
        for c, s, e in iv:
            by_c.setdefault(c, []).append((s, e))
        for c, lst in by_c.items():
            lst.sort()
            merged, cs, ce = [], None, None
            for s, e in lst:
                if cs is None:
                    cs, ce = s, e
                elif s <= ce:
                    ce = max(ce, e)
                else:
                    merged.append((cs, ce)); cs, ce = s, e
            if cs is not None:
                merged.append((cs, ce))
            per_chrom.setdefault(c, []).extend(merged)
    out = []
    for c, lst in per_chrom.items():
        events = []
        for s, e in lst:
            events.append((s, 1)); events.append((e, -1))
        events.sort()
        depth, start = 0, None
        for pos, d in events:
            prev = depth
            depth += d
            if prev < min_support <= depth:
                start = pos
            elif prev >= min_support > depth and start is not None:
                out.append((c, start, pos)); start = None
    return out


def to_anchors(segments, width=WIDTH):
    """Fixed-width windows centred on each segment midpoint."""
    a = []
    for c, s, e in segments:
        m = (s + e) // 2
        a.append((c, m - width // 2, m + width // 2, m))
    return a


def dedup(anchors, dist=DEDUP_DIST):
    by_c = {}
    for c, s, e, m in anchors:
        by_c.setdefault(c, []).append(m)
    out = []
    for c, ms in by_c.items():
        ms.sort()
        last = None
        for m in ms:
            if last is None or m - last >= dist:
                out.append((c, m - WIDTH // 2, m + WIDTH // 2, m))
                last = m
    return out


def main():
    stats = {}
    R("=" * 84)
    R("CROSS-MECHANISM ENHANCER TEST -- STEP 1: harmonization + quantification")
    R("=" * 84)
    R(f"anchor width {WIDTH} bp | min sample support {MIN_SUPPORT} | "
      f"dedup midpoint distance {DEDUP_DIST} bp | anchor space hg19")

    # ---------------- consensus per assay ----------------
    R("\n[1] CONSENSUS PEAK SETS (hg19)\n" + "-" * 84)
    atac_iv = []
    for arm, lst in GM21_ATAC.items():
        for gsm, tag in lst:
            p = DATA / "GSE206402" / "peaks" / f"{gsm}_{tag}_peaks.narrowPeak.gz"
            if not p.exists():
                cand = list((DATA / "GSE206402" / "peaks").glob(f"{gsm}_*.narrowPeak.gz"))
                p = cand[0] if cand else None
            if p is None:
                R(f"  MISSING narrowPeak for {gsm} {tag}"); continue
            iv = read_narrowpeak(p)
            atac_iv.append(iv)
            R(f"  GM21 ATAC  {arm:3s} {tag:20s} {len(iv):7,} peaks")
    k27_iv = []
    for arm, lst in GM21_K27.items():
        for gsm, tag in lst:
            p = DATA / "GSE205898" / "peaks" / f"{gsm}_{tag}_peaks.narrowPeak.gz"
            if not p.exists():
                R(f"  MISSING narrowPeak {p.name}"); continue
            iv = read_narrowpeak(p)
            k27_iv.append(iv)
            R(f"  GM21 K27ac {arm:3s} {tag:20s} {len(iv):7,} peaks")

    atac_seg = consensus(atac_iv, MIN_SUPPORT)
    k27_seg = consensus(k27_iv, MIN_SUPPORT)
    R(f"\n  GM21 ATAC  consensus segments (>= {MIN_SUPPORT} samples): {len(atac_seg):,}")
    R(f"  GM21 K27ac consensus segments (>= {MIN_SUPPORT} samples): {len(k27_seg):,}")
    stats["gm21_atac_consensus_segments"] = len(atac_seg)
    stats["gm21_k27ac_consensus_segments"] = len(k27_seg)

    anchors = dedup(to_anchors(atac_seg) + to_anchors(k27_seg))
    anchors.sort()
    R(f"  anchors after {WIDTH} bp recentering + dedup: {len(anchors):,}")
    stats["n_anchors"] = len(anchors)

    adf = pd.DataFrame(anchors, columns=["chrom_hg19", "start_hg19", "end_hg19", "mid_hg19"])
    # which assay supports each anchor
    def mark(segs, name):
        by_c = {}
        for c, s, e in segs:
            by_c.setdefault(c, []).append((s, e))
        for c in by_c:
            by_c[c].sort()
        flag = np.zeros(len(adf), bool)
        for i, (c, m) in enumerate(zip(adf["chrom_hg19"], adf["mid_hg19"])):
            lst = by_c.get(c)
            if not lst:
                continue
            lo, hi = 0, len(lst) - 1
            while lo <= hi:
                mid = (lo + hi) // 2
                s, e = lst[mid]
                if m < s:
                    hi = mid - 1
                elif m >= e:
                    lo = mid + 1
                else:
                    flag[i] = True; break
        adf[name] = flag
    mark(atac_seg, "in_gm21_atac")
    mark(k27_seg, "in_gm21_k27ac")
    R(f"  anchors in GM21 ATAC consensus : {int(adf['in_gm21_atac'].sum()):,}")
    R(f"  anchors in GM21 K27ac consensus: {int(adf['in_gm21_k27ac'].sum()):,}")
    R(f"  anchors in BOTH (ATAC & K27ac) : {int((adf['in_gm21_atac']&adf['in_gm21_k27ac']).sum()):,}")
    stats["anchors_atac"] = int(adf["in_gm21_atac"].sum())
    stats["anchors_k27ac"] = int(adf["in_gm21_k27ac"].sum())
    stats["anchors_both"] = int((adf["in_gm21_atac"] & adf["in_gm21_k27ac"]).sum())

    # ---------------- liftover hg19 -> hg38 ----------------
    R("\n[2] LIFTOVER hg19 -> hg38 (midpoints; span integrity checked separately)\n"
      + "-" * 84)
    from pyliftover import LiftOver
    lo = LiftOver(str(DATA / "annotation" / "hg19ToHg38.over.chain.gz"))
    mid38, ok = [], []
    for c, m in zip(adf["chrom_hg19"], adf["mid_hg19"]):
        r = lo.convert_coordinate(c, int(m))
        if r and r[0][0] == c:
            mid38.append(r[0][1]); ok.append(True)
        elif r:
            mid38.append(-1); ok.append(False)      # lifted to a different chromosome
        else:
            mid38.append(-1); ok.append(False)
    adf["mid_hg38"] = mid38
    adf["lifted"] = ok
    nfail = int((~adf["lifted"]).sum())
    R(f"  midpoints lifted: {len(adf)-nfail:,}/{len(adf):,} "
      f"({100*(len(adf)-nfail)/len(adf):.2f}%)   FAILED: {nfail:,} "
      f"({100*nfail/len(adf):.2f}%)")
    stats["liftover_midpoint_success_pct"] = float(100 * (len(adf) - nfail) / len(adf))
    stats["liftover_midpoint_failed"] = nfail

    # span integrity: lift both ends for a random sample and check the width survives
    rng = np.random.default_rng(0)
    samp = rng.choice(len(adf), size=min(5000, len(adf)), replace=False)
    good, bad, span_ok = 0, 0, 0
    for i in samp:
        c = adf["chrom_hg19"].iloc[i]
        s, e = int(adf["start_hg19"].iloc[i]), int(adf["end_hg19"].iloc[i])
        rs, re_ = lo.convert_coordinate(c, s), lo.convert_coordinate(c, e)
        if rs and re_ and rs[0][0] == c and re_[0][0] == c:
            good += 1
            if abs(abs(re_[0][1] - rs[0][1]) - WIDTH) <= 50:
                span_ok += 1
        else:
            bad += 1
    R(f"  span-integrity check on {len(samp):,} sampled anchors: both ends lift "
      f"{good:,} ({100*good/len(samp):.2f}%); of those, lifted span within 50 bp of "
      f"{WIDTH} bp: {span_ok:,} ({100*span_ok/max(good,1):.2f}%)")
    R("  (a midpoint that lifts while the span does not is the wide-interval failure mode "
      "the\n   99.9%-for-point-variants caveat refers to; quantified here rather than "
      "assumed away)")
    stats["span_integrity"] = dict(sampled=len(samp), both_ends_lift=good,
                                   span_preserved=span_ok,
                                   pct_span_preserved=float(100 * span_ok / max(good, 1)))

    # ---------------- WI-38 response from the deposited differential table ----------
    R("\n[3] WI-38 REPLICATIVE ATAC RESPONSE (GSE175533, hg38, deposited log2FC)\n"
      + "-" * 84)
    import openpyxl
    wb = openpyxl.load_workbook(DATA / "GSE175533" / "GSE175533_atac_peaks_sig.transitional.xlsx",
                                read_only=True)
    ws = wb["atac_peaks"]
    rows = ws.iter_rows(values_only=True)
    hdr = list(next(rows))
    wi = pd.DataFrame(list(rows), columns=hdr)
    R(f"  loaded {len(wi):,} regions x {len(hdr)} columns")
    OFFC = ["PDL25_v_htert2_log2FC", "PDL33_v_htert4_log2FC", "PDL37_v_htert5_log2FC"]
    ONC = "PDL50_v_htert7_log2FC"
    for c in OFFC + [ONC]:
        wi[c] = pd.to_numeric(wi[c], errors="coerce")
    # Difference of differences: each arm is already log2FC vs its own matched hTERT
    # control, so subtracting cancels the hTERT baseline. PDL45 is EXCLUDED as
    # transitional -- the prior audit placed the senescence inflection at PDL46-50.
    wi["wi38_response"] = wi[ONC] - wi[OFFC].mean(axis=1)
    wi = wi.dropna(subset=["wi38_response"])
    parts = wi["region"].str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
    wi["chrom38"] = parts[0]
    wi["start38"] = pd.to_numeric(parts[1])
    wi["end38"] = pd.to_numeric(parts[2])
    wi = wi.dropna(subset=["chrom38"])
    R(f"  usable regions with response: {len(wi):,}")
    R(f"  response = {ONC} - mean({', '.join(OFFC)})   (PDL45 excluded as transitional)")
    R(f"  response range [{wi['wi38_response'].min():.3f}, {wi['wi38_response'].max():.3f}] "
      f"sd {wi['wi38_response'].std():.3f}")

    # match lifted anchors into WI-38 peaks (midpoint inside peak interval)
    wi_by_c = {}
    for c, g in wi.groupby("chrom38"):
        g = g.sort_values("start38")
        wi_by_c[c] = (g["start38"].to_numpy(), g["end38"].to_numpy(),
                      g["wi38_response"].to_numpy(), g["annot"].to_numpy(),
                      g["gene_name"].to_numpy())
    resp, annot, gname = np.full(len(adf), np.nan), np.array([None]*len(adf), object), \
        np.array([None]*len(adf), object)
    for i, (c, m, okl) in enumerate(zip(adf["chrom_hg19"], adf["mid_hg38"], adf["lifted"])):
        if not okl or c not in wi_by_c:
            continue
        S, E, V, A, G = wi_by_c[c]
        j = np.searchsorted(S, m, side="right") - 1
        if 0 <= j < len(S) and S[j] <= m < E[j]:
            resp[i], annot[i], gname[i] = V[j], A[j], G[j]
    adf["wi38_response"] = resp
    adf["wi38_annot"] = annot
    adf["wi38_gene"] = gname
    nmatch = int(np.isfinite(resp).sum())
    R(f"  anchors matched into a WI-38 peak (lifted midpoint inside peak): "
      f"{nmatch:,}/{len(adf):,} ({100*nmatch/len(adf):.1f}%)")
    stats["anchors_matched_wi38"] = nmatch

    # ---------------- bigWig quantification ----------------
    R("\n[4] SIGNAL QUANTIFICATION (bigWig mean coverage over hg19 anchors)\n" + "-" * 84)
    import pyBigWig
    chroms = adf["chrom_hg19"].to_numpy()
    starts = adf["start_hg19"].to_numpy()
    ends = adf["end_hg19"].to_numpy()
    sig = {}

    def quantify(path, label):
        bw = pyBigWig.open(str(path))
        avail = set(bw.chroms().keys())
        pref = "" if "chr1" in avail else ("chr" if "1" in avail else None)
        vals = np.full(len(adf), np.nan)
        for i in range(len(adf)):
            c = chroms[i]
            cc = c if c in avail else (c[3:] if c[3:] in avail else None)
            if cc is None:
                continue
            try:
                v = bw.stats(cc, int(starts[i]), int(ends[i]), type="mean")[0]
            except Exception:
                v = None
            if v is not None:
                vals[i] = v
        bw.close()
        n_ok = int(np.isfinite(vals).sum())
        R(f"  {label:34s} covered {n_ok:,}/{len(adf):,} anchors  "
          f"mean {np.nanmean(vals):.4f}")
        return vals

    for gsm, fn in GM21_ATAC_BW.items():
        p = DATA / "GSE206402" / "bw" / fn
        if p.exists():
            sig[f"GM21_ATAC::{fn.split('_',1)[1].replace('.bw','')}"] = quantify(p, fn)
        else:
            R(f"  MISSING {fn}")
    for gsm, fn in GM21_K27_BW.items():
        p = DATA / "GSE205898" / "bw" / fn
        if p.exists():
            sig[f"GM21_K27ac::{fn.split('_',1)[1].replace('_dedup_blacklisted_merge.bw','')}"] = quantify(p, fn)
        else:
            R(f"  MISSING {fn}")
    for arm, fns in IMR90_K27_BW.items():
        for fn in fns:
            p = DATA / "GSE74238" / "bw" / fn
            if p.exists():
                sig[f"IMR90_K27ac::{arm}::{fn.split('_',1)[1].replace('.bigWig','')}"] = quantify(p, fn)
            else:
                R(f"  MISSING {fn}")

    sm = pd.DataFrame(sig)
    adf.to_csv(OUT / "anchors.csv", index=False)
    sm.to_csv(OUT / "signal_matrix.csv", index=False)
    json.dump(stats, open(OUT / "harmonization_stats.json", "w"), indent=2)
    (OUT / "STEP1_LOG.txt").write_text("\n".join(LOG))
    R(f"\n  wrote anchors.csv ({adf.shape}) and signal_matrix.csv ({sm.shape})")
    R("STEP1_DONE")


if __name__ == "__main__":
    main()
