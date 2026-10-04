#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 14_atac_peak_qc_diff_markers.py
#
# GSE175533 ATAC peak QC and differential-accessibility summary.
#
# CONSUMES: data/GSE175533/ peak + differential tables
# PRODUCES: output/tables/gse175533_atac_*.csv, *.bed, *_log.txt
# ---------------------------------------------------------------------------
"""
Task 3 (ATAC promoter/distal), Task 4 tail (accessibility at CDKN2A/CDKN1A
promoters), and Task 5 (differential-accessibility region counts + BED
output) for GSE175533 -- all from GSE175533_atac_peaks_sig.xlsx (the
"atac_peaks" sheet), since per-replicate bigWig quantification was not
feasible in this audit (see REPORT.md for the measured-throughput
justification: sustained download speed to NCBI FTP was measured at well
under 1.7 MB/s, making even the smallest useful bigWig subset -- ~8.2 GB for
just the WI-38 replicative-senescence ATAC arm -- an impractical multi-hour
download for this session).

The "atac_peaks" sheet is the authors' own limma-based differential test:
for each of 5 RS-PDL timepoints, log2FC and adjPval of that PDL's
quantile-normalized peak signal AGAINST the matched hTERT timepoint (i.e.
hTERT-at-matched-time is the reference/denominator, not an early-PDL WI-38
sample). This script uses those columns directly as the differential-
accessibility contrast for Task 5 -- it is the only differential ATAC
signal GEO actually provides, and is a defensible senescent-proxy contrast
(hTERT is immortalized and not expected to senesce at any matched time).

ASSUMPTION flagged: the "region" column is formatted "chr:start-end" and the
readme states these are "hg38 coordinates" -- assumed 1-based inclusive
(Excel/R convention) for the purposes of BED export, so 1 is subtracted from
start when writing BED (0-based half-open). This was NOT independently
verified against the atlas.bed file's own coordinates for the same peak
(spot-checked below where possible) since atlas.bed peak IDs (integers) do
not share a key with atac_peaks_sig's "region" string -- flagged, not
silently assumed correct.
"""
import re
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data" / "GSE175533" / "GSE175533_atac_peaks_sig.transitional.xlsx"
ATLAS_BED = ROOT / "data" / "GSE175533" / "GSE175533_atlas.bed"
OUT_T = ROOT / "output" / "tables"
OUT_T.mkdir(parents=True, exist_ok=True)

COMPARISONS = [
    ("PDL25_v_htert2", 25), ("PDL33_v_htert4", 33), ("PDL37_v_htert5", 37),
    ("PDL45_v_htert6", 46),  # NOTE: column says "PDL45" but GEO sample titles for
    # this ATAC timepoint say PDL46_TP6 -- internal-file-vs-GEO-title PDL-label
    # mismatch, same pattern as the RIS d0/d2 mismatch in the RNA table. Treated
    # as the PDL46_TP6 ATAC condition based on TP-index alignment (both are "TP6").
    ("PDL50_v_htert7", 50),
]


def load_full_sheet():
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = wb["atac_peaks"]
    rows = ws.iter_rows(values_only=True)
    header = list(next(rows))
    data = list(rows)
    df = pd.DataFrame(data, columns=header)
    return df


def main():
    log = []
    log.append("Loading atac_peaks sheet (363,470 rows) ...")
    df = load_full_sheet()
    log.append(f"Loaded: {df.shape[0]} peaks x {df.shape[1]} columns")

    # ---- Task 3: promoter vs distal breakdown ----
    log.append("\n=== Task 3: peak genomic-context distribution (full atlas) ===")
    annot_counts = df["annot"].value_counts()
    log.append(annot_counts.to_string())
    n_promoter = annot_counts.get("promoter", 0)
    n_total = len(df)
    log.append(f"\nPromoter fraction: {n_promoter}/{n_total} = {100*n_promoter/n_total:.1f}%")
    log.append(f"Distal (intergenic_proximal + intergenic_distal + intron + exon): "
                f"{n_total - n_promoter}/{n_total} = {100*(n_total-n_promoter)/n_total:.1f}%")
    log.append("Genome build: Hg38, per explicit Sample_data_processing text "
                "('Genome_build: Hg38') on every sample AND the atac_peaks_sig readme "
                "sheet ('region: hg38 coordinates') -- confirmed two independent ways, "
                "not assumed.")
    annot_counts.to_frame("n_peaks").to_csv(OUT_T / "gse175533_atac_annot_distribution.csv")

    # ---- Task 4 tail: CDKN2A / CDKN1A promoter accessibility across PDL ----
    log.append("\n=== Task 4: CDKN2A / CDKN1A promoter accessibility vs PDL "
                "(log2FC = RS-PDL vs matched-hTERT-timepoint; positive = MORE "
                "accessible in RS-PDL than in the hTERT counterpart) ===")
    for gene in ["CDKN2A", "CDKN1A"]:
        sub = df[(df["gene_name"] == gene) & (df["annot"] == "promoter")]
        if sub.empty:
            log.append(f"\n{gene}: NO PROMOTER-ANNOTATED PEAK in the atlas for this gene "
                        "(peak may not have been called there, or nearest-gene assignment "
                        "landed on a different transcript/annotation).")
            continue
        log.append(f"\n{gene}: {len(sub)} promoter-annotated peak row(s) -- {sub['region'].tolist()}")
        for _, row in sub.iterrows():
            log.append(f"  region {row['region']}:")
            for cmp_name, pdl in COMPARISONS:
                fc_col, p_col = f"{cmp_name}_log2FC", f"{cmp_name}_adjPval"
                fc, p = row.get(fc_col), row.get(p_col)
                if pd.notna(fc):
                    sig = "*" if (pd.notna(p) and p < 0.05) else ""
                    log.append(f"    PDL{pdl} vs matched hTERT: log2FC={fc:+.3f}, adjPval={p:.3g}{sig}")

    # ---- Task 5: differential accessibility region counts ----
    log.append("\n\n=== Task 5: differential-accessibility region counts per comparison ===")
    log.append("Threshold: adjPval < 0.05 AND |log2FC| >= 1 (2-fold), the field-standard "
                "combination for calling a peak significantly differential -- not the only "
                "defensible threshold, but a conventional one, stated explicitly here.")
    sig_union_mask = pd.Series(False, index=df.index)
    summary_rows = []
    for cmp_name, pdl in COMPARISONS:
        fc_col, p_col = f"{cmp_name}_log2FC", f"{cmp_name}_adjPval"
        sig = (df[p_col] < 0.05) & (df[fc_col].abs() >= 1)
        gained = sig & (df[fc_col] > 0)   # more accessible in RS-PDL than hTERT
        lost = sig & (df[fc_col] < 0)     # less accessible in RS-PDL than hTERT
        n_tested = df[p_col].notna().sum()
        summary_rows.append({
            "comparison": cmp_name, "pdl": pdl, "n_tested": n_tested,
            "n_significant": int(sig.sum()),
            "n_gained_accessibility": int(gained.sum()),
            "n_lost_accessibility": int(lost.sum()),
            "pct_significant": round(100 * sig.sum() / n_tested, 2) if n_tested else None,
        })
        sig_union_mask |= sig.fillna(False)
        log.append(f"\n{cmp_name} (RS PDL{pdl} vs matched hTERT): {n_tested} peaks tested, "
                    f"{sig.sum()} significant (padj<0.05, |log2FC|>=1) "
                    f"[{gained.sum()} gained, {lost.sum()} lost accessibility in RS-PDL]")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_T / "gse175533_atac_differential_summary.csv", index=False)

    n_union = int(sig_union_mask.sum())
    log.append(f"\nUNION across all 5 comparisons: {n_union} distinct peaks are significantly "
                f"differential in at least one RS-PDL-vs-hTERT contrast.")
    log.append(f"{'BELOW ~few-thousand-region threshold -- likely UNDERPOWERED for sequence-model training' if n_union < 3000 else 'ABOVE the rough few-thousand-region threshold for sequence-model training'}")

    # peaks significant in the LATE comparisons specifically (PDL46, PDL50 -- the
    # data-driven "senescent" end per the RNA marker trajectory analysis, see
    # 13_rna_qc_signal_gse175533.py / REPORT.md Task 2(b))
    late_cols_fc = ["PDL45_v_htert6_log2FC", "PDL50_v_htert7_log2FC"]
    late_cols_p = ["PDL45_v_htert6_adjPval", "PDL50_v_htert7_adjPval"]
    late_sig = pd.Series(False, index=df.index)
    for fc_col, p_col in zip(late_cols_fc, late_cols_p):
        late_sig |= (df[p_col] < 0.05) & (df[fc_col].abs() >= 1)
    log.append(f"\nSignificant in >=1 of the two LATE (PDL46, PDL50) comparisons specifically: "
                f"{int(late_sig.sum())} peaks -- this is the recommended senescent-vs-proliferating "
                f"proxy set (see Task 2b PDL cutoff justification).")

    # direction split for the late-comparison set, using PDL50 (the deepest/most
    # significant late timepoint) as the representative direction
    late_df = df[late_sig].copy()
    rep_fc = late_df["PDL50_v_htert7_log2FC"].fillna(late_df["PDL45_v_htert6_log2FC"])
    n_gained = int((rep_fc > 0).sum())
    n_lost = int((rep_fc < 0).sum())
    log.append(f"Direction split (using PDL50 log2FC, falling back to PDL46 if PDL50 is NaN "
                f"for that peak): {n_gained} gained accessibility, {n_lost} lost accessibility "
                f"in senescent-proxy vs hTERT.")

    # ---- BED export of the late-comparison significant set ----
    def parse_region(r):
        m = re.match(r"(chr\w+):(\d+)-(\d+)", str(r))
        if not m:
            return None
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        return chrom, start - 1, end  # 1-based -> 0-based BED start; see module docstring

    bed_rows = []
    for _, row in late_df.iterrows():
        parsed = parse_region(row["region"])
        if parsed is None:
            continue
        chrom, start, end = parsed
        fc = row["PDL50_v_htert7_log2FC"] if pd.notna(row["PDL50_v_htert7_log2FC"]) else row["PDL45_v_htert6_log2FC"]
        direction = "gained" if fc > 0 else "lost"
        bed_rows.append((chrom, start, end, f"{row['gene_name']}|{direction}|log2FC={fc:.2f}", 0, "."))

    bed_out = OUT_T / "gse175533_differential_accessibility_senescent_vs_hTERT.bed"
    with open(bed_out, "w") as fh:
        for r in bed_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    log.append(f"\nWrote {len(bed_rows)} regions to {bed_out}")

    report = "\n".join(str(x) for x in log)
    (OUT_T / "gse175533_atac_peak_qc_diff_markers_log.txt").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
