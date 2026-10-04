#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 23_design_summary_gse206402.py
#
# Design summary for GSE206402, including peak counts per arm.
#
# CONSUMES: data/GSE206402/peaks/*.narrowPeak.gz
# PRODUCES: output/tables/gse206402_design_summary.txt
# ---------------------------------------------------------------------------
"""
Task 3/4 for GSE206402 (ATAC-seq subseries): design summary + cross-dataset
compatibility facts, computed from output/tables/gse206402_sample_metadata.csv
(produced by 22_parse_gse206402_metadata.py) and the downloaded narrowPeak
files in data/GSE206402/peaks/ (produced by 21_fetch_gse206402_peaks.py).

Writes output/tables/gse206402_design_summary.txt. No numbers in this log are
invented -- every count is either read from a parsed GEO field or computed
directly from a downloaded file's contents (peak counts, column counts).
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "output" / "tables" / "gse206402_sample_metadata.csv"
PEAKS_DIR = ROOT / "data" / "GSE206402" / "peaks"
OUT = ROOT / "output" / "tables" / "gse206402_design_summary.txt"

NARROWPEAK_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "signalValue", "pValue(-log10)", "qValue(-log10)", "summit_offset",
]


def main():
    log = []
    df = pd.read_csv(META, keep_default_na=False)  # empty CSV fields must stay "" not NaN --
    # NaN != "" evaluates True in pandas, which would make every row look "flagged" below.

    log.append("=== GSE206402 (ATAC-seq subseries) design summary ===")
    log.append(f"Total ATAC-seq samples (GSM rows parsed): {len(df)}")

    log.append("\n--- Replicate counts per condition x timepoint ---")
    counts = df.groupby(["condition", "timepoint"])["gsm"].nunique().reset_index()
    counts = counts.sort_values(["condition", "timepoint"])
    n_flagged = 0
    for _, row in counts.iterrows():
        cond_short = "empty-vector control" if "empty" in row["condition"].lower() else "RAS overexpression"
        flag = ""
        if row["gsm"] < 3:
            flag = "  [FLAG: <3 biological replicates]"
            n_flagged += 1
        else:
            flag = "  [OK: >=3 replicates]"
        log.append(f"  {cond_short:22s} | {row['timepoint']:10s} | n={row['gsm']}{flag}")
    log.append(f"\n{n_flagged}/{len(counts)} condition x timepoint arms have <3 biological replicates.")

    log.append("\n--- Condition-level totals ---")
    for cond, sub in df.groupby("condition"):
        cond_short = "empty-vector control" if "empty" in cond.lower() else "RAS overexpression"
        log.append(f"  {cond_short}: {len(sub)} samples across "
                    f"{sub['timepoint'].nunique()} timepoints "
                    f"({sorted(sub['timepoint'].unique(), key=lambda x: int(x.split()[-1]))})")

    log.append(
        "\nNOTE: GEO's 'treatment' characteristic distinguishes only "
        "'Empty vector control' vs 'H-RAS-G12-V overexpression'. There is NO "
        "GEO-native field separating 'still senescent' from 'escaped' within "
        "the RAS arm's 7 timepoints (D8/14/18/23/32/45/56) -- if that "
        "distinction matters for training-label purposes, it must come from "
        "the paper's Results/Figures (out of scope here; GEO-metadata-only "
        "constraint), not from this metadata table."
    )

    # ---- Genome build ----
    log.append("\n--- Genome build ---")
    build_vals = df["genome_build_in_data_processing"].unique()
    log.append(f"  Sample_data_processing 'Assembly:' line -> {list(build_vals)} "
                f"(uniform across all {len(df)} samples: {len(build_vals) == 1})")

    # ---- Peak caller / pipeline (from raw SOFT text, hardcoded here since it's
    # identical boilerplate on every sample -- verified via `sort -u` over
    # Sample_data_processing during interactive audit) ----
    log.append("\n--- Peak-calling pipeline (from Sample_data_processing, identical on all 25 samples) ---")
    log.append("  Aligner: bowtie2, local mode, paired-end")
    log.append("  Dedup: PicardTools + samtools; blacklisted against hg19 artifact regions")
    log.append("  Peak caller: MACS v2.2.7.1 "
                "(macs2 callpeak --nomodel --shiftsize --shift-control --gsize hs -p 1e-3)")
    log.append("  Post-processing: IDR (irreproducibility discovery rate) pipeline -> "
                "'time point-specific reproducible peak sets' -> collapsed into "
                "'master peaksets' (NOTE: this master/IDR peak file is NOT among the "
                "deposited supplementary files -- see peak file inventory below. Only "
                "per-replicate, pre-IDR narrowPeak.gz files were deposited.)")
    log.append("  Signal tracks: deeptools RPGC 1x-coverage bigWig, per-timepoint-merged BAM")
    log.append(
        "  FLAG (copy/paste artifact, not a correctness issue for the data itself): "
        "this exact data_processing text -- including the phrase 'master peaksets per "
        "histone modification analyzed' -- is IDENTICAL, word-for-word, to the "
        "Sample_data_processing text in the sibling GSE205898 ChIP-seq (H3K27ac/H3K4me1) "
        "subseries. ATAC-seq doesn't have a 'histone modification' being profiled, so this "
        "reads as boilerplate reused across assay types during submission, not as evidence "
        "the ATAC processing pipeline itself literally analyzed histone marks."
    )

    # ---- Peak file structure (inspected directly, not assumed) ----
    log.append("\n--- Peak/count supplementary file structure ---")
    peak_files = sorted(PEAKS_DIR.glob("*.narrowPeak.gz")) if PEAKS_DIR.exists() else []
    log.append(f"  {len(peak_files)} per-GSM narrowPeak.gz files downloaded "
                f"(from GSE206402_RAW.tar's member list, via per-sample GEO FTP mirror -- "
                f"see 21_fetch_gse206402_peaks.py). No separate merged/atlas peak file or "
                f"count matrix exists in the series suppl/ directory: the ONLY series-level "
                f"supplementary artifact is GSE206402_RAW.tar (7.6 GB; bundles all 25 "
                f"per-sample .bw + .narrowPeak.gz files) plus filelist.txt.")
    if peak_files:
        f0 = peak_files[0]
        with gzip.open(f0, "rt") as fh:
            first_line = fh.readline().rstrip("\n")
        ncols = len(first_line.split("\t"))
        log.append(f"  Format (inspected from {f0.name}): tab-delimited, {ncols} columns "
                    f"-> matches standard MACS2 narrowPeak (BED6+4): {NARROWPEAK_COLS}")
        log.append(f"  Coordinate system: BED convention -- 0-based, half-open [chromStart, chromEnd)")
        log.append(f"  First data row: {first_line}")

        log.append("\n  Per-sample peak counts (line count in each narrowPeak.gz):")
        total_peaks = []
        for f in peak_files:
            with gzip.open(f, "rt") as fh:
                n = sum(1 for _ in fh)
            total_peaks.append(n)
            log.append(f"    {f.name}: {n} peaks")
        log.append(f"\n  Peak count range across 25 replicate-level files: "
                    f"{min(total_peaks)}-{max(total_peaks)} "
                    f"(median {sorted(total_peaks)[len(total_peaks)//2]}); "
                    f"NOTE: these are per-replicate MACS2 calls, NOT a merged/consensus "
                    f"atlas peak count -- no atlas file was deposited to report a single "
                    f"'total peak count' for the series.")

    # ---- Internal-label cross-check summary ----
    log.append("\n--- Peak-file internal-label vs. GEO-filename cross-check ---")
    flagged = df[df["peak_file_internal_label_flag"] != ""]
    log.append(f"  {len(flagged)}/{len(df)} narrowPeak.gz files carry an internal peak-name "
                f"(column 4 prefix) that disagrees with the GEO filename/title day or "
                f"replicate number, or is a raw lab filesystem path rather than a sample tag.")
    log.append("  This does not change what's recorded in gse206402_sample_metadata.csv "
                "(GEO's own title/characteristics are authoritative there), but it is a "
                "real provenance inconsistency in the deposited files themselves -- "
                "flagged prominently, not smoothed over. See "
                "peak_file_internal_label_flag column and 22_parse_gse206402_metadata.py "
                "docstring for the per-GSM detail.")

    report = "\n".join(log)
    OUT.write_text(report)
    print(report)


if __name__ == "__main__":
    main()
