#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 23_analyze_encode_imr90.py
#
# Audit the ENCODE IMR-90 ATAC experiment and its ChromBPNet annotation.
#
# CONSUMES: data/ENCSR200OML/experiment.json, data/ENCSR978WIX/annotation.json
# PRODUCES: output/tables/encsr*_*.csv, cdkn2a_cdkn1a_overlap.csv
# ---------------------------------------------------------------------------
"""
Parse ENCSR200OML (IMR-90 ATAC-seq experiment) and ENCSR978WIX (ChromBPNet
annotation derived from it) JSON records pulled by 23_fetch_encode_imr90.sh,
and analyze the downloaded processed peak/region files.

Outputs (all under output/):
  tables/encsr200oml_file_inventory.csv       - every file on the experiment
  tables/encsr200oml_replicate_qc.csv         - per-replicate depth/QC metrics
  tables/encsr200oml_audit_flags.csv          - ENCODE's own audit warnings
  tables/encsr978wix_file_inventory.csv       - every file on the annotation
  tables/encsr200oml_peak_summary.csv         - peak count + width percentiles
                                                 per downloaded peak set
  tables/cdkn2a_cdkn1a_overlap.csv            - peak overlap at the two loci
  figures/encsr200oml_peak_width_distribution.png

No network access in this script -- it only reads files fetched by
23_fetch_encode_imr90.sh.
"""
import csv
import gzip
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D1 = ROOT / "data" / "ENCSR200OML"
D2 = ROOT / "data" / "ENCSR978WIX"
TABLES = ROOT / "output" / "tables"
FIGURES = ROOT / "output" / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. ENCSR200OML file inventory
# ---------------------------------------------------------------------------
exp = json.loads((D1 / "experiment.json").read_text())

with open(TABLES / "encsr200oml_file_inventory.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["accession", "file_type", "output_type", "file_format",
                "assembly", "file_size_bytes", "status", "biological_replicates",
                "preferred_default", "analysis"])
    for f in exp["files"]:
        analyses = f.get("analyses")
        analysis_title = None
        if isinstance(analyses, list) and analyses:
            a0 = analyses[0]
            analysis_title = a0.get("title") if isinstance(a0, dict) else a0
        w.writerow([
            f.get("accession"), f.get("file_type"), f.get("output_type"),
            f.get("file_format"), f.get("assembly"), f.get("file_size"),
            f.get("status"), ";".join(str(x) for x in f.get("biological_replicates", [])),
            f.get("preferred_default"), analysis_title,
        ])
print("wrote encsr200oml_file_inventory.csv "
      f"({len(exp['files'])} released files; {len(exp.get('revoked_files', []))} revoked, "
      "not included)")

# ---------------------------------------------------------------------------
# 2. Per-replicate depth / QC metrics
# ---------------------------------------------------------------------------
files_by_acc = {f["accession"]: f for f in exp["files"]}

# fastq -> raw read pairs, per bio replicate
fastq_reads = {}  # bio_rep -> read_count (per mate; paired-end so this = read pairs)
for f in exp["files"]:
    if f.get("file_type") == "fastq" and f.get("read_count"):
        reps = f.get("biological_replicates", [])
        if reps:
            fastq_reads.setdefault(reps[0], []).append(f["read_count"])

replicate_rows = []
for r in exp["replicates"]:
    bio = r.get("biological_replicate_number")
    tech = r.get("technical_replicate_number")
    lib = r.get("library", {})
    bs = lib.get("biosample", {}) if isinstance(lib, dict) else {}
    frac = bs.get("subcellular_fraction_term_name", "unspecified")
    raw_pairs = fastq_reads.get(bio, [None])[0]

    # find the "unfiltered alignments" and "alignments" (filtered) bam for this rep
    unfiltered = next((f for f in exp["files"]
                        if f.get("output_type") == "unfiltered alignments"
                        and bio in f.get("biological_replicates", [])), None)
    filtered = next((f for f in exp["files"]
                      if f.get("output_type") == "alignments"
                      and bio in f.get("biological_replicates", [])), None)

    pct_mapped = usable_fragments = nrf = pbc1 = pbc2 = tss_enrich = frac_mito = None
    if unfiltered:
        for qm in unfiltered.get("quality_metrics", []):
            if "pct_mapped_reads" in qm:
                pct_mapped = qm["pct_mapped_reads"]
                frac_mito = qm.get("frac_mito_reads")
    if filtered:
        for qm in filtered.get("quality_metrics", []):
            if "usable_fragments" in qm:
                usable_fragments = qm["usable_fragments"]
            if "NRF" in qm:
                nrf, pbc1, pbc2 = qm["NRF"], qm["PBC1"], qm["PBC2"]
            if "tss_enrichment" in qm:
                tss_enrich = qm["tss_enrichment"]

    replicate_rows.append({
        "biological_replicate": bio,
        "technical_replicate": tech,
        "biosample_accession": bs.get("accession"),
        "subcellular_fraction": frac,
        "raw_read_pairs": raw_pairs,
        "pct_mapped_reads": pct_mapped,
        "frac_mito_reads": frac_mito,
        "usable_fragments_final": usable_fragments,
        "NRF": nrf,
        "PBC1": pbc1,
        "PBC2": pbc2,
        "TSS_enrichment": tss_enrich,
    })

with open(TABLES / "encsr200oml_replicate_qc.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(replicate_rows[0].keys()))
    w.writeheader()
    w.writerows(replicate_rows)
print("wrote encsr200oml_replicate_qc.csv")
for row in replicate_rows:
    print(" ", row)

# ---------------------------------------------------------------------------
# 3. ENCODE's own audit warnings (verbatim)
# ---------------------------------------------------------------------------
with open(TABLES / "encsr200oml_audit_flags.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["level", "category", "detail"])
    for level, items in exp.get("audit", {}).items():
        for it in items:
            w.writerow([level, it.get("category"), it.get("detail")])
print("wrote encsr200oml_audit_flags.csv")

# ---------------------------------------------------------------------------
# 4. ENCSR978WIX file inventory
# ---------------------------------------------------------------------------
annot = json.loads((D2 / "annotation.json").read_text())
with open(TABLES / "encsr978wix_file_inventory.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["accession", "file_type", "output_type", "file_format",
                "file_format_type", "file_size_bytes", "status",
                "output_category", "pipeline_step_title"])
    for f in annot["files"]:
        asv = f.get("analysis_step_version")
        step_title = None
        if isinstance(asv, dict):
            step = asv.get("analysis_step", {})
            step_title = step.get("title") if isinstance(step, dict) else None
        w.writerow([
            f.get("accession"), f.get("file_type"), f.get("output_type"),
            f.get("file_format"), f.get("file_format_type"), f.get("file_size"),
            f.get("status"), f.get("output_category"), step_title,
        ])
print(f"wrote encsr978wix_file_inventory.csv ({len(annot['files'])} files)")

# ---------------------------------------------------------------------------
# 5. Peak set summaries (count + width percentiles)
# ---------------------------------------------------------------------------
PEAK_SETS = {
    "pseudoreplicated_peaks_DEFAULT": (
        D1 / "ENCFF243NTP_pseudoreplicated_peaks.bed.gz",
        "ENCFF243NTP", "pseudoreplicated peaks (preferred_default=True)"),
    "conservative_idr_peaks": (
        D1 / "ENCFF982UNH_conservative_idr_peaks.bed.gz",
        "ENCFF982UNH", "conservative IDR thresholded peaks"),
    "idr_thresholded_peaks": (
        D1 / "ENCFF114GDS_idr_thresholded_peaks.bed.gz",
        "ENCFF114GDS", "IDR thresholded peaks (pooled pseudoreplicate IDR)"),
    "chrombpnet_selected_regions": (
        D2 / "ENCFF815MUM_selected_regions.bed.gz",
        "ENCFF815MUM", "ChromBPNet selected regions (bed3, fixed-width)"),
}


def load_bed_intervals(path):
    intervals = []  # (chrom, start, end)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals.append((chrom, start, end))
    return intervals


def percentile(sorted_vals, p):
    if not sorted_vals:
        return None
    k = (len(sorted_vals) - 1) * p / 100.0
    f = int(k)
    c = min(f + 1, len(sorted_vals) - 1)
    if f == c:
        return sorted_vals[f]
    return sorted_vals[f] + (sorted_vals[c] - sorted_vals[f]) * (k - f)


peak_data = {}         # label -> deduplicated interval list (used for all downstream analysis)
summary_rows = []
for label, (path, acc, desc) in PEAK_SETS.items():
    if not path.exists():
        print(f"SKIP {label}: file not found ({path})")
        continue
    raw_intervals = load_bed_intervals(path)
    n_raw = len(raw_intervals)
    # IMPORTANT: this ENCODE4 ATAC pipeline output format is NOT deduplicated
    # by genomic interval -- the same (chrom,start,end) footprint can appear
    # as multiple rows (different Peak_ID/summit/score from different
    # replicate/pseudoreplicate support comparisons that both cleared the
    # peak-calling threshold for that locus). A raw line count therefore
    # overstates the number of distinct accessible loci. We dedupe by exact
    # (chrom,start,end) and use that for peak counts / width stats / the
    # CDKN2A-CDKN1A overlap check below.
    intervals = sorted(set(raw_intervals))
    n_dedup = len(intervals)
    widths = sorted(e - s for (_, s, e) in intervals)
    peak_data[label] = intervals
    row = {
        "peak_set": label,
        "file_accession": acc,
        "description": desc,
        "n_lines_raw": n_raw,
        "n_distinct_intervals": n_dedup,
        "duplication_ratio": round(n_raw / n_dedup, 3) if n_dedup else None,
        "width_min": widths[0] if widths else None,
        "width_p10": percentile(widths, 10),
        "width_p25": percentile(widths, 25),
        "width_median": percentile(widths, 50),
        "width_mean": round(statistics.mean(widths), 1) if widths else None,
        "width_p75": percentile(widths, 75),
        "width_p90": percentile(widths, 90),
        "width_p99": percentile(widths, 99),
        "width_max": widths[-1] if widths else None,
    }
    summary_rows.append(row)
    print(f"{label}: n_raw_lines={n_raw} n_distinct={n_dedup} "
          f"(dup ratio {row['duplication_ratio']}) widths min/median/max = "
          f"{row['width_min']}/{row['width_median']}/{row['width_max']}")

with open(TABLES / "encsr200oml_peak_summary.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
    w.writeheader()
    w.writerows(summary_rows)
print("wrote encsr200oml_peak_summary.csv")

# ---------------------------------------------------------------------------
# 6. CDKN2A / CDKN1A promoter overlap check
# ---------------------------------------------------------------------------
# Gene-body coordinates as given in the task brief (GRCh38):
#   CDKN2A: chr9:21,967,752-21,995,043  (minus strand; TSS of the canonical
#           p16-INK4a transcript is at the HIGH-coordinate end of this span,
#           ~21,994,490-21,995,043; this locus also has an alternate ARF
#           promoter further downstream/high-coordinate. Strand orientation
#           is analyst knowledge, not sourced from the ENCODE JSON.)
#   CDKN1A: chr6:36,676,460-36,685,703  (plus strand; TSS at the LOW-coordinate
#           end, ~36,676,460.)
# We pad the given span +/-2kb on both sides and report ANY peak overlap
# inside that padded window ("rough overlap check", per task instructions),
# plus separately flag whether the overlap falls specifically in the last
# 2kb proximal-promoter slice nearest the TSS.
LOCI = {
    "CDKN2A_p16": {"chrom": "chr9", "start": 21_967_752, "end": 21_995_043,
                   "strand": "-", "tss_end": "high"},
    "CDKN1A_p21": {"chrom": "chr6", "start": 36_676_460, "end": 36_685_703,
                   "strand": "+", "tss_end": "low"},
}
PAD = 2000

overlap_rows = []
for gene, loc in LOCI.items():
    win_start = loc["start"] - PAD
    win_end = loc["end"] + PAD
    if loc["tss_end"] == "high":
        promoter_start, promoter_end = loc["end"] - 2000, loc["end"] + PAD
    else:
        promoter_start, promoter_end = loc["start"] - PAD, loc["start"] + 2000

    for peak_set, intervals in peak_data.items():
        hits = [(c, s, e) for (c, s, e) in intervals
                if c == loc["chrom"] and s < win_end and e > win_start]
        promoter_hits = [(c, s, e) for (c, s, e) in hits
                          if s < promoter_end and e > promoter_start]
        overlap_rows.append({
            "gene": gene,
            "peak_set": peak_set,
            "window_chrom": loc["chrom"],
            "window_start_padded": win_start,
            "window_end_padded": win_end,
            "n_peaks_in_window": len(hits),
            "peaks_in_window": ";".join(f"{c}:{s}-{e}" for c, s, e in hits),
            "n_peaks_in_proximal_promoter_2kb": len(promoter_hits),
            "peaks_in_proximal_promoter": ";".join(f"{c}:{s}-{e}" for c, s, e in promoter_hits),
        })
        print(f"{gene} / {peak_set}: {len(hits)} peak(s) in padded gene window, "
              f"{len(promoter_hits)} in proximal-promoter slice")

with open(TABLES / "cdkn2a_cdkn1a_overlap.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(overlap_rows[0].keys()))
    w.writeheader()
    w.writerows(overlap_rows)
print("wrote cdkn2a_cdkn1a_overlap.csv")

# ---------------------------------------------------------------------------
# 7. Figure: peak width distribution (default peak set + conservative IDR)
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
colors = {"pseudoreplicated_peaks_DEFAULT": "#2b6cb0",
          "conservative_idr_peaks": "#c05621"}
labels = {"pseudoreplicated_peaks_DEFAULT": "Default peak set (ENCFF243NTP, n={n})",
          "conservative_idr_peaks": "Conservative IDR (ENCFF982UNH, n={n})"}
bins = np.logspace(np.log10(50), np.log10(20000), 60)
for key in ["pseudoreplicated_peaks_DEFAULT", "conservative_idr_peaks"]:
    if key not in peak_data:
        continue
    widths = [e - s for (_, s, e) in peak_data[key]]
    ax.hist(widths, bins=bins, alpha=0.55, color=colors[key],
            label=labels[key].format(n=len(widths)))
ax.set_xscale("log")
ax.set_xlabel("Peak width (bp, log scale)")
ax.set_ylabel("Number of peaks")
ax.set_title("ENCSR200OML (IMR-90 ATAC-seq) peak width distribution")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES / "encsr200oml_peak_width_distribution.png")
print("wrote figures/encsr200oml_peak_width_distribution.png")

print("\nDone.")
