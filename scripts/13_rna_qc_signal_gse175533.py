#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 13_rna_qc_signal_gse175533.py
#
# GSE175533 RNA QC, PCA and senescence-marker check.
#
# CONSUMES: data/GSE175533/ TPM tables
# PRODUCES: output/tables/gse175533_rna_*.csv, gse175533_rna_qc_signal_log.txt
# ---------------------------------------------------------------------------
"""
Task 3 (RNA QC) + Task 4 (RNA signal check + marker panel) for GSE175533.

Source: GSE175533_TPM_table.transitional.xlsx (Strict-OOXML-fixed copy of the
deposited GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx -- see 12_fix_strict_ooxml.sh).

ASSUMPTIONS flagged explicitly:
  1. Only Salmon TPM (already per-million-normalized) is deposited for bulk
     RNA -- no raw/pre-normalization counts. "Library size" in the literal
     raw-read-count sense (Task 3) is therefore NOT computable from what GEO
     provides; used instead: (a) each sample's TPM column sum as a sanity
     check (should be ~1e6 by definition of TPM -- a real deviation would
     flag a broken column), and (b) detected-gene count (TPM>0), which IS a
     meaningful QC metric independent of that limitation.
  2. Task 4 asks for "CPM-normalize, log-transform" -- since only TPM (not
     raw counts) exists, this script log2(TPM+1)-transforms directly rather
     than re-deriving CPM from TPM (which would be circular/meaningless).
     TPM and CPM are both per-million normalizations; log2(TPM+1) is the
     standard substitute when only TPM is available. Flagged as a deliberate
     adaptation to available data, not a silent substitution.
  3. RS_TC_raw_TPM (not RS_TC_batch_corrected_TPM) is used for this script's
     independent QC/PCA -- the workbook's own readme states the
     batch-corrected sheet was "used only for visualization *Figure 1E"; an
     independent audit should run PCA on the uncorrected data to see the
     actual raw structure (including any batch effects) rather than
     inheriting the authors' correction.
  4. RIS_TPM's "d0_A/B/C" columns are assumed to be the same 3 samples GEO
     titles this series' record as "RIS_d2_no_xray_A/B/C" (the mock-IR
     control) -- day-numbering differs between the GEO sample titles
     (absolute culture day) and this internal table (day relative to
     IR treatment, so the pre-treatment baseline is "day 0"). Both give
     exactly 3 samples occupying the same "control before the IR time
     course" structural position, which is why this mapping is treated as
     near-certain rather than a guess -- but it is NOT independently
     confirmed via a shared sample identifier, so it is flagged here rather
     than asserted silently.
"""
import re
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "data" / "GSE175533" / "GSE175533_TPM_table.transitional.xlsx"
OUT_T = ROOT / "output" / "tables"
OUT_F = ROOT / "output" / "figures"
OUT_T.mkdir(parents=True, exist_ok=True)
OUT_F.mkdir(parents=True, exist_ok=True)

MARKERS = {
    "CDKN2A": "up (p16)", "CDKN1A": "up (p21)", "LMNB1": "down",
    "MKI67": "down", "IL6": "up (SASP)", "CXCL8": "up (SASP, aka IL8)",
    "SERPINE1": "up (SASP, PAI-1)",
}


def read_sheet(sheet):
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = ws.iter_rows(values_only=True)
    header = next(rows)
    data = list(rows)
    df = pd.DataFrame(data, columns=header)
    df = df.set_index("gene")
    return df.apply(pd.to_numeric, errors="coerce")


def pdl_of(col):
    m = re.match(r"PDL(\d+)", col)
    return int(m.group(1)) if m else None


def main():
    log = []
    rs = read_sheet("RS_TC_raw_TPM")
    log.append(f"RS_TC_raw_TPM: {rs.shape[0]} genes x {rs.shape[1]} samples")

    # ---- Task 3: RNA QC ----
    colsum = rs.sum(axis=0)
    detected = (rs > 0).sum(axis=0)
    qc = pd.DataFrame({
        "sample": rs.columns, "pdl": [pdl_of(c) for c in rs.columns],
        "tpm_column_sum": colsum.values,
        "detected_genes_TPM>0": detected.values,
        "detected_genes_pct": (detected.values / rs.shape[0] * 100).round(2),
    })
    qc.to_csv(OUT_T / "gse175533_rna_qc_RS.csv", index=False)
    log.append("\n=== Task 3: RNA QC (RS arm) ===")
    log.append(qc.to_string(index=False))
    log.append(f"\nTPM column sums range {colsum.min():.0f}-{colsum.max():.0f} "
                f"(expected ~1e6 by TPM definition -- {'OK, all close to 1e6' if (colsum.between(9e5,1.1e6)).all() else 'DEVIATION -- investigate'}).")
    log.append("No raw read counts are deposited for this series' bulk RNA -- true sequencing "
                "library size is not computable from GEO supplementary files; "
                "'Average read depth across samples was 50 million paired-end reads' is stated "
                "as a SERIES-WIDE average in Sample_data_processing, not a per-sample value.")

    # ---- Task 4: PCA + Spearman correlation, RS arm ----
    log2tpm = np.log2(rs + 1)
    X = log2tpm.T.values
    X = X - X.mean(axis=0)
    pca = PCA(n_components=5)
    scores = pca.fit_transform(X)
    varexp = pca.explained_variance_ratio_ * 100
    pdls = [pdl_of(c) for c in rs.columns]
    pca_df = pd.DataFrame(scores[:, :3], columns=["PC1", "PC2", "PC3"])
    pca_df.insert(0, "pdl", pdls)
    pca_df.insert(0, "sample", rs.columns)
    pca_df.to_csv(OUT_T / "gse175533_rna_pca_RS.csv", index=False)
    log.append(f"\n=== Task 4: RNA PCA (RS arm, log2(TPM+1)) ===")
    log.append(f"Variance explained: PC1={varexp[0]:.1f}%, PC2={varexp[1]:.1f}%, PC3={varexp[2]:.1f}%")
    log.append(pca_df.to_string(index=False))
    corr_pc1_pdl = np.corrcoef(pca_df["PC1"], pca_df["pdl"])[0, 1]
    log.append(f"\nPearson correlation of PC1 with PDL: r={corr_pc1_pdl:.3f}")

    corr, _ = spearmanr(log2tpm.values)
    corr_df = pd.DataFrame(corr, index=rs.columns, columns=rs.columns)
    corr_df.to_csv(OUT_T / "gse175533_rna_spearman_RS.csv")

    # replicate tightness: mean within-PDL correlation vs mean across all pairs
    within = []
    for pdl in sorted(set(pdls)):
        cols = [c for c, p in zip(rs.columns, pdls) if p == pdl]
        sub = corr_df.loc[cols, cols].values
        n = len(cols)
        if n > 1:
            within.append(sub[np.triu_indices(n, k=1)].mean())
    all_pairs = corr_df.values[np.triu_indices(len(rs.columns), k=1)]
    log.append(f"\nMean within-PDL replicate Spearman rho: {np.mean(within):.4f}")
    log.append(f"Mean across-all-samples Spearman rho: {all_pairs.mean():.4f}")

    # ---- plots ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    sc = axes[0].scatter(pca_df["PC1"], pca_df["PC2"], c=pca_df["pdl"], cmap="viridis", s=70)
    plt.colorbar(sc, ax=axes[0], label="PDL")
    axes[0].set_xlabel(f"PC1 ({varexp[0]:.1f}%)"); axes[0].set_ylabel(f"PC2 ({varexp[1]:.1f}%)")
    axes[0].set_title("GSE175533 RS RNA-seq PCA (log2 TPM)")
    im = axes[1].imshow(corr_df.values, cmap="viridis")
    axes[1].set_title("Spearman correlation (RS RNA, log2 TPM)")
    plt.colorbar(im, ax=axes[1], fraction=0.046)
    plt.tight_layout()
    fig.savefig(OUT_F / "gse175533_rna_pca_correlation.png", dpi=150)
    plt.close(fig)

    # ---- marker panel across PDL course ----
    log.append("\n=== Senescence marker panel (RS arm, mean TPM per PDL) ===")
    marker_rows = []
    pdl_order = sorted(set(pdls))
    mean_by_pdl = pd.DataFrame({pdl: rs.loc[:, [c for c, p in zip(rs.columns, pdls) if p == pdl]].mean(axis=1)
                                 for pdl in pdl_order})
    for gene, expectation in MARKERS.items():
        if gene not in mean_by_pdl.index:
            log.append(f"{gene}: NOT FOUND in TPM table")
            marker_rows.append({"gene": gene, "expected": expectation, "status": "NOT FOUND"})
            continue
        series = mean_by_pdl.loc[gene]
        first, last = series.iloc[0], series.iloc[-1]
        fc = (last + 0.01) / (first + 0.01)
        log2fc = np.log2(fc)
        direction = "UP" if log2fc > 0 else "DOWN"
        expect_dir = "UP" if expectation.startswith("up") else "DOWN"
        match = direction == expect_dir
        log.append(f"{gene}: PDL{pdl_order[0]}={first:.2f} TPM -> PDL{pdl_order[-1]}={last:.2f} TPM, "
                    f"log2FC(last/first)={log2fc:+.2f}, direction={direction}, expected={expect_dir} "
                    f"[{'OK' if match else '*** MISMATCH ***'}]")
        marker_rows.append({
            "gene": gene, "expected": expectation,
            f"TPM_PDL{pdl_order[0]}": round(first, 2), f"TPM_PDL{pdl_order[-1]}": round(last, 2),
            "log2FC_last_vs_first": round(log2fc, 3), "observed_direction": direction,
            "matches_expected": match,
            "full_trajectory_TPM": [round(mean_by_pdl.loc[gene, p], 2) for p in pdl_order],
        })
    marker_df = pd.DataFrame(marker_rows)
    marker_df.to_csv(OUT_T / "gse175533_rna_markers_RS.csv", index=False)
    log.append(f"\nFull PDL order used: {pdl_order}")

    # ---- DESeq2 time-as-numeric result for the same markers (RS_deseq2) ----
    deseq = read_sheet("RS_deseq2")
    log.append("\n=== RS_deseq2 (Wald test, time-as-numeric-covariate) for marker genes ===")
    for gene in MARKERS:
        if gene in deseq.index:
            r = deseq.loc[gene]
            log.append(f"{gene}: log2FC(per unit time)={r['log2FoldChange']:+.3f}, "
                        f"padj={r['padj']:.3g}")
        else:
            log.append(f"{gene}: not in RS_deseq2 table")

    report = "\n".join(str(x) for x in log)
    (OUT_T / "gse175533_rna_qc_signal_log.txt").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
