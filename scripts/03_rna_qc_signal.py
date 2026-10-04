#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 03_rna_qc_signal.py
#
# RNA QC: library size, detected genes, PCA, marker check.
#
# CONSUMES: data/<ACC>/ RNA count tables
# PRODUCES: output/tables/rna_qc_*.csv, rna_pca_scores.csv, rna_senescence_markers.csv, rna_qc_signal_log.txt
# ---------------------------------------------------------------------------
"""
Task 3 (processed-data QC) and Task 4 (signal check) for GSE220545 RNA-seq.

Input: GSE220545_72h_vs_30h_DESEQ.txt -- the ONLY supplementary file GEO has
for this series. It is a DESeq2 results table for the 72h-vs-30h contrast,
but its last 6 columns (A30h, A72h, B30h, B72h, C30h, C72h) carry DESeq2
size-factor-normalized per-sample counts, so a per-sample expression matrix
can be reconstructed from it even though no standalone counts/TPM matrix
was deposited.

ASSUMPTIONS (stated explicitly per the task instructions):
  1. Column letters A/B/C map 1:1 to "Replicate A/B/C" in Sample_title -- this
     is the only reasonable reading and is internally consistent (3 replicates
     x 2 timepoints = 6 columns = 6 GSMs).
  2. These are DESeq2 median-of-ratios normalized counts, NOT raw counts.
     "Library size" computed from them is therefore each sample's total
     normalized signal, not its raw sequencing depth -- flagged in the output.
  3. Column headers "ko.sd" / "wt.sd" appear to be leftover generic DESeq2
     template names (knockout/wild-type) reused for 72h/30h respectively
     (they sit immediately after the "72h mean (log2)" / "30h mean (log2)"
     columns) -- not used in this analysis, noted here for transparency only.
  4. No ATAC-seq matrix exists for GSE254358 (accession does not exist), so
     Task 3/4's ATAC arm cannot be executed at all -- not attempted below.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
DESEQ = ROOT / "data" / "GSE220545" / "GSE220545_72h_vs_30h_DESEQ.txt"
OUT_TABLES = ROOT / "output" / "tables"
OUT_FIG = ROOT / "output" / "figures"
OUT_TABLES.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SAMPLE_COLS = ["A30h", "B30h", "C30h", "A72h", "B72h", "C72h"]
GSM_MAP = {
    "A30h": "GSM6806685", "B30h": "GSM6806686", "C30h": "GSM6806687",
    "A72h": "GSM6806688", "B72h": "GSM6806689", "C72h": "GSM6806690",
}
STATUS_MAP = {c: ("uncommitted_30h" if "30h" in c else "committed_72h") for c in SAMPLE_COLS}

MARKERS = {
    "CDKN2A": "up (p16, senescence effector)",
    "CDKN1A": "up (p21, senescence effector)",
    "LMNB1": "down (lamin B1 loss is a core senescence marker)",
    "MKI67": "down (proliferation marker)",
    "IL6": "up (core SASP factor)",
    "CXCL8": "up (SASP factor, aka IL8)",
}


def load_matrix():
    df = pd.read_csv(DESEQ, sep="\t")
    df = df.rename(columns={df.columns[0]: "row_id"})
    # de-dup gene symbols (keep highest baseMean row per symbol) for marker lookups
    df["Gene"] = df["Gene"].astype(str)
    return df


def main():
    log = []

    df = load_matrix()
    log.append(f"Loaded {DESEQ.name}: {df.shape[0]} rows x {df.shape[1]} columns")
    missing = [c for c in SAMPLE_COLS if c not in df.columns]
    if missing:
        sys.exit(f"FATAL: expected sample columns missing: {missing}")

    counts = df.set_index("EnsemblID")[SAMPLE_COLS].apply(pd.to_numeric, errors="coerce")
    counts = counts.dropna(how="all")
    log.append(f"Count matrix (genes x samples): {counts.shape}")

    # ---- Task 3: library size + detected-gene QC ----
    lib_size = counts.sum(axis=0)
    detected = (counts > 0).sum(axis=0)
    qc = pd.DataFrame({
        "gsm": [GSM_MAP[c] for c in SAMPLE_COLS],
        "column": SAMPLE_COLS,
        "status": [STATUS_MAP[c] for c in SAMPLE_COLS],
        "normalized_library_size": lib_size[SAMPLE_COLS].values,
        "detected_genes_count>0": detected[SAMPLE_COLS].values,
        "detected_genes_pct_of_total": (detected[SAMPLE_COLS].values / counts.shape[0] * 100).round(2),
    })
    qc.to_csv(OUT_TABLES / "rna_qc_library_detected_genes.csv", index=False)
    log.append("\n=== RNA QC: normalized library size & detected genes ===")
    log.append(qc.to_string(index=False))
    log.append("NOTE: 'normalized_library_size' is sum of DESeq2 size-factor-normalized "
                "counts, NOT raw read depth -- GEO does not provide raw counts or FASTQ/"
                "SRA-derived depth for this series (per task constraints, SRA was not queried).")

    # ---- Task 4: CPM-normalize (on top of the already-normalized counts, for a
    # comparable per-million scale), log2-transform, PCA, Spearman correlation ----
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    log2cpm = np.log2(cpm + 1)

    # keep expressed genes only (avoid PCA being dominated by all-zero rows)
    expressed = log2cpm[(counts > 0).sum(axis=1) >= 3]
    log.append(f"\nGenes retained for PCA/correlation (detected in >=3/6 samples): {expressed.shape[0]}")

    X = expressed.T.values  # samples x genes
    X = X - X.mean(axis=0)  # center genes
    pca = PCA(n_components=min(5, X.shape[0] - 1))
    scores = pca.fit_transform(X)
    var_explained = pca.explained_variance_ratio_ * 100

    pca_df = pd.DataFrame(scores[:, :3], columns=["PC1", "PC2", "PC3"])
    pca_df.insert(0, "status", [STATUS_MAP[c] for c in SAMPLE_COLS])
    pca_df.insert(0, "gsm", [GSM_MAP[c] for c in SAMPLE_COLS])
    pca_df.insert(0, "sample_column", SAMPLE_COLS)
    pca_df.to_csv(OUT_TABLES / "rna_pca_scores.csv", index=False)

    log.append("\n=== RNA PCA ===")
    log.append(f"Variance explained: PC1={var_explained[0]:.1f}%, PC2={var_explained[1]:.1f}%, "
                f"PC3={var_explained[2]:.1f}%")
    log.append(pca_df.to_string(index=False))

    pc1_by_status = pca_df.groupby("status")["PC1"].agg(["mean", "std"])
    sep_pc1 = abs(pc1_by_status.loc["uncommitted_30h", "mean"] - pc1_by_status.loc["committed_72h", "mean"])
    within_std = pc1_by_status["std"].mean()
    log.append(f"\nPC1 separation: |mean diff| = {sep_pc1:.2f}, mean within-group std = {within_std:.2f} "
                f"-> separation/within-group-spread ratio = {sep_pc1/within_std:.2f}")
    separates_pc1 = sep_pc1 > 2 * within_std
    log.append(f"Conclusion: samples {'DO' if separates_pc1 else 'DO NOT clearly'} separate by status on PC1.")

    # Spearman correlation heatmap data
    corr, pval = spearmanr(expressed.values)
    corr_df = pd.DataFrame(corr, index=SAMPLE_COLS, columns=SAMPLE_COLS)
    corr_df.to_csv(OUT_TABLES / "rna_spearman_correlation.csv")
    log.append("\n=== Spearman correlation matrix (log2 CPM, expressed genes) ===")
    log.append(corr_df.round(3).to_string())

    within_30h = corr_df.loc[["A30h", "B30h", "C30h"], ["A30h", "B30h", "C30h"]].values
    within_72h = corr_df.loc[["A72h", "B72h", "C72h"], ["A72h", "B72h", "C72h"]].values
    between = corr_df.loc[["A30h", "B30h", "C30h"], ["A72h", "B72h", "C72h"]].values
    within_30h_mean = within_30h[np.triu_indices(3, k=1)].mean()
    within_72h_mean = within_72h[np.triu_indices(3, k=1)].mean()
    between_mean = between.mean()
    log.append(f"\nMean within-30h replicate rho = {within_30h_mean:.4f}")
    log.append(f"Mean within-72h replicate rho = {within_72h_mean:.4f}")
    log.append(f"Mean between-group rho        = {between_mean:.4f}")

    # ---- plots ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"uncommitted_30h": "#4C72B0", "committed_72h": "#C44E52"}
    for status, sub in pca_df.groupby("status"):
        axes[0].scatter(sub["PC1"], sub["PC2"], label=status, s=80, color=colors[status])
        for _, r in sub.iterrows():
            axes[0].annotate(r["sample_column"], (r["PC1"], r["PC2"]), fontsize=8,
                              xytext=(4, 4), textcoords="offset points")
    axes[0].set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    axes[0].set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
    axes[0].set_title("GSE220545 RNA-seq PCA (log2 CPM)")
    axes[0].legend()

    im = axes[1].imshow(corr_df.values, vmin=corr_df.values.min(), vmax=1.0, cmap="viridis")
    axes[1].set_xticks(range(6)); axes[1].set_xticklabels(SAMPLE_COLS, rotation=45, ha="right")
    axes[1].set_yticks(range(6)); axes[1].set_yticklabels(SAMPLE_COLS)
    axes[1].set_title("Spearman correlation (log2 CPM)")
    plt.colorbar(im, ax=axes[1], fraction=0.046)
    plt.tight_layout()
    fig.savefig(OUT_FIG / "rna_pca_and_correlation.png", dpi=150)
    plt.close(fig)
    log.append(f"\nSaved figure: {OUT_FIG / 'rna_pca_and_correlation.png'}")

    # ---- Task 4: senescence marker check ----
    log.append("\n=== Senescence marker check (72h/committed vs 30h/uncommitted) ===")
    marker_rows = []
    for gene, expectation in MARKERS.items():
        sub = df[df["Gene"] == gene]
        if sub.empty:
            marker_rows.append({"gene": gene, "expected": expectation, "status": "NOT FOUND IN TABLE"})
            log.append(f"{gene}: NOT FOUND in supplementary table")
            continue
        # if multiple Ensembl IDs map to the same symbol, report the one with highest baseMean
        sub = sub.sort_values("baseMean", ascending=False)
        r = sub.iloc[0]
        lfc = r["log2FoldChange"]
        padj = r["padj"]
        direction = "UP" if lfc > 0 else ("DOWN" if lfc < 0 else "NO CHANGE")
        expected_dir = "UP" if "up" in expectation.split()[0] else "DOWN"
        matches = (direction == expected_dir)
        fold = 2 ** lfc
        marker_rows.append({
            "gene": gene, "ensembl_id": r["EnsemblID"], "expected": expectation,
            "log2FC_72h_vs_30h": round(lfc, 3), "fold_change": round(fold, 3),
            "padj": r["padj"], "observed_direction": direction,
            "matches_expected_direction": matches,
        })
        flag = "OK" if matches else "*** MISMATCH ***"
        log.append(f"{gene}: log2FC={lfc:+.3f} (fold={fold:.2f}x), padj={padj:.3g}, "
                    f"direction={direction}, expected={expected_dir} [{flag}]")

    marker_df = pd.DataFrame(marker_rows)
    marker_df.to_csv(OUT_TABLES / "rna_senescence_markers.csv", index=False)

    n_mismatch = (~marker_df.get("matches_expected_direction", pd.Series(dtype=bool))).sum() \
        if "matches_expected_direction" in marker_df else None
    if "matches_expected_direction" in marker_df:
        total = marker_df["matches_expected_direction"].notna().sum()
        n_ok = marker_df["matches_expected_direction"].sum()
        log.append(f"\nMarker directions matching expectation: {n_ok}/{total}")

    report_text = "\n".join(str(x) for x in log)
    (OUT_TABLES / "rna_qc_signal_log.txt").write_text(report_text)
    print(report_text)


if __name__ == "__main__":
    main()
