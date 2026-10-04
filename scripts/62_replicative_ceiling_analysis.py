#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 62_replicative_ceiling_analysis.py
#
# The reproducibility matrix: R1-R3, X1-X4, A1-A2, I1-I3 at four thresholds.
#
# CONSUMES: output/replicative_ceiling_test/anchors_with_hg18.csv, signal_matrix_replicative.csv
# PRODUCES: output/replicative_ceiling_test/replicative_correlations.csv, locus_sanity_checks.csv, batch_pca_stats.csv, step2_stats.json, STEP2_LOG.txt
# ---------------------------------------------------------------------------
"""
WITHIN-REPLICATIVE CEILING TEST -- STEP 2: responses, batch check, correlations, loci.

Normalization identical to the prior enhancer test so numbers are comparable: each sample
divided by its own mean over covered anchors (mean 1), pseudocount 0.10, log2 ratio.
Thresholds identical: |log2FC| >= 0 / 0.25 / 0.5 / 1.0 in at least one axis, 0.5 primary.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
OUT = ROOT / "output" / "replicative_ceiling_test"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PSEUDO = 0.10
THRESHOLDS = [0.0, 0.25, 0.5, 1.0]
PRIMARY = 0.5
LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    adf = pd.read_csv(OUT / "anchors_with_hg18.csv")
    sm = pd.read_csv(OUT / "signal_matrix_replicative.csv")
    R("=" * 84)
    R("WITHIN-REPLICATIVE CEILING TEST -- STEP 2")
    R("=" * 84)
    R(f"anchors {len(adf):,}  |  new signal columns {sm.shape[1]}")

    norm = sm.copy()
    for c in norm.columns:
        mu = np.nanmean(norm[c].to_numpy())
        norm[c] = norm[c] / mu if mu and np.isfinite(mu) else np.nan

    def lfc(on, off):
        a = norm[on].mean(axis=1).to_numpy()
        b = norm[off].mean(axis=1).to_numpy()
        return np.log2((a + PSEUDO) / (b + PSEUDO))

    resp = pd.DataFrame(index=adf.index)
    # NEW replicative datasets
    resp["IMR90_rep_k27ac"] = lfc(["IMR90rep_Sen_R1", "IMR90rep_Sen_R2"],
                                  ["IMR90rep_Young_R1", "IMR90rep_Young_R2"])
    resp["BJ_rep_k27ac"] = lfc(["BJrep_Sen"], ["BJrep_Young"])
    resp["Sen2019_rep_k27ac"] = lfc(["Sen2019_Sen"], ["Sen2019_Pro"])
    # same-lab OIS and the IR arm
    resp["IMR90_OIS_k27ac_samelab"] = lfc(["IMR90ois_RAS"], ["IMR90ois_GFP"])
    resp["Sen2019_IR_k27ac"] = lfc(["Sen2019_IR_R1", "Sen2019_IR_R2"], ["Sen2019_Pro"])
    # carried over from the prior test (identical anchors)
    resp["GM21_OIS_k27ac"] = adf["gm21_k27ac"].to_numpy()
    resp["IMR90_OIS_k27ac_prior"] = adf["imr90_k27ac"].to_numpy()
    resp["GM21_OIS_atac"] = adf["gm21_atac"].to_numpy()
    resp["WI38_rep_atac"] = adf["wi38_atac"].to_numpy()

    N_PER_ARM = {
        "IMR90_rep_k27ac": "n=2 vs 2", "BJ_rep_k27ac": "n=1 vs 1",
        "Sen2019_rep_k27ac": "n=1 vs 1", "IMR90_OIS_k27ac_samelab": "n=1 vs 1",
        "Sen2019_IR_k27ac": "n=2 vs 1", "GM21_OIS_k27ac": "n=2 vs 2",
        "IMR90_OIS_k27ac_prior": "n=2 vs 2", "GM21_OIS_atac": "n=4 vs 4",
        "WI38_rep_atac": "deposited log2FC",
    }
    R("\n[1] RESPONSES (log2 senescent/proliferating)\n" + "-" * 84)
    for c in resp.columns:
        v = resp[c].dropna()
        R(f"  {c:26s} {N_PER_ARM[c]:18s} n={len(v):7,}  mean {v.mean():+.4f}  "
          f"sd {v.std():.4f}")

    # ---------------- batch check ----------------
    R("\n[2] BATCH CHECK -- PCA on the harmonized region x sample matrix\n" + "-" * 84)
    cols = [c for c in sm.columns]
    M = np.log2(norm[cols].to_numpy() + PSEUDO)
    keep = np.isfinite(M).all(axis=1)
    M = M[keep]
    R(f"  matrix {M.shape[0]:,} anchors x {M.shape[1]} samples (complete cases)")
    DS = {"IMR90rep": "GSE146585", "BJrep": "GSE146585", "IMR90ois": "GSE146585",
          "Sen2019": "GSE106146"}
    meta = pd.DataFrame([{
        "sample": c,
        "study": next(v for k, v in DS.items() if c.startswith(k)),
        "status": ("SEN" if ("_Sen" in c or "_RAS" in c or "_IR_" in c) else "PRO"),
    } for c in cols])

    def r2(y, lab):
        y = np.asarray(y, float); gt = y.mean()
        sst = ((y - gt) ** 2).sum()
        if sst == 0:
            return 0.0
        return float(sum((lab == l).sum() * (y[lab == l].mean() - gt) ** 2
                         for l in pd.unique(lab)) / sst)

    def pca(X, tag):
        Z = X.T - X.T.mean(0, keepdims=True)
        U, S, Vt = np.linalg.svd(Z, full_matrices=False)
        fr = S ** 2 / (S ** 2).sum(); sc = U * S
        rows = []
        for k in range(min(3, sc.shape[1])):
            a, b = r2(sc[:, k], meta["study"].to_numpy()), r2(sc[:, k], meta["status"].to_numpy())
            rows.append(dict(space=tag, pc=f"PC{k+1}", var=float(fr[k]), r2_study=a, r2_status=b))
            R(f"  {tag:20s} PC{k+1}  var {100*fr[k]:5.1f}%  R2(study) {a:.3f}  "
              f"R2(status) {b:.3f}  -> {'STUDY' if a>b else 'STATUS'} dominates")
        return sc, fr, rows

    sc_r, fr_r, r1 = pca(M, "raw pooled")
    Mc = M.copy()
    for st in meta["study"].unique():
        m = (meta["study"] == st).to_numpy()
        Mc[:, m] = M[:, m] - M[:, m].mean(axis=1, keepdims=True)
    R("")
    sc_c, fr_c, r2rows = pca(Mc, "per-study centered")
    pd.DataFrame(r1 + r2rows).to_csv(OUT / "batch_pca_stats.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (sc, fr, t) in zip(axes, [(sc_r, fr_r, "raw pooled"), (sc_c, fr_c, "per-study centered")]):
        for st, mk in zip(meta["study"].unique(), ["o", "s"]):
            for sen, col in [("PRO", "#1f77b4"), ("SEN", "#d62728")]:
                m = ((meta["study"] == st) & (meta["status"] == sen)).to_numpy()
                if m.any():
                    ax.scatter(sc[m, 0], sc[m, 1], marker=mk, c=col, s=110, edgecolors="k",
                               linewidths=.5, label=f"{st} {sen}" if t == "raw pooled" else None)
        ax.set_title(f"{t}\nPC1 {100*fr[0]:.1f}%  PC2 {100*fr[1]:.1f}%")
        ax.axhline(0, lw=.5, c="#ccc"); ax.axvline(0, lw=.5, c="#ccc")
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    axes[0].legend(fontsize=7)
    fig.suptitle("Replicative H3K27ac: batch structure (shape = study, colour = status)")
    fig.tight_layout(); fig.savefig(FIG / "replicative_batch_pca.png", dpi=150)
    R(f"  figure -> {FIG/'replicative_batch_pca.png'}")

    # ---------------- correlations ----------------
    R("\n[3] CORRELATIONS\n" + "-" * 84)
    COMPS = [
        ("R1", "IMR90 replicative K27ac (n=2v2)", "IMR90_rep_k27ac",
         "GSE106146 replicative K27ac (n=1v1)", "Sen2019_rep_k27ac",
         "*** WITHIN-REPLICATIVE CEILING -- cross-lab, cross-build. Analog of OIS +0.653"),
        ("R2", "IMR90 replicative K27ac (n=2v2)", "IMR90_rep_k27ac",
         "BJ replicative K27ac (n=1v1)", "BJ_rep_k27ac",
         "within-replicative, SAME lab/pipeline/build, different cell line (upper bound)"),
        ("R3", "BJ replicative K27ac (n=1v1)", "BJ_rep_k27ac",
         "GSE106146 replicative K27ac (n=1v1)", "Sen2019_rep_k27ac",
         "within-replicative, cross-lab"),
        ("X1", "IMR90 replicative K27ac (n=2v2)", "IMR90_rep_k27ac",
         "IMR90 OIS K27ac SAME LAB (n=1v1)", "IMR90_OIS_k27ac_samelab",
         "CROSS-MECHANISM, same lab/cell line/pipeline/build -- cleanest possible"),
        ("X2", "IMR90 replicative K27ac (n=2v2)", "IMR90_rep_k27ac",
         "GM21 OIS K27ac (n=2v2)", "GM21_OIS_k27ac", "CROSS-MECHANISM, cross-lab"),
        ("X3", "IMR90 replicative K27ac (n=2v2)", "IMR90_rep_k27ac",
         "IMR90 OIS K27ac Tasdemir (n=2v2)", "IMR90_OIS_k27ac_prior",
         "CROSS-MECHANISM, same cell line, different lab"),
        ("X4", "GSE106146 replicative K27ac (n=1v1)", "Sen2019_rep_k27ac",
         "GM21 OIS K27ac (n=2v2)", "GM21_OIS_k27ac", "CROSS-MECHANISM, cross-lab"),
        ("A1", "WI-38 replicative ATAC", "WI38_rep_atac",
         "IMR90 replicative K27ac (n=2v2)", "IMR90_rep_k27ac",
         "within-replicative, CROSS-ASSAY (analog of OIS +0.607)"),
        ("A2", "WI-38 replicative ATAC", "WI38_rep_atac",
         "GSE106146 replicative K27ac (n=1v1)", "Sen2019_rep_k27ac",
         "within-replicative, cross-assay"),
        ("I1", "GSE106146 IR senescence K27ac (n=2v1)", "Sen2019_IR_k27ac",
         "GSE106146 replicative K27ac (n=1v1)", "Sen2019_rep_k27ac",
         "IR vs replicative, SAME study/pipeline (shares the proliferating control)"),
        ("I2", "GSE106146 IR senescence K27ac (n=2v1)", "Sen2019_IR_k27ac",
         "GM21 OIS K27ac (n=2v2)", "GM21_OIS_k27ac", "IR vs OIS, cross-lab"),
        ("I3", "GSE106146 IR senescence K27ac (n=2v1)", "Sen2019_IR_k27ac",
         "IMR90 replicative K27ac (n=2v2)", "IMR90_rep_k27ac", "IR vs replicative, cross-lab"),
    ]
    rows = []
    for tag, n1, c1, n2, c2, note in COMPS:
        R(f"\n  ({tag}) {n1}  vs  {n2}")
        R(f"       {note}")
        for t in THRESHOLDS:
            m = resp[[c1, c2]].dropna()
            if t > 0:
                m = m[(m[c1].abs() >= t) | (m[c2].abs() >= t)]
            if len(m) < 20:
                R(f"       |lfc|>={t:<4} n={len(m):<7,} TOO FEW REGIONS")
                continue
            sp = stats.spearmanr(m[c1], m[c2]); pe = stats.pearsonr(m[c1], m[c2])
            star = "  <-- PRIMARY" if t == PRIMARY else ""
            R(f"       |lfc|>={t:<4} n={len(m):<8,} Spearman {sp.statistic:+.4f}   "
              f"Pearson {pe.statistic:+.4f}{star}")
            rows.append(dict(comparison=tag, label=f"{n1} vs {n2}", note=note,
                             threshold=t, n=int(len(m)),
                             spearman=float(sp.statistic), spearman_p=float(sp.pvalue),
                             pearson=float(pe.statistic)))
    cdf = pd.DataFrame(rows)
    cdf.to_csv(OUT / "replicative_correlations.csv", index=False)

    # scatter grid at primary threshold
    key = ["R1", "R2", "X1", "X2", "A1", "I1"]
    sel = [c for c in COMPS if c[0] in key]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for ax, (tag, n1, c1, n2, c2, note) in zip(axes.ravel(), sel):
        m = resp[[c1, c2]].dropna()
        m = m[(m[c1].abs() >= PRIMARY) | (m[c2].abs() >= PRIMARY)]
        if len(m) >= 20:
            ax.scatter(m[c1], m[c2], s=3, alpha=.12, edgecolors="none", c="#222")
            sp = stats.spearmanr(m[c1], m[c2]).statistic
            ax.set_title(f"({tag}) rho={sp:+.3f}  n={len(m):,}", fontsize=10)
        else:
            ax.set_title(f"({tag}) too few regions", fontsize=10)
        ax.set_xlabel(n1, fontsize=7); ax.set_ylabel(n2, fontsize=7)
        ax.axhline(0, lw=.5, c="#c33"); ax.axvline(0, lw=.5, c="#c33")
    fig.suptitle(f"Replicative-ceiling test: per-region H3K27ac response, "
                 f"|log2FC| >= {PRIMARY} in at least one axis")
    fig.tight_layout(); fig.savefig(FIG / "replicative_ceiling_scatter.png", dpi=150)
    R(f"\n  figure -> {FIG/'replicative_ceiling_scatter.png'}")

    # ---------------- loci ----------------
    R("\n[4] SANITY -- CDKN2A / CDKN1A promoter vs distal (H3K27ac)\n" + "-" * 84)
    LOCI = {"CDKN2A": ("chr9", 21_967_751, 21_995_300),
            "CDKN1A": ("chr6", 36_644_237, 36_655_116)}
    FLANK = 200_000
    srows = []
    COLS = [("IMR90_rep_k27ac", "IMR90 rep"), ("BJ_rep_k27ac", "BJ rep"),
            ("Sen2019_rep_k27ac", "GSE106146 rep"), ("Sen2019_IR_k27ac", "GSE106146 IR"),
            ("IMR90_OIS_k27ac_samelab", "IMR90 OIS(same lab)"),
            ("GM21_OIS_k27ac", "GM21 OIS"), ("IMR90_OIS_k27ac_prior", "IMR90 OIS(prior)")]
    for g, (c, s, e) in LOCI.items():
        sel_m = ((adf["chrom_hg19"] == c) & (adf["mid_hg19"] >= s - FLANK) &
                 (adf["mid_hg19"] <= e + FLANK)).to_numpy()
        d = adf[sel_m]; rr = resp[sel_m]
        dist = np.where(d["mid_hg19"] < s, s - d["mid_hg19"],
                        np.where(d["mid_hg19"] > e, d["mid_hg19"] - e, 0))
        prox = dist <= 2000
        R(f"\n  {g} ({c}:{s:,}-{e:,} hg19 +/-{FLANK//1000} kb): {int(sel_m.sum())} anchors "
          f"(promoter <=2 kb: {int(prox.sum())}, distal: {int((~prox).sum())})")
        for lab, m in [("promoter", prox), ("distal", ~prox)]:
            if m.sum() == 0:
                R(f"    {lab:9s} none"); continue
            row = dict(locus=g, region_class=lab, n=int(m.sum()))
            parts = []
            for col, nm in COLS:
                v = rr.loc[m, col].dropna()
                row[f"{col}_mean"] = float(v.mean()) if len(v) else None
                row[f"{col}_n"] = int(len(v))
                parts.append(f"{nm} {('%+.3f' % v.mean()) if len(v) else ' n/a'}")
            srows.append(row)
            R(f"    {lab:9s} " + " | ".join(parts))
    pd.DataFrame(srows).to_csv(OUT / "locus_sanity_checks.csv", index=False)

    full = pd.concat([adf.reset_index(drop=True), resp.reset_index(drop=True)], axis=1)
    full.to_csv(OUT / "region_responses_replicative.csv", index=False)
    json.dump({"pseudocount": PSEUDO, "primary_threshold": PRIMARY,
               "n_per_arm": N_PER_ARM}, open(OUT / "step2_stats.json", "w"), indent=2)
    (OUT / "STEP2_LOG.txt").write_text("\n".join(LOG))
    R(f"\n  wrote region_responses_replicative.csv {full.shape}")
    R("STEP2_DONE")


if __name__ == "__main__":
    main()
