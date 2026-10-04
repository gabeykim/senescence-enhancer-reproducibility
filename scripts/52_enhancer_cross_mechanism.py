#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 52_enhancer_cross_mechanism.py
#
# Cross-mechanism correlations; ATAC-vs-H3K27ac concordance; batch PCA.
#
# CONSUMES: output/enhancer_test/anchors.csv, signal_matrix.csv
# PRODUCES: output/enhancer_test/cross_mechanism_correlations.csv, step2_stats.json, locus_sanity_checks.csv, batch_pca_stats.csv, STEP2_LOG.txt
# ---------------------------------------------------------------------------
"""
Cross-mechanism enhancer kill test -- STEP 2: responses, batch check, correlations.

NORMALIZATION, stated: each bigWig sample's per-anchor mean coverage is divided by that
sample's mean across all covered anchors, so every sample has mean 1 over the anchor set.
This is a within-sample depth/scale normalization; it cannot fix genuine pipeline
differences between studies, which is why the PCA batch check in step [2] is run before
any correlation is believed. A pseudocount of PSEUDO is added before the log ratio so that
anchors with near-zero coverage do not produce unbounded fold changes.

Responses (all log2, senescent vs proliferating):
  WI-38 replicative ATAC : deposited log2FC, PDL50_v_hTERT - mean(PDL25/33/37_v_hTERT)
  GM21 OIS ATAC          : RAS D18+D23 vs pBABE D8+D32
  GM21 OIS H3K27ac       : RAS D18 vs EV
  IMR90 OIS H3K27ac      : senescent vs proliferating (2 reps each)
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
OUT = ROOT / "output" / "enhancer_test"
FIG = ROOT / "output" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

PSEUDO = 0.10                       # after per-sample mean-1 normalization
THRESHOLDS = [0.0, 0.25, 0.5, 1.0]  # |log2FC| in >=1 study; 0.5 is primary
PRIMARY_T = 0.5

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    adf = pd.read_csv(OUT / "anchors.csv")
    sm = pd.read_csv(OUT / "signal_matrix.csv")
    R("=" * 84)
    R("CROSS-MECHANISM ENHANCER TEST -- STEP 2")
    R("=" * 84)
    R(f"anchors {len(adf):,} | signal columns {sm.shape[1]}")
    for c in sm.columns:
        R(f"   {c}")

    # ---------------- normalize ----------------
    norm = sm.copy()
    for c in norm.columns:
        mu = np.nanmean(norm[c].to_numpy())
        norm[c] = norm[c] / mu if mu and np.isfinite(mu) else np.nan
    R(f"\n[1] normalization: each sample divided by its mean over anchors "
      f"(all samples now mean 1); pseudocount {PSEUDO}")

    def grp(prefix, key=None):
        cols = [c for c in norm.columns if c.startswith(prefix)
                and (key is None or key in c)]
        return cols

    atac_ev = grp("GM21_ATAC::", "pBABE")
    atac_sen = grp("GM21_ATAC::", "RAS")
    k27_ev = grp("GM21_K27ac::", "EV")
    k27_sen = grp("GM21_K27ac::", "RAS")
    imr_pro = grp("IMR90_K27ac::PRO")
    imr_sen = grp("IMR90_K27ac::SEN")
    R(f"\n  GM21 ATAC   EV {len(atac_ev)}  SEN {len(atac_sen)}")
    R(f"  GM21 K27ac  EV {len(k27_ev)}  SEN {len(k27_sen)}")
    R(f"  IMR90 K27ac PRO {len(imr_pro)}  SEN {len(imr_sen)}")

    def lfc(on_cols, off_cols):
        on = norm[on_cols].mean(axis=1).to_numpy()
        off = norm[off_cols].mean(axis=1).to_numpy()
        return np.log2((on + PSEUDO) / (off + PSEUDO))

    resp = pd.DataFrame(index=adf.index)
    resp["gm21_atac"] = lfc(atac_sen, atac_ev)
    resp["gm21_k27ac"] = lfc(k27_sen, k27_ev)
    resp["imr90_k27ac"] = lfc(imr_sen, imr_pro)
    resp["wi38_atac"] = adf["wi38_response"].to_numpy()
    # a response is only defined where the underlying anchor is in that assay's consensus
    resp.loc[~adf["in_gm21_atac"].to_numpy(), "gm21_atac"] = np.nan
    resp.loc[~adf["in_gm21_k27ac"].to_numpy(), "gm21_k27ac"] = np.nan

    R("\n  response summary (log2 senescent/proliferating):")
    for c in resp.columns:
        v = resp[c].dropna()
        R(f"    {c:14s} n={len(v):7,}  mean {v.mean():+.4f}  sd {v.std():.4f}  "
          f"range [{v.min():+.3f}, {v.max():+.3f}]")

    # ---------------- [2] BATCH CHECK ----------------
    R("\n[2] BATCH CHECK -- PCA on the harmonized region x sample matrix\n" + "-" * 84)
    mat_cols = atac_ev + atac_sen + k27_ev + k27_sen + imr_pro + imr_sen
    M = np.log2(norm[mat_cols].to_numpy() + PSEUDO)
    keep = np.isfinite(M).all(axis=1)
    M = M[keep]
    R(f"  matrix {M.shape[0]:,} anchors x {M.shape[1]} samples (complete cases)")
    meta = []
    for c in mat_cols:
        st = ("GSE206402" if c.startswith("GM21_ATAC") else
              "GSE205898" if c.startswith("GM21_K27ac") else "GSE74238")
        sen = ("SEN" if ("RAS" in c or "::SEN::" in c) else "PRO")
        meta.append((c, st, sen))
    meta = pd.DataFrame(meta, columns=["sample", "study", "status"])

    def anova_r2(y, lab):
        y = np.asarray(y, float); gt = y.mean()
        sst = ((y - gt) ** 2).sum()
        if sst == 0:
            return 0.0
        ssb = sum((lab == l).sum() * (y[lab == l].mean() - gt) ** 2 for l in pd.unique(lab))
        return float(ssb / sst)

    def pca(X, tag):
        Z = X.T - X.T.mean(0, keepdims=True)
        U, S, Vt = np.linalg.svd(Z, full_matrices=False)
        fr = S ** 2 / (S ** 2).sum()
        sc = U * S
        rows = []
        for k in range(min(3, sc.shape[1])):
            r_st = anova_r2(sc[:, k], meta["study"].to_numpy())
            r_se = anova_r2(sc[:, k], meta["status"].to_numpy())
            rows.append(dict(space=tag, pc=f"PC{k+1}", var=float(fr[k]),
                             r2_study=r_st, r2_status=r_se))
            R(f"  {tag:18s} PC{k+1}  var {100*fr[k]:5.1f}%  R2(study) {r_st:.3f}  "
              f"R2(status) {r_se:.3f}  -> {'STUDY' if r_st>r_se else 'STATUS'} dominates")
        return sc, fr, rows

    sc_raw, fr_raw, rows1 = pca(M, "raw pooled")
    Mc = M.copy()
    for st in meta["study"].unique():
        cols = (meta["study"] == st).to_numpy()
        Mc[:, cols] = M[:, cols] - M[:, cols].mean(axis=1, keepdims=True)
    R("")
    sc_c, fr_c, rows2 = pca(Mc, "per-study centered")
    pd.DataFrame(rows1 + rows2).to_csv(OUT / "batch_pca_stats.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (sc, fr, t) in zip(axes, [(sc_raw, fr_raw, "raw pooled"),
                                      (sc_c, fr_c, "per-study centered")]):
        for st, mk in zip(meta["study"].unique(), ["o", "s", "^"]):
            for sen, col in [("PRO", "#1f77b4"), ("SEN", "#d62728")]:
                m = ((meta["study"] == st) & (meta["status"] == sen)).to_numpy()
                if m.any():
                    ax.scatter(sc[m, 0], sc[m, 1], marker=mk, c=col, s=90,
                               edgecolors="k", linewidths=.5,
                               label=f"{st} {sen}" if t == "raw pooled" else None)
        ax.set_title(f"{t}\nPC1 {100*fr[0]:.1f}%  PC2 {100*fr[1]:.1f}%")
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        ax.axhline(0, lw=.5, c="#ccc"); ax.axvline(0, lw=.5, c="#ccc")
    axes[0].legend(fontsize=7)
    fig.suptitle("Enhancer signal: batch structure (shape = study, colour = senescence status)")
    fig.tight_layout(); fig.savefig(FIG / "enhancer_batch_pca.png", dpi=150)
    R(f"  figure -> {FIG/'enhancer_batch_pca.png'}")

    # ---------------- [3] CORRELATIONS ----------------
    R("\n[3] CROSS-MECHANISM CORRELATIONS\n" + "-" * 84)
    COMPARISONS = [
        ("a", "WI-38 replicative ATAC", "wi38_atac", "GM21 OIS ATAC", "gm21_atac",
         "CROSS-MECHANISM (direct analog of the expression test)"),
        ("b", "GM21 OIS H3K27ac", "gm21_k27ac", "IMR90 OIS H3K27ac", "imr90_k27ac",
         "WITHIN-MECHANISM ceiling (different labs & cell lines)"),
        ("c", "WI-38 replicative ATAC", "wi38_atac", "GM21 OIS H3K27ac", "gm21_k27ac",
         "CROSS-MECHANISM, CROSS-ASSAY"),
        ("d", "GM21 OIS ATAC", "gm21_atac", "GM21 OIS H3K27ac", "gm21_k27ac",
         "WITHIN-STUDY, CROSS-ASSAY"),
    ]
    corr_rows = []
    for tag, n1, c1, n2, c2, note in COMPARISONS:
        R(f"\n  ({tag}) {n1}  vs  {n2}")
        R(f"       {note}")
        for t in THRESHOLDS:
            m = resp[[c1, c2]].dropna()
            if t > 0:
                m = m[(m[c1].abs() >= t) | (m[c2].abs() >= t)]
            if len(m) < 20:
                R(f"       |lfc|>={t:<4} n={len(m):<7,} TOO FEW REGIONS")
                corr_rows.append(dict(comparison=tag, threshold=t, n=len(m),
                                      spearman=None, pearson=None))
                continue
            sp = stats.spearmanr(m[c1], m[c2])
            pe = stats.pearsonr(m[c1], m[c2])
            star = "  <-- PRIMARY" if t == PRIMARY_T else ""
            R(f"       |lfc|>={t:<4} n={len(m):<7,} Spearman {sp.statistic:+.4f} "
              f"(p={sp.pvalue:.2e})   Pearson {pe.statistic:+.4f}{star}")
            corr_rows.append(dict(comparison=tag, label=f"{n1} vs {n2}", threshold=t,
                                  n=int(len(m)), spearman=float(sp.statistic),
                                  spearman_p=float(sp.pvalue), pearson=float(pe.statistic)))
    cdf = pd.DataFrame(corr_rows)
    cdf.to_csv(OUT / "cross_mechanism_correlations.csv", index=False)

    # scatter plots at the primary threshold
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, (tag, n1, c1, n2, c2, note) in zip(axes, COMPARISONS):
        m = resp[[c1, c2]].dropna()
        m = m[(m[c1].abs() >= PRIMARY_T) | (m[c2].abs() >= PRIMARY_T)]
        if len(m) >= 20:
            ax.scatter(m[c1], m[c2], s=3, alpha=.15, edgecolors="none", c="#333")
            sp = stats.spearmanr(m[c1], m[c2]).statistic
            ax.set_title(f"({tag}) rho={sp:+.3f}  n={len(m):,}", fontsize=10)
        else:
            ax.set_title(f"({tag}) too few regions", fontsize=10)
        ax.set_xlabel(n1, fontsize=8); ax.set_ylabel(n2, fontsize=8)
        ax.axhline(0, lw=.5, c="#c33"); ax.axvline(0, lw=.5, c="#c33")
    fig.suptitle(f"Per-region senescence response, |log2FC| >= {PRIMARY_T} in at least one axis")
    fig.tight_layout(); fig.savefig(FIG / "enhancer_cross_mechanism_scatter.png", dpi=150)
    R(f"\n  figure -> {FIG/'enhancer_cross_mechanism_scatter.png'}")

    # ---------------- [4] ACTIVE vs MERELY OPEN ----------------
    R("\n[4] ACTIVE vs MERELY OPEN (GM21, both assays in the same study)\n" + "-" * 84)
    both = resp[["gm21_atac", "gm21_k27ac"]].dropna()
    act = {}
    for t in [0.5, 1.0]:
        da = both[both["gm21_atac"].abs() >= t]
        conc = (da["gm21_k27ac"].abs() >= t)
        same = conc & (np.sign(da["gm21_k27ac"]) == np.sign(da["gm21_atac"]))
        act[str(t)] = dict(n_diff_atac=int(len(da)),
                           n_also_diff_k27ac=int(conc.sum()),
                           frac_also_diff=float(conc.mean()) if len(da) else None,
                           n_same_direction=int(same.sum()),
                           frac_same_direction=float(same.mean()) if len(da) else None)
        R(f"  |log2FC| >= {t}: {int(conc.sum()):,}/{len(da):,} "
          f"({100*conc.mean():.1f}%) of differential-ACCESSIBILITY anchors are also "
          f"differential in H3K27ac; {100*same.mean():.1f}% agree in direction")
    R(f"  anchors with both assays measured: {len(both):,}")

    # ---------------- [5] SANITY: CDKN2A / CDKN1A ----------------
    R("\n[5] SANITY CHECK -- CDKN2A (chr9) and CDKN1A (chr6) loci\n" + "-" * 84)
    R("  Promoter-proximal vs distal, using GSE175533's own peak annotation where the")
    R("  anchor matched a WI-38 peak. Prior work: both PROMOTERS show no significant")
    R("  accessibility change -- the enhancer premise is that DISTAL regions do.")
    LOCI = {"CDKN2A": ("chr9", 21_967_751, 21_995_300),      # hg19
            "CDKN1A": ("chr6", 36_644_237, 36_655_116)}      # hg19
    FLANK = 200_000
    sane = []
    for g, (c, s, e) in LOCI.items():
        sel = ((adf["chrom_hg19"] == c) &
               (adf["mid_hg19"] >= s - FLANK) & (adf["mid_hg19"] <= e + FLANK)).to_numpy()
        d = adf[sel].copy()
        rr = resp[sel].copy()
        dist = np.where(d["mid_hg19"] < s, s - d["mid_hg19"],
                        np.where(d["mid_hg19"] > e, d["mid_hg19"] - e, 0))
        rr["dist_to_gene"] = dist
        prox = dist <= 2000
        R(f"\n  {g}  ({c}:{s:,}-{e:,} hg19 +/-{FLANK//1000} kb)  anchors: {len(d):,} "
          f"(promoter-proximal <=2 kb: {int(prox.sum())}, distal: {int((~prox).sum())})")
        for lab, m in [("promoter-proximal", prox), ("distal", ~prox)]:
            if m.sum() == 0:
                R(f"    {lab:18s} none"); continue
            row = dict(locus=g, region_class=lab, n=int(m.sum()))
            for col, nm in [("wi38_atac", "WI-38 ATAC"), ("gm21_atac", "GM21 ATAC"),
                            ("gm21_k27ac", "GM21 K27ac"), ("imr90_k27ac", "IMR90 K27ac")]:
                v = rr.loc[m, col].dropna()
                row[f"{col}_n"] = int(len(v))
                row[f"{col}_mean"] = float(v.mean()) if len(v) else None
                row[f"{col}_max_abs"] = float(v.abs().max()) if len(v) else None
            sane.append(row)
            parts = []
            for col, nm in [("wi38_atac", "WI-38 ATAC"), ("gm21_atac", "GM21 ATAC"),
                            ("gm21_k27ac", "GM21 K27ac"), ("imr90_k27ac", "IMR90 K27ac")]:
                v = rr.loc[m, col].dropna()
                parts.append(f"{nm} {('%+.3f' % v.mean()) if len(v) else '  n/a'} "
                             f"(n={len(v)})")
            R(f"    {lab:18s} " + " | ".join(parts))
    pd.DataFrame(sane).to_csv(OUT / "locus_sanity_checks.csv", index=False)

    # ---------------- save ----------------
    full = pd.concat([adf.reset_index(drop=True), resp.reset_index(drop=True)], axis=1)
    full.to_csv(OUT / "region_responses.csv", index=False)
    json.dump({"active_vs_open": act,
               "primary_threshold": PRIMARY_T, "pseudocount": PSEUDO},
              open(OUT / "step2_stats.json", "w"), indent=2)
    (OUT / "STEP2_LOG.txt").write_text("\n".join(LOG))
    R(f"\n  wrote region_responses.csv {full.shape}")
    R("STEP2_DONE")


if __name__ == "__main__":
    main()
