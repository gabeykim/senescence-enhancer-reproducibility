#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 71_ois_trainset.py
#
# Quantify targets, batch PCA, window covering, and emit the OIS AnnData (111,671 regions).
#
# CONSUMES: output/ois_enhancer_trainset/regions_hg19_hg38.csv, signal_matrix_ois.csv
# PRODUCES: output/ois_enhancer_trainset/ois_enhancer_trainset.h5ad, step2_stats.json, window_covering_per_chrom.csv, batch_pca_stats.csv, locus_sanity_checks.csv, STEP2_LOG.txt
# ---------------------------------------------------------------------------
"""
OIS-scoped enhancer training set -- STEP 2: targets, batch check, AnnData.

THE SPLIT: train on GM21 (GSE205898), validate on IMR90 (GSE74238). Cross-lab and
cross-tissue (skin vs lung fibroblast), which no prior stage of this project has had --
earlier work used chromosome holdouts inside a single dataset. chr9 and chr6 are ADDITIONALLY
held out inside the training source so CDKN2A and CDKN1A stay clean benchmark loci.

Note on what the cross-study split does and does not do: the holdout is on DATA SOURCE, not
on regions. The same genomic regions are scored in both studies -- that is the point, since
the question is whether a model fit to GM21's response predicts IMR90's response at the same
elements. Region-level leakage is therefore not the relevant hazard; the relevant hazard is
that region DEFINITION came from GM21 peaks, which is stated and quantified.
"""
import json
from pathlib import Path

import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
OUT = ROOT / "output" / "ois_enhancer_trainset"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SEQ_LEN = 524_288
PRED_SPAN = 196_608          # Borzoi's predicted output span within the 524 kb input
PSEUDO = 0.10
TEST_CHROMS = ["chr9", "chr6"]
ATAC_LFC_THRESH = 0.5        # optional accessibility filter, DEFAULT OFF

N_ARM = {"GM21_RAS_vs_EV": "n=2 RAS vs 2 EV",
         "IMR90_SEN_vs_PRO": "n=2 SEN vs 2 PRO",
         "IMR90_SEN_vs_QUI": "n=2 SEN vs 2 QUI",
         "GM21_ATAC_RAS_vs_EV": "n=4 RAS vs 4 EV"}
LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    rdf = pd.read_csv(OUT / "regions_hg19_hg38.csv")
    sm = pd.read_csv(OUT / "signal_matrix_ois.csv")
    R("=" * 84)
    R("OIS ENHANCER TRAINING SET -- STEP 2")
    R("=" * 84)
    R(f"regions {len(rdf):,} | signal columns {sm.shape[1]}")

    norm = sm.copy()
    for c in norm.columns:
        mu = np.nanmean(norm[c].to_numpy())
        norm[c] = norm[c] / mu if mu and np.isfinite(mu) else np.nan

    def lfc(on, off):
        return np.log2((norm[on].mean(axis=1).to_numpy() + PSEUDO) /
                       (norm[off].mean(axis=1).to_numpy() + PSEUDO))

    R("\n[1] TARGETS (log2 fold change, per dataset, own arms)\n" + "-" * 84)
    tgt = pd.DataFrame(index=rdf.index)
    tgt["GM21_RAS_vs_EV"] = lfc(["GM21_K27_RAS_R1", "GM21_K27_RAS_R2"],
                                ["GM21_K27_EV_R1", "GM21_K27_EV_R2"])
    tgt["IMR90_SEN_vs_PRO"] = lfc(["IMR90_K27_SEN_R1", "IMR90_K27_SEN_R2"],
                                  ["IMR90_K27_PRO_R1", "IMR90_K27_PRO_R2"])
    tgt["IMR90_SEN_vs_QUI"] = lfc(["IMR90_K27_SEN_R1", "IMR90_K27_SEN_R2"],
                                  ["IMR90_K27_QUI_R1", "IMR90_K27_QUI_R2"])
    tgt["GM21_ATAC_RAS_vs_EV"] = lfc(
        ["GM21_ATAC_RAS_R1", "GM21_ATAC_RAS_R2", "GM21_ATAC_RAS_R3", "GM21_ATAC_RAS_R4"],
        ["GM21_ATAC_EV_R1", "GM21_ATAC_EV_R2", "GM21_ATAC_EV_R3", "GM21_ATAC_EV_R4"])
    for c in tgt.columns:
        v = tgt[c].dropna()
        R(f"  {c:22s} {N_ARM[c]:18s} n={len(v):7,}  mean {v.mean():+.4f}  sd {v.std():.4f}  "
          f"range [{v.min():+.2f}, {v.max():+.2f}]")

    # quiescent-vs-proliferating, to show the quiescent arm behaves as growth arrest
    qvp = lfc(["IMR90_K27_QUI_R1", "IMR90_K27_QUI_R2"],
              ["IMR90_K27_PRO_R1", "IMR90_K27_PRO_R2"])
    R(f"\n  supplementary: IMR90 QUI vs PRO  mean {np.nanmean(qvp):+.4f} "
      f"sd {np.nanstd(qvp):.4f}")
    r_sq = stats.spearmanr(tgt["IMR90_SEN_vs_PRO"], tgt["IMR90_SEN_vs_QUI"],
                           nan_policy="omit").statistic
    r_qp = stats.spearmanr(tgt["IMR90_SEN_vs_PRO"], qvp, nan_policy="omit").statistic
    R(f"  Spearman(SEN_vs_PRO, SEN_vs_QUI) = {r_sq:+.4f}   "
      f"Spearman(SEN_vs_PRO, QUI_vs_PRO) = {r_qp:+.4f}")
    R("  A senescence-specific programme should give a HIGH first value (the senescence")
    R("  response is similar whichever arrested control is used) and a LOW second value")
    R("  (quiescence is not just a weaker senescence).")

    # ---------------- reproduce the prior +0.653 ----------------
    R("\n[2] HARMONIZATION CHECK -- does GM21 vs IMR90 reproduce the prior +0.653?\n" + "-" * 84)
    m = tgt[["GM21_RAS_vs_EV", "IMR90_SEN_vs_PRO"]].dropna()
    for t in [0.0, 0.25, 0.5, 1.0]:
        mm = m[(m["GM21_RAS_vs_EV"].abs() >= t) | (m["IMR90_SEN_vs_PRO"].abs() >= t)] if t else m
        sp = stats.spearmanr(mm.iloc[:, 0], mm.iloc[:, 1]).statistic
        pe = stats.pearsonr(mm.iloc[:, 0], mm.iloc[:, 1]).statistic
        star = "  <-- comparable to prior +0.653" if t == 0.5 else ""
        R(f"  |lfc|>={t:<4} n={len(mm):<7,} Spearman {sp:+.4f}  Pearson {pe:+.4f}{star}")
    prior = 0.653
    got = stats.spearmanr(
        *[m[(m['GM21_RAS_vs_EV'].abs() >= .5) | (m['IMR90_SEN_vs_PRO'].abs() >= .5)][c]
          for c in m.columns]).statistic
    R(f"\n  prior test: +{prior:.3f} (1 kb midpoint-centred anchors from combined ATAC+K27ac "
      f"consensus)")
    R(f"  this build: {got:+.3f} (1 kb SUMMIT-centred anchors from K27ac peaks only)")
    R(f"  difference: {got-prior:+.3f}")
    R("  The two region definitions are not identical -- this build centres on H3K27ac")
    R("  summits and drops ATAC-only regions -- so exact equality is not expected. What")
    R("  matters is that the cross-study OIS agreement is reproduced at the same order.")

    # ---------------- batch check ----------------
    R("\n[3] BATCH CHECK\n" + "-" * 84)
    k27 = [c for c in norm.columns if "K27" in c]
    M = np.log2(norm[k27].to_numpy() + PSEUDO)
    keep = np.isfinite(M).all(axis=1)
    M = M[keep]
    meta = pd.DataFrame([{"sample": c,
                          "study": "GSE205898" if c.startswith("GM21") else "GSE74238",
                          "status": ("SEN" if ("RAS" in c or "SEN" in c) else
                                     "QUI" if "QUI" in c else "PRO")} for c in k27])
    R(f"  matrix {M.shape[0]:,} regions x {M.shape[1]} H3K27ac samples")

    def r2(y, lab):
        y = np.asarray(y, float); gt = y.mean(); sst = ((y - gt) ** 2).sum()
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
            a = r2(sc[:, k], meta["study"].to_numpy())
            b = r2(sc[:, k], meta["status"].to_numpy())
            rows.append(dict(space=tag, pc=f"PC{k+1}", var=float(fr[k]),
                             r2_study=a, r2_status=b))
            R(f"  {tag:20s} PC{k+1}  var {100*fr[k]:5.1f}%  R2(study) {a:.3f}  "
              f"R2(status) {b:.3f}  -> {'STUDY' if a > b else 'STATUS'} dominates")
        return sc, fr, rows

    sc_r, fr_r, r1 = pca(M, "raw pooled")
    Mc = M.copy()
    for s in meta["study"].unique():
        mm = (meta["study"] == s).to_numpy()
        Mc[:, mm] = M[:, mm] - M[:, mm].mean(axis=1, keepdims=True)
    R("")
    sc_c, fr_c, r2r = pca(Mc, "per-study centered")
    pd.DataFrame(r1 + r2r).to_csv(OUT / "batch_pca_stats.csv", index=False)

    R("\n  TRADEOFF (batch correction vs the cross-study holdout):")
    R("    The split is cross-study BY DESIGN, so study and validation-fold are perfectly")
    R("    confounded. Any correction fitted across both studies would use validation-set")
    R("    information at training time and would also remove exactly the between-study")
    R("    variation the holdout exists to test.")
    R("  RECOMMENDATION: do NOT apply a pooled batch correction. The targets are already")
    R("    per-study log2 fold changes -- within-study contrasts in which the study-level")
    R("    offset cancels by construction. That is the correct 'correction' here: it is")
    R("    computed inside each study independently and never mixes them. The raw-signal PCA")
    R("    above is reported only to show what pooling raw coverage WOULD have done.")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, (sc, fr, t) in zip(axes, [(sc_r, fr_r, "raw pooled"),
                                      (sc_c, fr_c, "per-study centered")]):
        for s, mk in zip(meta["study"].unique(), ["o", "s"]):
            for st, col in [("PRO", "#1f77b4"), ("SEN", "#d62728"), ("QUI", "#2ca02c")]:
                mm = ((meta["study"] == s) & (meta["status"] == st)).to_numpy()
                if mm.any():
                    ax.scatter(sc[mm, 0], sc[mm, 1], marker=mk, c=col, s=120,
                               edgecolors="k", linewidths=.5,
                               label=f"{s} {st}" if t == "raw pooled" else None)
        ax.set_title(f"{t}\nPC1 {100*fr[0]:.1f}%  PC2 {100*fr[1]:.1f}%")
        ax.axhline(0, lw=.5, c="#ccc"); ax.axvline(0, lw=.5, c="#ccc")
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    axes[0].legend(fontsize=7)
    fig.suptitle("OIS H3K27ac: batch structure (shape = study, colour = status)")
    fig.tight_layout(); fig.savefig(FIG / "ois_batch_pca.png", dpi=150)
    R(f"  figure -> {FIG/'ois_batch_pca.png'}")

    # ---------------- accessibility filter ----------------
    R("\n[4] ACCESSIBILITY FILTER (optional, DEFAULT OFF)\n" + "-" * 84)
    both = tgt[["GM21_RAS_vs_EV", "GM21_ATAC_RAS_vs_EV"]].dropna()
    for t in [0.5, 1.0]:
        dk = both[both["GM21_RAS_vs_EV"].abs() >= t]
        also = (dk["GM21_ATAC_RAS_vs_EV"].abs() >= t)
        same = also & (np.sign(dk["GM21_ATAC_RAS_vs_EV"]) == np.sign(dk["GM21_RAS_vs_EV"]))
        R(f"  |log2FC| >= {t}: {int(also.sum()):,}/{len(dk):,} ({100*also.mean():.1f}%) of "
          f"differential-H3K27ac regions are also differential in ATAC; "
          f"{100*same.mean():.1f}% agree in direction")
    acc_flag = (tgt["GM21_ATAC_RAS_vs_EV"].abs() >= ATAC_LFC_THRESH).fillna(False).to_numpy()
    R(f"  regions passing the ATAC filter (|lfc|>={ATAC_LFC_THRESH}): "
      f"{int(acc_flag.sum()):,}/{len(tgt):,}")
    R(f"  DEFAULT IS OFF -- the flag is stored in .var['passes_atac_filter'] so a training")
    R(f"  run can subset without a rebuild.")

    # ---------------- intervals ----------------
    R("\n[5] BORZOI INTERVALS (hg38, 524,288 bp centred on each region)\n" + "-" * 84)
    sizes = {}
    for line in open(ROOT / "data" / "annotation" / "hg38.chrom.sizes"):
        a, b = line.split()
        sizes[a] = int(b)
    ok = rdf["lifted"].to_numpy()
    mid = rdf["summit_hg38"].to_numpy()
    half = SEQ_LEN // 2
    start = mid - half
    end = start + SEQ_LEN
    csz = np.array([sizes.get(c, 0) for c in rdf["chrom"]])
    valid = ok & (start >= 0) & (end <= csz) & (csz > 0)
    R(f"  regions with a valid 524 kb window: {int(valid.sum()):,}/{len(rdf):,}")
    R(f"    lost -- liftover failed        : {int((~ok).sum()):,}")
    R(f"    lost -- window before chrom 0  : {int((ok & (start < 0)).sum()):,}")
    R(f"    lost -- window past chrom end  : {int((ok & (end > csz) & (csz > 0)).sum()):,}")

    # ---------------- distinct windows (GPU budget) ----------------
    R("\n[6] DISTINCT 524 kb WINDOWS REQUIRED (this sets the GPU budget)\n" + "-" * 84)
    R(f"  Borzoi predicts a {PRED_SPAN:,} bp span from a {SEQ_LEN:,} bp input, so one")
    R(f"  forward pass can serve every region falling inside that central span. Minimum")
    R(f"  windows computed by greedy interval covering per chromosome (measured, not")
    R(f"  estimated).")
    vv = rdf[valid].copy()
    total_win = 0
    per_chrom = []
    for c, g in vv.groupby("chrom"):
        pos = np.sort(g["summit_hg38"].to_numpy())
        n_w, i = 0, 0
        while i < len(pos):
            anchor = pos[i] + PRED_SPAN // 2 - 1        # window covers [pos[i], pos[i]+SPAN)
            n_w += 1
            i = np.searchsorted(pos, pos[i] + PRED_SPAN, side="left")
        total_win += n_w
        per_chrom.append(dict(chrom=c, regions=len(g), windows=n_w))
    R(f"  regions: {int(valid.sum()):,}")
    R(f"  distinct windows needed: {total_win:,}")
    R(f"  compression: {valid.sum()/max(total_win,1):.1f}x fewer forward passes")
    prior_rate = 0.2018   # s/pass measured on the RTX 4090 in the generation probe
    R(f"  at the measured {prior_rate:.4f} s/forward-pass (fwd only): "
      f"{total_win*prior_rate/3600:.2f} h; with reverse-complement too: "
      f"{2*total_win*prior_rate/3600:.2f} h")
    R(f"  naive one-window-per-region would cost "
      f"{2*valid.sum()*prior_rate/3600:.1f} h -- the covering saves "
      f"{2*(valid.sum()-total_win)*prior_rate/3600:.1f} h")
    pd.DataFrame(per_chrom).to_csv(OUT / "window_covering_per_chrom.csv", index=False)

    # ---------------- AnnData ----------------
    R("\n[7] ANNDATA\n" + "-" * 84)
    v = rdf[valid].reset_index(drop=True)
    T = tgt[valid].reset_index(drop=True)
    accv = acc_flag[valid]
    var = pd.DataFrame({
        "chrom": v["chrom"], "start": (v["summit_hg38"] - half).astype(int),
        "end": (v["summit_hg38"] - half + SEQ_LEN).astype(int),
        "region_center_hg38": v["summit_hg38"].astype(int),
        "region_start_hg38": (v["summit_hg38"] - 500).astype(int),
        "region_end_hg38": (v["summit_hg38"] + 500).astype(int),
        "summit_hg19": v["summit_hg19"].astype(int),
        "design_anchor_hg38": v["atac_summit_hg38"].astype(int),
        "design_anchor_dist": v["atac_summit_dist"].astype(int),
        "n_samples_supporting": v["n_samples_supporting"].astype(int),
        "passes_atac_filter": accv,
    })
    var.index = (v["chrom"] + ":" + v["summit_hg38"].astype(str))
    var.index.name = "region"
    var["split"] = np.where(var["chrom"].isin(TEST_CHROMS), "test", "train")

    TASKS = ["GM21_RAS_vs_EV", "IMR90_SEN_vs_PRO", "IMR90_SEN_vs_QUI"]
    X = np.vstack([T[t].to_numpy() for t in TASKS]).astype(np.float32)
    X = np.nan_to_num(X, nan=0.0)
    obs = pd.DataFrame({
        "task": TASKS,
        "study": ["GSE205898", "GSE74238", "GSE74238"],
        "cell_line": ["GM21", "IMR90", "IMR90"],
        "tissue": ["skin fibroblast", "lung fibroblast", "lung fibroblast"],
        "mechanism": ["OIS_ERRAS", "OIS_HRasV12", "OIS_HRasV12"],
        "control_arm": ["empty vector", "proliferating", "quiescent"],
        "n_replicates": ["2 vs 2", "2 vs 2", "2 vs 2"],
        "role": ["train", "validation", "validation"],
    }).set_index("task")

    ad = anndata.AnnData(X=X, obs=obs, var=var)
    ad.layers  # noqa
    ad.uns["contract"] = {
        "seq_len": SEQ_LEN, "bin_size": 32, "label_len": PRED_SPAN, "genome": "hg38",
        "input_shape": "(4, 524288) float32 channels-first",
        "trunk_output": "(batch, 1920, 6144)",
        "target_per_example": "(n_tasks, 1) float32",
        "head": "ConvHead act_func=None (targets are SIGNED log2 fold changes; a softplus "
                "head cannot represent the negative half)",
    }
    ad.uns["target_transform"] = ("per-region log2 fold change of H3K27ac, computed WITHIN "
                                 "each study (mean-1 normalized bigWig coverage over 1 kb "
                                 f"summit-centred windows, pseudocount {PSEUDO})")
    ad.uns["split_strategy"] = (
        "CROSS-STUDY: train on GM21 (GSE205898), validate on IMR90 (GSE74238) -- cross-lab "
        "and cross-tissue (skin vs lung fibroblast). chr9 and chr6 are additionally held out "
        "within the training source so CDKN2A and CDKN1A remain clean benchmark loci. The "
        "holdout is on DATA SOURCE, not on regions: the same regions are scored in both "
        "studies by design.")
    ad.uns["scope"] = ("OIS only. Replicative senescence excluded: within-replicative H3K27ac "
                       "reproducibility +0.054/+0.144/+0.268 vs OIS +0.653.")
    ad.uns["accessibility_filter"] = {
        "default": "OFF", "threshold": ATAC_LFC_THRESH,
        "column": "var['passes_atac_filter']",
        "n_passing": int(accv.sum()), "n_total": int(len(var))}
    ad.uns["window_covering"] = {"regions": int(valid.sum()), "distinct_windows": total_win,
                                 "predicted_span": PRED_SPAN}
    ad.uns["provenance"] = {"built": pd.Timestamp.today().strftime("%Y-%m-%d"),
                            "scripts": "70_build_ois_regions.py; 71_ois_trainset.py"}
    p = OUT / "ois_enhancer_trainset.h5ad"
    ad.write_h5ad(p)
    R(f"  wrote {p.name}  shape {ad.shape} (n_tasks, n_intervals)")
    R(f"  X dtype {ad.X.dtype}  range [{ad.X.min():.3f}, {ad.X.max():.3f}]")
    R(f"  split counts: {var['split'].value_counts().to_dict()}")
    R(f"  interval widths distinct: {set((var['end']-var['start']).unique())}")

    # ---------------- sanity ----------------
    R("\n[8] SANITY CHECKS\n" + "-" * 84)
    full = pd.concat([v.reset_index(drop=True), T.reset_index(drop=True)], axis=1)
    full["passes_atac_filter"] = accv
    full.to_csv(OUT / "region_responses_ois.csv", index=False)

    LOCI = {"CDKN2A": ("chr9", 21_967_751, 21_995_300),
            "CDKN1A": ("chr6", 36_644_237, 36_655_116)}
    srows = []
    for g, (c, s, e) in LOCI.items():
        selm = ((full["chrom"] == c) & (full["summit_hg19"] >= s - 200_000) &
                (full["summit_hg19"] <= e + 200_000)).to_numpy()
        d = full[selm]
        dist = np.where(d["summit_hg19"] < s, s - d["summit_hg19"],
                        np.where(d["summit_hg19"] > e, d["summit_hg19"] - e, 0))
        prox = dist <= 2000
        R(f"\n  {g}: {int(selm.sum())} regions (promoter<=2kb {int(prox.sum())}, "
          f"distal {int((~prox).sum())})")
        for lab, mm in [("promoter", prox), ("distal", ~prox)]:
            if mm.sum() == 0:
                R(f"    {lab:9s} none"); continue
            row = dict(locus=g, region_class=lab, n=int(mm.sum()))
            parts = []
            for t in ["GM21_RAS_vs_EV", "IMR90_SEN_vs_PRO", "IMR90_SEN_vs_QUI"]:
                x = d.loc[mm, t].dropna()
                row[t] = float(x.mean()) if len(x) else None
                parts.append(f"{t.split('_')[0]} {('%+.3f' % x.mean()) if len(x) else 'n/a'}")
            srows.append(row)
            R(f"    {lab:9s} " + " | ".join(parts))
    pd.DataFrame(srows).to_csv(OUT / "locus_sanity_checks.csv", index=False)

    R("\n  SASP loci (regions within 50 kb of the gene body, hg19):")
    SASP = {"IL6": ("chr7", 22_766_766, 22_771_621),
            "CXCL8": ("chr4", 74_606_223, 74_609_433),
            "IL1A": ("chr2", 113_531_492, 113_542_971),
            "MMP3": ("chr11", 102_706_532, 102_714_342)}
    sasp_rows = []
    for g, (c, s, e) in SASP.items():
        mm = ((full["chrom"] == c) & (full["summit_hg19"] >= s - 50_000) &
              (full["summit_hg19"] <= e + 50_000)).to_numpy()
        if mm.sum() == 0:
            R(f"    {g:6s} no regions in window"); continue
        d = full[mm]
        row = dict(gene=g, n_regions=int(mm.sum()))
        parts = []
        for t in ["GM21_RAS_vs_EV", "IMR90_SEN_vs_PRO", "IMR90_SEN_vs_QUI"]:
            x = d[t].dropna()
            row[t + "_mean"] = float(x.mean()) if len(x) else None
            row[t + "_max"] = float(x.max()) if len(x) else None
            parts.append(f"{t.split('_')[0]} mean {('%+.3f' % x.mean()) if len(x) else 'n/a'} "
                         f"max {('%+.3f' % x.max()) if len(x) else 'n/a'}")
        sasp_rows.append(row)
        R(f"    {g:6s} n={int(mm.sum()):3d}  " + " | ".join(parts))
    pd.DataFrame(sasp_rows).to_csv(OUT / "sasp_checks.csv", index=False)

    R("\n  Split integrity:")
    R(f"    chr9/chr6 regions marked 'test': "
      f"{int((var['split']=='test').sum()):,}; none are used for training.")
    R(f"    training task rows: {list(obs[obs['role']=='train'].index)}")
    R(f"    validation task rows: {list(obs[obs['role']=='validation'].index)}")
    R("    The cross-study holdout is on DATA SOURCE; the same regions appear under both")
    R("    training and validation tasks by design, so region-level overlap is expected and")
    R("    is not leakage. Sample-level overlap between studies: none (disjoint GSMs).")

    json.dump({"pseudocount": PSEUDO, "test_chroms": TEST_CHROMS,
               "atac_threshold": ATAC_LFC_THRESH,
               "spearman_gm21_vs_imr90_primary": float(got),
               "distinct_windows": total_win, "regions": int(valid.sum())},
              open(OUT / "step2_stats.json", "w"), indent=2)
    (OUT / "STEP2_LOG.txt").write_text("\n".join(LOG))
    R("\nSTEP2_DONE")


if __name__ == "__main__":
    main()
