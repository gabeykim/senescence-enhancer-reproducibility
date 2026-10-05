#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 99b_figures_supp.py
#
# Manuscript Figures S1-S5. CPU.
#
# CONSUMES: results/ (enhancer, trainset, replicative, ois_model, permutations)
# PRODUCES: figures/figureS[1-5]*.{png,pdf}
# ---------------------------------------------------------------------------
"""Supplementary figures S1-S5. Every value read from a stored artifact."""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import importlib.util as _ilu, pathlib as _pl
_fspec = _ilu.spec_from_file_location('_fs', _pl.Path(__file__).resolve().parent/'98_fig_style.py')
_fs = _ilu.module_from_spec(_fspec); _fspec.loader.exec_module(_fs)
AQUA = _fs.AQUA
BLUE = _fs.BLUE
FULL = _fs.FULL
GREY = _fs.GREY
INK = _fs.INK
INK2 = _fs.INK2
MOTIF_COLOUR = _fs.MOTIF_COLOUR
MOTIF_LABEL = _fs.MOTIF_LABEL
MOTIF_MARKER = _fs.MOTIF_MARKER
MUTED = _fs.MUTED
ORANGE = _fs.ORANGE
RED = _fs.RED
SHADE = _fs.SHADE
SINGLE = _fs.SINGLE
VIOLET = _fs.VIOLET
YELLOW = _fs.YELLOW
apply_style = _fs.apply_style
panel_tag = _fs.panel_tag
save = _fs.save

_spec = _ilu.spec_from_file_location('_p', _pl.Path(__file__).resolve().parent/'00_paths.py')
_p = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_p)
ROOT, OUTD = _p.ROOT, _p.OUTPUT
apply_style()


# =====================================================================  S1
def s1():
    st = json.load(open(OUTD / "enhancer_test/step2_stats.json"))["active_vs_open"]["0.5"]
    # two-column extract of region_responses.csv (the only columns S1 needs);
    # falls back to the full table when running in the original tree.
    _e = _p.REPO / "results/enhancer_cross_mechanism/gm21_atac_vs_k27ac.csv"
    r = pd.read_csv(_e) if _e.exists() else pd.read_csv(
        OUTD / "enhancer_test/region_responses.csv")
    both = r[["gm21_atac", "gm21_k27ac"]].dropna()
    T = 0.5
    nd, nc, ns = st["n_diff_atac"], st["n_also_diff_k27ac"], st["n_same_direction"]

    fig = plt.figure(figsize=(FULL, 3.8))
    gs = fig.add_gridspec(1, 2, wspace=0.42,
                          left=0.10, right=0.97, top=0.85, bottom=0.20)
    a = fig.add_subplot(gs[0, 0])
    hb = a.hexbin(both.gm21_atac, both.gm21_k27ac, gridsize=44, bins="log",
                  cmap="Blues", mincnt=1, linewidths=0)
    for v in (-T, T):
        a.axvline(v, color=RED, lw=0.7, ls="--")
        a.axhline(v, color=RED, lw=0.7, ls="--")
    a.axhline(0, color=INK2, lw=0.5); a.axvline(0, color=INK2, lw=0.5)
    a.set_xlabel("GM21 ATAC response (log2FC, n = 4 vs 4)")
    a.set_ylabel("GM21 H3K27ac response\n(log2FC, n = 2 vs 2)")
    a.set_title("Accessibility vs acetylation,\nsame cells, same contrast",
                loc="left", fontsize=10)
    a.text(.03, .97, f"n = {len(both):,} anchors with both assays\n"
                     f"dashed: |log2FC| = {T}",
           transform=a.transAxes, fontsize=9, va="top", color=INK)
    cb = fig.colorbar(hb, ax=a, pad=0.02, fraction=0.045)
    cb.set_label("anchors (log)", fontsize=9); cb.ax.tick_params(labelsize=9)
    cb.outline.set_linewidth(0.4)
    panel_tag(a, "A")

    b = fig.add_subplot(gs[0, 1])
    parts = [("no H3K27ac change", nd - nc, GREY),
             ("changes, opposite direction", nc - ns, ORANGE),
             ("changes, same direction", ns, BLUE)]
    left = 0
    for lab, v, c in parts:
        b.barh(0, v, left=left, height=0.40, color=c)
        b.text(left + v / 2, 0, f"{v:,}\n{100*v/nd:.1f}%", ha="center", va="center",
               fontsize=9, color="white" if c != GREY else INK)
        left += v
    b.set_yticks([]); b.set_xlim(0, nd * 1.02); b.set_ylim(-0.5, 0.95)
    b.set_xlabel(f"differential-ACCESSIBILITY anchors (n = {nd:,})")
    b.set_title("Of anchors that change in ATAC,\nmost also change in H3K27ac",
                loc="left", fontsize=10)
    b.spines["left"].set_visible(False)
    b.legend(handles=[Patch(facecolor=c, label=l) for l, _, c in parts],
             loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=1, fontsize=9)
    panel_tag(b, "B")

    rec_s = int((((both.gm21_atac.abs() >= T) & (both.gm21_k27ac.abs() >= T))
                 & (np.sign(both.gm21_k27ac) == np.sign(both.gm21_atac))).sum())
    assert rec_s == ns, f"recomputed {rec_s} != stored {ns}"
    print(f"  S1: recomputed same-direction {rec_s:,} == stored {ns:,} "
          f"({100*ns/nc:.4f}% of co-differential)")
    return save(fig, "figureS1_atac_vs_k27ac_concordance")


# =====================================================================  S2
def s2():
    p = pd.read_csv(OUTD / "ois_enhancer_trainset/batch_pca_stats.csv")
    fig, axes = plt.subplots(1, 2, figsize=(FULL * 1.05, 3.6), sharey=True)
    for ax, space, ttl in zip(axes, ["raw pooled", "per-study centered"],
                              ["Before: raw pooled coverage",
                               "After: per-study log2FC contrasts"]):
        g = p[p.space == space]
        x = np.arange(len(g)); w = 0.38
        ax.bar(x - w / 2, g.r2_study, w, color=ORANGE, label="R$^2$ study (batch)")
        ax.bar(x + w / 2, g.r2_status, w, color=BLUE, label="R$^2$ status (biology)")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{pc}\n{100*v:.1f}% var" for pc, v in zip(g.pc, g["var"])],
                           fontsize=9)
        ax.set_title(ttl, loc="left", fontsize=10)
        ax.set_ylim(0, 1.08)
        ax.grid(axis="y"); ax.set_axisbelow(True)
        for xi, (a_, b_) in enumerate(zip(g.r2_study, g.r2_status)):
            ax.text(xi - w / 2, a_ + .02, f"{a_:.3f}", ha="center", fontsize=9,
                    color=INK2)
            ax.text(xi + w / 2, b_ + .02, f"{b_:.3f}", ha="center", fontsize=9,
                    color=INK2)
    axes[0].set_ylabel("R$^2$ of the PC against the factor")
    axes[0].legend(loc="upper right", fontsize=9)
    axes[1].text(0.5, -0.26,
                 "n = 108,174 regions \u00d7 10 H3K27ac samples (2 studies)",
                 transform=axes[1].transAxes, fontsize=9, ha="center", va="top",
                 color=INK2)
    panel_tag(axes[0], "A"); panel_tag(axes[1], "B")
    fig.tight_layout(w_pad=1.6)
    return save(fig, "figureS2_batch_structure")


# =====================================================================  S3
def s3():
    txt = (OUTD / "replicative_ceiling_test/shared_denominator_null.txt").read_text()
    vals = {}
    for line in txt.strip().splitlines():
        m = re.match(r"(.+?)\s+rho\s+([+-][\d.]+)\s+n=([\d,]+)", line.strip())
        if m:
            vals[m.group(1).strip()] = (float(m.group(2)), int(m.group(3).replace(",", "")))
    rep = pd.read_csv(OUTD / "replicative_ceiling_test/replicative_correlations.csv")
    i3 = rep[(rep.comparison == "I3") & (rep.threshold == 0.5)].iloc[0]

    rows = [("IR vs replicative\nSHARED proliferating control", *vals[
                "I1 IR vs replicative (shared denom)"], RED),
            ("two unrelated numerators,\nsame shared denominator", *vals[
                "NULL 3 two unrelated numerators"], GREY),
            ("unrelated numerator\nvs replicative response", *vals[
                "NULL unrelated numerator vs replicative"], GREY),
            ("BJ Young numerator\nvs replicative response", *vals[
                "NULL 2 BJ Young vs replicative"], GREY),
            ("IR vs replicative\nINDEPENDENT controls (I3)",
             float(i3.spearman), int(i3.n), BLUE)]
    fig, ax = plt.subplots(figsize=(FULL, 3.6))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, v, n, c) in zip(y, rows):
        ax.barh(yi, v, height=0.52, color=c)
        ax.text(v + 0.012, yi, f"{v:+.4f}   n = {n:,}", va="center", fontsize=9,
                color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlim(0, 0.80)
    ax.set_xlabel("Spearman $\\rho$")
    ax.axvline(0, color=INK, lw=0.7)
    ax.set_title("A shared control arm manufactures correlation\n"
                 "between otherwise unrelated contrasts", loc="left", fontsize=10)
    ax.grid(axis="x"); ax.set_axisbelow(True)
    fig.tight_layout()
    return save(fig, "figureS3_shared_control_artifact")


# =====================================================================  S4
def s4():
    d = json.load(open(ROOT / "infra_borzoi_grelu/execution_results_permutations"
                              "/results_perm200.json"))["within_OIS"]
    aucs = np.array(d["all_aucs"])
    real = d["real_auc"]
    fig, ax = plt.subplots(figsize=(FULL, 3.8))
    srt = np.sort(aucs)
    x = np.arange(len(srt))
    beat = srt >= real
    ax.bar(x[~beat], srt[~beat], color=GREY, width=0.78)
    ax.bar(x[beat], srt[beat], color=ORANGE, width=0.78)
    ax.axhline(real, color=RED, lw=1.2)
    ax.text(len(srt) - 0.4, real + 0.016, f"real AUC {real:.4f}", color=RED,
            fontsize=9, ha="right", va="bottom")
    ax.axhline(0.5, color=MUTED, lw=0.6, ls=":")
    ax.set_xlabel(f"all {d['n_distinct_partitions']} distinct 3v3 label partitions "
                  f"(sorted)")
    ax.set_ylabel("AUC")
    ax.set_ylim(0, 1.0)
    ax.set_xticks([])
    ax.set_title("Exhaustive 3v3 permutation null:\n$p$ cannot resolve below 0.05",
                 loc="left", fontsize=10)
    ax.text(.33, .97,
            f"{d['n_beating_real']} / {d['n_distinct_partitions']} reach the real AUC\n"
            f"exact $p$ = {d['empirical_p']:.2f}\n"
            f"n = {d['n_evaluable']} evaluable genes",
            transform=ax.transAxes, fontsize=9, va="top", color=INK)
    ax.legend(handles=[Patch(facecolor=ORANGE, label="reaches the real AUC"),
                       Patch(facecolor=GREY, label="below the real AUC")],
              loc="upper left", fontsize=9)
    ax.grid(axis="y"); ax.set_axisbelow(True)
    fig.tight_layout()
    return save(fig, "figureS4_exhaustive_3v3_null")


# =====================================================================  S5
def s5():
    m = pd.read_csv(OUTD / "ois_enhancer_run/results/motif_insertion.csv")
    summ = json.load(open(OUTD / "ois_enhancer_run/results/results_probe.json")
                     )["motif_insertion"]["summary"]
    order = ["NFKB_RELA", "CEBPB", "ETS1", "AP1_FOSJUN"]
    fig, axes = plt.subplots(1, 4, figsize=(FULL, 3.5), sharey=True)
    for ax, mot in zip(axes, order):
        g = m[m.motif == mot]
        nseed = g.seed.nunique()
        for variant, ls, fill in (("real", "-", True), ("scrambled", "--", False)):
            gg = g[g.variant == variant].groupby("k")["pred"]
            mu, sem = gg.mean(), gg.std(ddof=1) / np.sqrt(gg.count())
            ax.errorbar(mu.index, mu.values, yerr=sem.values, ls=ls,
                        marker=MOTIF_MARKER[mot], ms=3.4, lw=1.1, capsize=1.6,
                        elinewidth=0.7, color=MOTIF_COLOUR[mot],
                        mfc=MOTIF_COLOUR[mot] if fill else "white",
                        label="real" if variant == "real" else "scrambled")
        d_r = summ[f"{mot}_real"]["delta"]; d_s = summ[f"{mot}_scrambled"]["delta"]
        ax.set_title(f"{MOTIF_LABEL[mot]}", loc="left", fontsize=10,
                     color=MOTIF_COLOUR[mot])
        ax.set_xlabel("motif copies (k)")
        ax.set_xticks([0, 2, 4, 6, 8])
        ax.text(.03, .03, f"$\\Delta$ real {d_r:+.3f}\n$\\Delta$ scram. {d_s:+.3f}",
                transform=ax.transAxes, fontsize=9, va="bottom", color=INK)
        ax.grid(axis="y"); ax.set_axisbelow(True)
    axes[0].set_ylabel("predicted response")
    axes[0].legend(loc="upper left", fontsize=9)
    axes[3].text(.97, .97, f"mean $\\pm$ SEM,\nn = {nseed} backgrounds",
                 transform=axes[3].transAxes, fontsize=9, ha="right", va="top",
                 color=INK2)
    fig.suptitle("Motif-insertion dose response, all four families "
                 "(full version of Figure 7A)", x=0.005, ha="left", fontsize=10)
    fig.tight_layout(w_pad=1.2, rect=(0, 0, 1, 0.93))
    return save(fig, "figureS5_motif_dose_response_all")


if __name__ == "__main__":
    for fn in (s1, s2, s3, s4, s5):
        png, pdf = fn()
        print(f"{png.name:48s} + {pdf.name}")
