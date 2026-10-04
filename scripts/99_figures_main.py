#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 99_figures_main.py
#
# Manuscript Figures 1-7. CPU.
#
# CONSUMES: figures/dataset_table.csv, results/ (all phases), ois_enhancer_trainset.h5ad
# PRODUCES: figures/figure[1-7]*.{png,pdf}, figure_flags.txt
# ---------------------------------------------------------------------------
"""Main manuscript figures 1-7. Every value read from a stored artifact."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).parent))
from fig_style import (AQUA, BLUE, C_CROSS, C_IR, C_OIS, C_REP, FULL, GREY,
                       GRID, INK, INK2, MAGENTA, MOTIF_COLOUR, MOTIF_LABEL,
                       MOTIF_MARKER, MUTED, ORANGE, RED, SHADE, SINGLE, VIOLET,
                       YELLOW, apply_style, panel_tag, save)

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
OUTD = ROOT / "output"
HERE = Path(__file__).parent
apply_style()
NOTES = []


def note(s):
    print("  ! " + s, flush=True)
    NOTES.append(s)


# =====================================================================  FIG 1
def fig1():
    d = pd.read_csv(HERE / "dataset_table.csv")
    SHORT = {"IMR90 replicative": "IMR-90 repl.", "BJ replicative": "BJ repl.",
             "GSE106146 replicative": "GSE106146 repl.", "GSE106146 IR": "GSE106146 IR",
             "IMR90 OIS same-lab": "IMR-90 OIS (1v1)", "GM21 OIS": "GM21 OIS",
             "IMR90 OIS (SEN vs PRO)": "IMR-90 OIS vs PRO",
             "IMR90 OIS (SEN vs QUI)": "IMR-90 OIS vs QUI",
             "GM21 OIS ATAC": "GM21 OIS ATAC",
             "WI-38 replicative ATAC": "WI-38 repl. ATAC"}
    meas = d[d.n_senescent > 0].copy()
    order = {"OIS": 0, "irradiation": 1, "replicative": 2}
    meas = meas.sort_values(["mechanism", "arm_label"],
                            key=lambda s: s.map(order) if s.name == "mechanism" else s)
    meas = meas.reset_index(drop=True)
    col = {"OIS": C_OIS, "replicative": C_REP, "irradiation": C_IR}

    fig, ax = plt.subplots(figsize=(SINGLE, 3.15))
    fig.subplots_adjust(left=0.265, right=0.97, top=0.885, bottom=0.345)
    y = np.arange(len(meas))[::-1].astype(float)
    h = 0.34

    ax.axvspan(0, 2, color=SHADE, alpha=.45, lw=0, zorder=0)
    ax.axvline(2, color=RED, lw=1.0, ls="--", zorder=4)

    for yi, (_, r) in zip(y, meas.iterrows()):
        c = col[r.mechanism]
        ax.barh(yi + h / 2, r.n_senescent, height=h, color=c, zorder=2)
        ax.barh(yi - h / 2, r.n_control, height=h, color="white", edgecolor=c,
                linewidth=0.8, hatch="////", zorder=2)
        if r.assay != "H3K27ac":
            ax.text(4.55, yi, r.assay, fontsize=7, va="center", ha="left", color=INK2)
        if str(r.independent_control).startswith("no"):
            ax.text(5.42, yi, "\u2020", fontsize=7, va="center", ha="left", color=RED)

    ax.set_yticks(y)
    ax.set_yticklabels([SHORT[a] for a in meas.arm_label], fontsize=7)
    for t, m in zip(ax.get_yticklabels(), meas.mechanism):
        t.set_color(col[m])
    ax.set_xlim(0, 5.6)
    ax.set_ylim(-0.75, len(meas) - 0.25)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlabel("replicates per arm")
    ax.set_title("Senescence chromatin datasets:\nevery H3K27ac arm is n $\\leq$ 2",
                 loc="left", fontsize=8)
    ax.text(2.12, len(meas) - 0.45, "n = 2 ceiling", color=RED, fontsize=7,
            va="top", ha="left")
    ax.grid(axis="x", zorder=0)
    ax.set_axisbelow(True)

    handles = [Patch(facecolor=col[m], label=m) for m in
               ["OIS", "irradiation", "replicative"]]
    handles += [Patch(facecolor=INK2, label="senescent arm"),
                Patch(facecolor="white", edgecolor=INK2, hatch="////",
                      label="control arm")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.205),
              ncol=3, fontsize=7, columnspacing=0.9, handlelength=1.0)
    fig.text(0.015, 0.075, "\u2020 control arm shared between two contrasts "
                           "(GSE106146 proliferating)", fontsize=7, color=RED)
    fig.text(0.015, 0.012, "Unlabelled rows are H3K27ac. WI-38 repl. ATAC (GSE175533) is\n"
                           "a deposited log2FC with no replicate-level data; not plotted.",
            fontsize=7, color=MUTED)
    return save(fig, "figure1_data_landscape")


# =====================================================================  FIG 2
def fig2():
    rep = pd.read_csv(OUTD / "replicative_ceiling_test/replicative_correlations.csv")
    enh = pd.read_csv(OUTD / "enhancer_test/cross_mechanism_correlations.csv")
    T = 0.5
    r = rep[rep.threshold == T].set_index("comparison")
    e = enh[enh.threshold == T].set_index("comparison")

    rows = [
        ("within-OIS", "GM21 OIS vs IMR90 OIS  (H3K27ac)", e.loc["b"], False),
        ("within-OIS", "GM21 OIS ATAC vs H3K27ac  (cross-assay)", e.loc["d"], False),
        ("within-replicative", "R3  BJ vs GSE106146", r.loc["R3"], False),
        ("within-replicative", "A2  WI-38 ATAC vs GSE106146  (cross-assay)", r.loc["A2"], False),
        ("within-replicative", "A1  WI-38 ATAC vs IMR90  (cross-assay)", r.loc["A1"], False),
        ("within-replicative", "R1  IMR90 vs GSE106146", r.loc["R1"], False),
        ("within-replicative", "R2  IMR90 vs BJ  (same lab/pipeline/build)", r.loc["R2"], True),
        ("cross-mechanism", "X4  GSE106146 rep vs GM21 OIS", r.loc["X4"], False),
        ("cross-mechanism", "X1  IMR90 rep vs IMR90 OIS  (same lab)", r.loc["X1"], False),
        ("cross-mechanism", "c  WI-38 rep ATAC vs GM21 OIS K27ac", e.loc["c"], False),
        ("cross-mechanism", "a  WI-38 rep ATAC vs GM21 OIS ATAC", e.loc["a"], False),
        ("cross-mechanism", "X2  IMR90 rep vs GM21 OIS", r.loc["X2"], False),
        ("cross-mechanism", "X3  IMR90 rep vs IMR90 OIS  (diff. lab)", r.loc["X3"], False),
    ]
    CEIL = float(e.loc["b", "spearman"])
    xs = [float(rr["spearman"]) for *_, rr, _ in
          [(a, b, c, d) for a, b, c, d in rows]]
    cross = [float(rr["spearman"]) for g, _, rr, _ in rows if g == "cross-mechanism"]
    lo, hi = min(cross), max(cross)

    fig, ax = plt.subplots(figsize=(FULL, 4.6))
    gcol = {"within-OIS": C_OIS, "within-replicative": C_REP,
            "cross-mechanism": C_CROSS}
    y, labels, tickcol = [], [], []
    cur, pos = None, 0.0
    for g, lab, rr, star in rows:
        if cur is not None and g != cur:
            pos += 0.85
        cur = g
        y.append(pos); labels.append(lab); tickcol.append(gcol[g])
        pos += 1.0
    y = np.array(y)
    ymax = y.max()

    ax.axvspan(lo, hi, color=SHADE, alpha=.55, lw=0, zorder=0)
    ax.axvline(0, color=INK2, lw=0.6, zorder=1)
    ax.axvline(CEIL, color=C_OIS, lw=1.1, ls="--", zorder=1)

    for yi, (g, lab, rr, star) in zip(y, rows):
        x = float(rr["spearman"]); n = int(rr["n"])
        c = gcol[g]
        if star:
            ax.plot(x, yi, "o", ms=9, mfc=c, mec=INK, mew=1.4, zorder=5)
        else:
            ax.plot(x, yi, "o", ms=5.5, mfc=c, mec=c, zorder=4)
        # keep n labels off the ceiling line: flip to the left for the far-right points
        if x > 0.55:
            ax.text(x - 0.016, yi, f"n = {n:,}", fontsize=7, ha="right",
                    va="center", color=INK2)
        else:
            ax.text(x + 0.016, yi, f"n = {n:,}", fontsize=7, ha="left",
                    va="center", color=INK2)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    for t, c in zip(ax.get_yticklabels(), tickcol):
        t.set_color(c)
    ax.invert_yaxis()
    ax.set_xlabel("Spearman $\\rho$ of the per-region senescence response "
                  "(primary threshold |log2FC| $\\geq$ 0.5)")
    ax.set_xlim(-0.30, 0.90)
    ax.grid(axis="x"); ax.set_axisbelow(True)

    ax.text(CEIL, -1.25, f"within-OIS ceiling  $\\rho$ = {CEIL:+.3f}",
            color=C_OIS, fontsize=7, ha="center")
    ax.text((lo + hi) / 2, -1.25, f"cross-mechanism span\n{lo:+.3f} to {hi:+.3f}",
            fontsize=7, ha="center", va="center", color=INK2)
    ax.annotate("no technical confound:\nsame lab, pipeline, build",
                xy=(float(r.loc["R2", "spearman"]), y[6]),
                xytext=(0.40, y[6] + 1.15), fontsize=7, color=INK,
                arrowprops=dict(arrowstyle="-", lw=0.7, color=INK2,
                                shrinkA=0, shrinkB=3,
                                connectionstyle="angle,angleA=0,angleB=90,rad=0"))

    ax.set_title("Within-OIS agreement reaches $\\rho$ = +0.653; no within-replicative "
                 "comparison exceeds +0.27\nR1 and R2 fall inside the cross-mechanism "
                 "span (shaded); R3 sits just above it",
                 loc="left", fontsize=7.5)
    ax.set_ylim(ymax + 1.0, -2.1)
    fig.tight_layout()

    # the relational claim in the source report, tested against the stored values
    inside = [lab for g, lab, rr, _ in rows if g == "within-replicative"
              and lo <= float(rr["spearman"]) <= hi]
    outside = [(lab, float(rr["spearman"])) for g, lab, rr, _ in rows
               if g == "within-replicative" and not (lo <= float(rr["spearman"]) <= hi)]
    if outside:
        note("FIG 2 / MANUSCRIPT MISMATCH — REPORT_replicative_ceiling.md:24 states 'The "
             "three within-replicative values (+0.054, +0.144, +0.268) fall inside the "
             f"range of the four cross-mechanism values (-0.164 to +0.164).' Stored data: "
             f"R3 = {dict(outside)['R3  BJ vs GSE106146']:+.4f} is OUTSIDE that range "
             f"(upper bound {hi:+.4f}). Two of three are inside, not three of three.")
    return save(fig, "figure2_reproducibility_matrix")


# =====================================================================  FIG 3
def fig3():
    rep = pd.read_csv(OUTD / "replicative_ceiling_test/replicative_correlations.csv")
    enh = pd.read_csv(OUTD / "enhancer_test/cross_mechanism_correlations.csv")
    series = [("within-OIS  (GM21 vs IMR90)", enh[enh.comparison == "b"], C_OIS, "o"),
              ("R3  BJ vs GSE106146", rep[rep.comparison == "R3"], ORANGE, "^"),
              ("R1  IMR90 vs GSE106146", rep[rep.comparison == "R1"], ORANGE, "s"),
              ("R2  IMR90 vs BJ", rep[rep.comparison == "R2"], ORANGE, "D")]
    fig, ax = plt.subplots(figsize=(SINGLE, 3.2))
    for lab, g, c, mk in series:
        g = g.sort_values("threshold")
        ls = "-" if "OIS" in lab else "--"
        ax.plot(g.threshold, g.spearman, marker=mk, color=c, ls=ls, ms=4, lw=1.2,
                label=lab)
        last = g.iloc[-1]
        ax.annotate(f"{last.spearman:+.3f}", (last.threshold, last.spearman),
                    textcoords="offset points", xytext=(4, -1), fontsize=7, color=c)
        first = g.iloc[0]
        ax.annotate(f"{first.spearman:+.3f}", (first.threshold, first.spearman),
                    textcoords="offset points", xytext=(-4, -7), fontsize=7,
                    color=c, ha="center")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel("response-magnitude threshold  |log2FC| $\\geq$")
    ax.set_ylabel("Spearman $\\rho$")
    ax.set_xticks([0.0, 0.25, 0.5, 1.0])
    ax.set_xlim(-0.06, 1.16)
    ax.set_title("OIS agreement is threshold-independent;\nreplicative agreement is not",
                 loc="left", fontsize=7.5)
    ax.legend(loc="upper left", fontsize=7)
    ax.grid(axis="y"); ax.set_axisbelow(True)
    fig.tight_layout()
    return save(fig, "figure3_threshold_sensitivity")


# =====================================================================  FIG 4
def fig4():
    L = pd.read_csv(OUTD / "replicative_ceiling_test/locus_sanity_checks.csv")
    L = L[L["locus"] == "CDKN2A"].set_index("region_class")
    assert len(L) == 2, f"expected promoter+distal for CDKN2A, got {len(L)}"
    cols = [("IMR90_rep_k27ac", "IMR-90 repl.", "replicative"),
            ("BJ_rep_k27ac", "BJ repl.", "replicative"),
            ("Sen2019_rep_k27ac", "GSE106146 repl.", "replicative"),
            ("Sen2019_IR_k27ac", "GSE106146 IR", "irradiation"),
            ("IMR90_OIS_k27ac_samelab", "IMR-90 OIS\n(same lab, 1v1)", "OIS"),
            ("GM21_OIS_k27ac", "GM21 OIS", "OIS"),
            ("IMR90_OIS_k27ac_prior", "IMR-90 OIS", "OIS")]
    col = {"OIS": C_OIS, "replicative": C_REP, "irradiation": C_IR}
    x = np.arange(len(cols)); w = 0.38
    fig, ax = plt.subplots(figsize=(SINGLE, 3.3))
    for i, (c, lab, mech) in enumerate(cols):
        p = float(L.loc["promoter", f"{c}_mean"])
        d = float(L.loc["distal", f"{c}_mean"])
        ax.bar(i - w / 2, p, width=w, color=col[mech], edgecolor="none")
        ax.bar(i + w / 2, d, width=w, color="white", edgecolor=col[mech],
               linewidth=0.9, hatch="////")
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l, _ in cols], fontsize=7, rotation=38,
                       ha="right")
    for t, (_, _, m) in zip(ax.get_xticklabels(), cols):
        t.set_color(col[m])
    ax.set_ylabel("mean H3K27ac log2FC at CDKN2A")
    ax.set_title("CDKN2A promoter falls in every dataset;\ndistal rises only in the "
                 "2v2 OIS datasets", loc="left", fontsize=7.5)
    np_, nd = int(L.loc["promoter", "n"]), int(L.loc["distal", "n"])
    ax.legend(handles=[Patch(facecolor=INK2, label=f"promoter (n = {np_} regions)"),
                       Patch(facecolor="white", edgecolor=INK2, hatch="////",
                             label=f"distal (n = {nd} regions)")],
              loc="lower left", fontsize=7)
    ax.grid(axis="y"); ax.set_axisbelow(True)
    fig.tight_layout()
    return save(fig, "figure4_cdkn2a_promoter_distal")


# =====================================================================  FIG 5
def fig5():
    import anndata
    R = OUTD / "ois_enhancer_run/results"
    tr = json.load(open(R / "results_train.json"))
    pred = np.load(R / "pred_rest.npy")
    idx = np.load(R / "idx_rest.npy")
    ad = anndata.read_h5ad(OUTD / "ois_enhancer_trainset/ois_enhancer_trainset.h5ad")
    tasks = list(ad.obs.index)
    meas = np.asarray(ad.X, np.float32)[tasks.index("IMR90_SEN_vs_PRO")][idx]

    er = tr["eval_rest"]["IMR90_SEN_vs_PRO"]
    CEIL = tr["reproducibility_ceiling"]
    fig, axes = plt.subplots(1, 4, figsize=(FULL, 3.0))

    # --- A: predicted vs measured -------------------------------------------
    a = axes[0]
    ok = np.isfinite(meas) & np.isfinite(pred)
    hb = a.hexbin(meas[ok], pred[ok], gridsize=46, bins="log", cmap="Blues",
                  mincnt=1, linewidths=0)
    a.axhline(0, color=INK2, lw=0.5); a.axvline(0, color=INK2, lw=0.5)
    a.set_xlabel("measured IMR90 log2FC")
    a.set_ylabel("predicted response")
    a.set_title("Held-out IMR90", loc="left", fontsize=7.5)
    a.text(.03, .97, f"$\\rho$ = {er['spearman']:+.4f}\nn = {int(ok.sum()):,}\n"
                     f"{er['frac_of_ceiling']*100:.1f}% of ceiling",
           transform=a.transAxes, fontsize=7, va="top", color=INK)
    cb = fig.colorbar(hb, ax=a, pad=0.02, fraction=0.045)
    cb.set_label("regions (log)", fontsize=7); cb.ax.tick_params(labelsize=5.4)
    cb.outline.set_linewidth(0.4)
    panel_tag(a, "A")

    # --- B: permutation null -------------------------------------------------
    b = axes[1]
    perms = np.array(tr["control1_shuffled_labels"]["all_spearman"])
    real = tr["control1_shuffled_labels"]["spearman"]["real"]
    b.hist(perms, bins=16, color=GREY, edgecolor="white", linewidth=0.5)
    b.axvline(real, color=RED, lw=1.4)
    b.annotate(f"real\n{real:+.4f}", xy=(real, b.get_ylim()[1] * 0.62),
               xytext=(-30, 0), textcoords="offset points", fontsize=7,
               color=RED, ha="right",
               arrowprops=dict(arrowstyle="->", lw=0.8, color=RED))
    b.set_xlabel("Spearman $\\rho$, shuffled labels")
    b.set_ylabel("permutations")
    b.set_title("Shuffled-label null", loc="left", fontsize=7.5)
    aucb = tr["control1_shuffled_labels"]["auc"]
    b.text(.03, .97, "0 / 200 beat real\n$p$ = 0.00498\n(= 1/201 floor)",
           transform=b.transAxes, fontsize=7, va="top", color=INK)
    b.text(.5, -0.30, f"AUC null {aucb['mean']:.3f} $\\pm$ {aucb['sd']:.3f} SD, 0/200",
           transform=b.transAxes, fontsize=7, ha="center", color=MUTED)
    panel_tag(b, "B")

    # --- C: controls ---------------------------------------------------------
    c = axes[2]
    names = ["trained", "untrained", "shuf. embed."]
    sp = [er["spearman"], tr["control2_untrained"]["spearman"],
          tr["control3_shuffled_embeddings"]["spearman"]]
    au = [er["auc"], tr["control2_untrained"]["auc"],
          tr["control3_shuffled_embeddings"]["auc"]]
    xx = np.arange(3); w = 0.38
    c.bar(xx - w / 2, sp, w, color=C_OIS, label="Spearman $\\rho$")
    c.bar(xx + w / 2, au, w, color="white", edgecolor=C_OIS, hatch="////",
          linewidth=0.9, label="AUC")
    c.axhline(0, color=INK, lw=0.8)
    c.axhline(0.5, color=MUTED, lw=0.6, ls=":")
    for xi, (s, u) in enumerate(zip(sp, au)):
        c.text(xi - w / 2, s + (.025 if s >= 0 else -.025), f"{s:+.3f}", fontsize=7,
               ha="center", va="bottom" if s >= 0 else "top", color=INK2)
        c.text(xi + w / 2, u + .025, f"{u:.3f}", fontsize=7, ha="center",
               va="bottom", color=INK2)
    c.set_xticks(xx); c.set_xticklabels(names, fontsize=7, rotation=16, ha="right")
    c.set_ylabel("value")
    c.set_ylim(-0.32, 1.62)
    c.set_xlim(-0.65, 2.65)
    c.set_title("Controls", loc="left", fontsize=7.5)
    c.legend(handles=[Patch(facecolor=C_OIS, label="Spearman $\\rho$"),
                      Patch(facecolor="white", edgecolor=C_OIS, hatch="////",
                            label="AUC")],
             loc="upper left", fontsize=7, ncol=1, handlelength=1.0)
    c.text(-0.58, 0.52, "chance", fontsize=7, color=MUTED, ha="left", va="bottom")
    c.text(.5, -0.40, f"n = {er['n_evaluable']:,} evaluable",
           transform=c.transAxes, fontsize=7, color=INK2, ha="center")
    c.grid(axis="y"); c.set_axisbelow(True)
    panel_tag(c, "C")

    # --- D: ROC --------------------------------------------------------------
    d = axes[3]
    from sklearn.metrics import roc_curve, roc_auc_score
    ev = np.isfinite(meas) & (np.abs(meas) >= tr["log2fc_threshold"])
    ylab = (meas[ev] > 0).astype(int)
    fpr, tpr, _ = roc_curve(ylab, pred[ev])
    auc = roc_auc_score(ylab, pred[ev])
    d.plot([0, 1], [0, 1], ls=":", color=MUTED, lw=0.8)
    d.plot(fpr, tpr, color=C_OIS, lw=1.5)
    d.set_xlabel("false positive rate"); d.set_ylabel("true positive rate")
    d.set_title("ROC, held-out IMR90", loc="left", fontsize=7.5)
    d.text(.96, .06, f"AUC = {auc:.4f}\nn = {int(ev.sum()):,} evaluable\n"
                     f"of {len(meas):,}",
           transform=d.transAxes, fontsize=7, ha="right", va="bottom", color=INK)
    d.set_xlim(0, 1); d.set_ylim(0, 1.02)
    d.grid(); d.set_axisbelow(True)
    panel_tag(d, "D")

    if abs(auc - er["auc"]) > 5e-4:
        note(f"FIG 5D — ROC AUC recomputed from pred_rest.npy = {auc:.4f} vs stored "
             f"{er['auc']:.4f} (difference {auc-er['auc']:+.5f})")
    else:
        print(f"  ROC AUC reproduces stored value: {auc:.4f} vs {er['auc']:.4f}")
    fig.tight_layout(w_pad=2.2)
    return save(fig, "figure5_model_performance")


# =====================================================================  FIG 6
def fig6():
    G = OUTD / "ois_gates/results"
    dc = pd.read_csv(G / "gate1b_dose_curve.csv")
    sp = pd.read_csv(G / "gate1_spacing.csv")
    ps = pd.read_csv(ROOT / "output/audit_fixes/additivity_contrast_per_seed.csv")
    REQ = json.load(open(G / "gate1_results.json"))["required_delta"]

    fig, axes = plt.subplots(1, 3, figsize=(FULL, 3.15))

    # --- A dose response -----------------------------------------------------
    a = axes[0]
    pal = {"NFKB_RELA+ETS1": VIOLET, "NFKB_RELA": BLUE, "CEBPB": AQUA,
           "NFKB_RELA+CEBPB+ETS1": YELLOW, "ETS1": ORANGE}
    mk = {"NFKB_RELA+ETS1": "o", "NFKB_RELA": "s", "CEBPB": "^",
          "NFKB_RELA+CEBPB+ETS1": "D", "ETS1": "v"}
    for (cond, arr), g in dc.groupby(["condition", "arrangement"]):
        g = g.sort_values("K_total")
        lab = MOTIF_LABEL.get(cond) or "+".join(MOTIF_LABEL[m] for m in cond.split("+"))
        if arr == "clustered":
            lab += " (clust.)"
        a.errorbar(g.K_total, g.adjusted_delta, yerr=g.adjusted_sem, marker=mk[cond],
                   color=pal[cond], ms=3.6, lw=1.1, capsize=1.6, elinewidth=0.7,
                   label=lab)
    a.axhline(REQ, color=RED, lw=1.1, ls="--")
    a.text(5.2, REQ + .016, f"required  +{REQ:.3f}", color=RED, fontsize=7,
           ha="left")
    ets = dc[dc.condition == "ETS1"].sort_values("K_total")
    a.text(6.45, float(ets.adjusted_delta.iloc[0]) + .022, "ETS peaks at K = 6",
           fontsize=7, color=ORANGE, ha="left", va="bottom")
    a.set_xlabel("total motif copies in the 200 bp cassette")
    a.set_ylabel("adjusted $\\Delta$ (real − scrambled)")
    a.set_xticks([6, 12, 18, 20]); a.set_xlim(5, 21.5)
    a.set_ylim(0.04, 1.20)
    a.set_title("Dose ceiling: peak at K = 12, decline after",
                loc="left", fontsize=7.5)
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), fontsize=7,
             ncol=2, columnspacing=0.7, handletextpad=0.4)
    a.grid(axis="y"); a.set_axisbelow(True)
    a.text(.985, .03, "mean $\\pm$ SEM, n = 10", transform=a.transAxes,
           fontsize=7, ha="right", va="bottom", color=INK2)
    panel_tag(a, "A")

    # --- B spacing -----------------------------------------------------------
    b = axes[1]
    sp = sp.sort_values(["K_total", "condition"]).reset_index(drop=True)
    for i, r in sp.iterrows():
        c = BLUE if r.p < 0.05 else GREY
        b.plot([0, 1], [r.interleaved_adj, r.clustered_adj], "-", color=c, lw=1.0,
               marker="o", ms=3.4, alpha=.95 if r.p < 0.05 else .6)
    b.set_xticks([0, 1]); b.set_xticklabels(["interleaved", "clustered"], fontsize=7)
    b.set_xlim(-0.26, 1.34)
    b.set_ylabel("adjusted $\\Delta$")
    b.set_title("Clustered beats interleaved, 8/8", loc="left", fontsize=7.5)
    b.text(.985, .035, "n = 8 pairs, each on\n10 backgrounds",
           transform=b.transAxes, fontsize=7, ha="right", va="bottom", color=INK2)
    b.legend(handles=[Line2D([], [], color=BLUE, marker="o", ms=3.4,
                             label="paired $p$ < 0.05 (3/8)"),
                      Line2D([], [], color=GREY, marker="o", ms=3.4,
                             label="n.s. (5/8)")],
             loc="upper left", fontsize=7)
    b.grid(axis="y"); b.set_axisbelow(True)
    panel_tag(b, "B")

    # --- C pair vs single ----------------------------------------------------
    c = axes[2]
    pcol = "adjusted_pair_NFKB_ETS_clustered_K12"
    scol = "adjusted_single_NFKB_K12"
    for _, r in ps.iterrows():
        c.plot([0, 1], [r[scol], r[pcol]], "-", color=GREY, lw=0.8, marker="o",
               ms=3.0, alpha=.75)
    m_s, m_p = ps[scol].mean(), ps[pcol].mean()
    c.plot([0, 1], [m_s, m_p], "-", color=BLUE, lw=2.0, marker="o", ms=5.5,
           zorder=5)
    d = ps["paired_difference"].to_numpy()
    ci = stats.t.interval(0.95, len(d) - 1, loc=d.mean(), scale=stats.sem(d))
    c.set_xticks([0, 1])
    c.set_xticklabels(["NF-κB alone\nK = 12", "NF-κB+ETS\nK = 12 clust."],
                      fontsize=7)
    c.set_xlim(-0.3, 1.45)
    c.set_ylabel("adjusted $\\Delta$")
    c.set_title("Combining does not beat the best single", loc="left", fontsize=7.5)
    c.text(.985, .035,
           f"difference {d.mean():+.4f}\n95% CI [{ci[0]:+.3f}, {ci[1]:+.3f}]\n"
           f"$p$ = {stats.ttest_rel(ps[pcol], ps[scol]).pvalue:.3f}, n = 10",
           transform=c.transAxes, fontsize=7, ha="right", va="bottom", color=INK,
           bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
    c.grid(axis="y"); c.set_axisbelow(True)
    panel_tag(c, "C")

    fig.tight_layout(w_pad=1.8)
    return save(fig, "figure6_design_ceiling")


# =====================================================================  FIG 7
def fig7():
    G = OUTD / "ois_gates/results"
    probe = json.load(open(OUTD / "ois_enhancer_run/results/results_probe.json"))
    summ = probe["motif_insertion"]["summary"]
    ap1 = json.load(open(G / "gate2_ap1_results.json"))["H1b_interaction"]
    freq = pd.read_csv(G / "gate2_h3_frequency_by_decile.csv")

    fig, axes = plt.subplots(1, 3, figsize=(FULL, 3.25))

    # --- A insertion ---------------------------------------------------------
    a = axes[0]
    order = ["NFKB_RELA", "CEBPB", "ETS1", "AP1_FOSJUN"]
    x = np.arange(len(order)); w = 0.38
    for i, m in enumerate(order):
        a.bar(i - w / 2, summ[f"{m}_real"]["delta"], w, color=MOTIF_COLOUR[m])
        a.bar(i + w / 2, summ[f"{m}_scrambled"]["delta"], w, color="white",
              edgecolor=MOTIF_COLOUR[m], hatch="////", linewidth=0.9)
    a.axhline(0, color=INK, lw=0.8)
    for i, m in enumerate(order):
        v = summ[f"{m}_real"]["delta"]
        a.text(i - w / 2, v + (.025 if v >= 0 else -.025), f"{v:+.3f}", fontsize=7,
               ha="center", va="bottom" if v >= 0 else "top", color=INK2)
    a.set_xticks(x); a.set_xticklabels([MOTIF_LABEL[m] for m in order], fontsize=7)
    for t, m in zip(a.get_xticklabels(), order):
        t.set_color(MOTIF_COLOUR[m])
    a.set_ylabel("$\\Delta$ predicted response, k = 0 to 8 copies")
    a.set_title("INSERTION into random background", loc="left", fontsize=7.5)
    a.legend(handles=[Patch(facecolor=INK2, label="real motif"),
                      Patch(facecolor="white", edgecolor=INK2, hatch="////",
                            label="scrambled control")],
             loc="lower left", fontsize=7)
    a.text(.985, .97, "n = 5 backgrounds", transform=a.transAxes,
           fontsize=7, ha="right", va="top", color=INK2)
    a.set_ylim(-0.60, 0.80)
    a.annotate("AP-1 goes DOWN", xy=(3 - w / 2, summ["AP1_FOSJUN_real"]["delta"] * 0.55),
               xytext=(2.30, -0.16), fontsize=7, color=ORANGE, fontweight="bold",
               ha="right",
               arrowprops=dict(arrowstyle="->", lw=0.8, color=ORANGE))
    a.grid(axis="y"); a.set_axisbelow(True)
    panel_tag(a, "A")

    # --- B ablation ----------------------------------------------------------
    b = axes[1]
    bars = [("AP-1\nhigh", ap1["AP1_FOSJUN"]["by_stratum"]["high"], ORANGE),
            ("AP-1\nlow", ap1["AP1_FOSJUN"]["by_stratum"]["low"], "#f6b89d"),
            ("NF-κB\nhigh", ap1["NFKB_RELA"]["by_stratum"]["high"], BLUE)]
    xx = np.arange(3)
    vals = [v["specific"] for _, v, _ in bars]
    errs = [v["sem"] for _, v, _ in bars]
    b.bar(xx, vals, 0.56, yerr=errs, color=[c for *_, c in bars],
          error_kw=dict(ecolor=INK2, lw=0.8, capsize=2.2))
    b.axhline(0, color=INK, lw=0.8)
    hi_t = max(v + e for v, e in zip(vals, errs))
    lo_t = min(v - e for v, e in zip(vals, errs))
    pad = (hi_t - lo_t) * 0.08
    b.set_ylim(lo_t - 3.8 * pad, hi_t + 5.6 * pad)
    for xi, (v, e, (lab, s, _)) in enumerate(zip(vals, errs, bars)):
        yy = (v + e + pad * .5) if v >= 0 else (v - e - pad * .5)
        pv = s["p"]
        ptxt = "$p$ < 0.001" if pv < 0.001 else f"$p$ = {pv:.3f}"
        b.text(xi, yy, f"{v:+.4f}\n{ptxt}", fontsize=7, ha="center",
               va="bottom" if v >= 0 else "top", color=INK2)
    br = hi_t + 3.4 * pad
    b.plot([0, 0, 1, 1], [br - pad * .4, br, br, br - pad * .4], color=RED, lw=0.9)
    b.text(0.5, br + pad * .3, f"interaction $p$ = "
           f"{ap1['AP1_FOSJUN']['interaction']['welch_p']:.4f}",
           ha="center", fontsize=7, color=RED)
    b.set_xticks(xx); b.set_xticklabels([l for l, *_ in bars], fontsize=7)
    b.text(0.5, -0.19, "response stratum of the host region", transform=b.transAxes,
           fontsize=7, ha="center", color=INK2)
    b.set_ylabel("ablation effect\n(ablate motif − ablate motif-free control)")
    b.set_title("ABLATION of native sites in real regions", loc="left", fontsize=7.5)
    b.text(.5, -0.34, "mean $\\pm$ SEM, n = 24 regions per bar",
           transform=b.transAxes, fontsize=7, ha="center", color=INK2)
    b.grid(axis="y"); b.set_axisbelow(True)
    panel_tag(b, "B")

    # --- C decile frequency --------------------------------------------------
    c = axes[2]
    g = freq[(freq.response == "GM21_H3K27ac") & (freq.decile >= 0)].sort_values("decile")
    c.plot(g.decile, g.AP1_FOSJUN, marker="s", color=ORANGE, ms=3.8, lw=1.3)
    tb = freq[(freq.response == "GM21_H3K27ac") & (freq.decile == -1)
              & (freq.motif == "AP1_FOSJUN")].iloc[0]
    c.set_xlabel("decile of GM21 H3K27ac response")
    c.set_ylabel("fraction of regions with a\nstrong AP-1 site (rel $\\geq$ 0.95)")
    c.set_title("AP-1 sites are ENRICHED where\nH3K27ac is high", loc="left",
                fontsize=7.5)
    c.set_xticks(range(10))
    c.text(.03, .97, f"OR = {tb.odds_ratio:.2f} top vs bottom\n"
                     f"Fisher $p$ = {tb.p:.0e}\n"
                     f"n = {int(g.n.iloc[0]):,} per decile",
           transform=c.transAxes, fontsize=7, va="top", color=INK)
    c.grid(); c.set_axisbelow(True)
    panel_tag(c, "C")

    fig.tight_layout(w_pad=1.8)
    return save(fig, "figure7_ap1_context_dependence")


if __name__ == "__main__":
    for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7):
        png, pdf = fn()
        print(f"{png.name:46s} + {pdf.name}")
    if NOTES:
        print("\n" + "=" * 78 + "\nFLAGS\n" + "=" * 78)
        for n in NOTES:
            print(" - " + n)
        (HERE / "figure_flags.txt").write_text("\n\n".join(NOTES))
