#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 98b_gate_figures.py
#
# Working figures for the design-gate report (superseded by 99_figures_main.py).
#
# CONSUMES: results/gate1*.csv, gate2*.csv/json
# PRODUCES: ois_gates/figures/fig[1-5]*.png
# ---------------------------------------------------------------------------
"""Figures for the two gates. Palette: dataviz default categorical slots 1-4,
validated colourblind-safe (worst adjacent pair dE 9.1 protan / 22.9 normal).
Low-contrast slots carry direct labels, as the validator's relief rule requires."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/ois_gates/scripts")
from gate_common import OUT

RES = Path(OUT) / "results"
FIG = Path(OUT) / "figures"
FIG.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SURF = "#fcfcfb"
INK, INK2, MUTED = "#1a1a19", "#4a4a47", "#8a8a85"
REQ = 0.8612759113311768

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.8,
    "font.size": 9, "axes.titlesize": 10, "axes.spines.top": False,
    "axes.spines.right": False, "grid.color": "#e6e6e3", "grid.linewidth": 0.7,
})

MCOL = {"NFKB_RELA": BLUE, "AP1_FOSJUN": ORANGE, "CEBPB": AQUA, "ETS1": YELLOW}
NICE = {"NFKB_RELA": "NF-κB", "AP1_FOSJUN": "AP-1", "CEBPB": "C/EBPβ",
        "ETS1": "ETS"}


def short(c):
    return (c.replace("NFKB_RELA", "NF-κB").replace("CEBPB", "C/EBPβ")
             .replace("ETS1", "ETS"))


# ---------------------------------------------------------------- Fig 1: Gate 1
s = pd.read_csv(RES / "gate1_additivity_summary.csv")
d = s[s.K_total == 12].sort_values("adjusted_delta")
fig, ax = plt.subplots(figsize=(7.4, 4.6))
lab = [f"{short(r.condition)}" + ("" if r.arrangement == "single"
                                  else f"  ({r.arrangement[:5]}.)")
       for _, r in d.iterrows()]
cols = [BLUE if r.n_members == 1 else (ORANGE if r.n_members == 2 else AQUA)
        for _, r in d.iterrows()]
y = np.arange(len(d))
ax.barh(y, d.adjusted_delta, xerr=d.adjusted_sem, color=cols, height=0.68,
        error_kw=dict(ecolor=INK2, lw=1.0, capsize=2.5))
ax.axvline(REQ, color="#d03b3b", lw=1.6, ls="--")
ax.text(REQ, len(d) - 0.2, f"  required  +{REQ:.3f}", color="#d03b3b", fontsize=8.5,
        va="top")
ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=8.5)
ax.set_xlabel("adjusted Δ (real − scrambled), paired within background")
ax.set_title("Gate 1  ·  motif combinations at matched total load (K = 12 copies)\n"
             "no condition reaches the delta needed to move median → 95th percentile",
             loc="left")
ax.set_xlim(0, max(REQ, d.adjusted_delta.max()) * 1.18)
ax.grid(axis="x", alpha=.9); ax.set_axisbelow(True)
for h, (_, r) in zip(y, d.iterrows()):
    ax.text(r.adjusted_delta + r.adjusted_sem + .012, h,
            f"{r.adjusted_delta:+.3f}  ({r.frac_of_required:.2f}×)",
            va="center", fontsize=7.8, color=INK2)
hs = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (BLUE, ORANGE, AQUA)]
ax.legend(hs, ["single motif", "pair", "all three"], frameon=False, fontsize=8.5,
          loc="lower right")
fig.tight_layout(); fig.savefig(FIG / "fig1_gate1_additivity.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------ Fig 2: dose curve
c = pd.read_csv(RES / "gate1b_dose_curve.csv")
fig, ax = plt.subplots(figsize=(7.0, 4.6))
palette = [BLUE, ORANGE, AQUA, YELLOW, "#4a3aa7"]
for i, ((cond, arr), g) in enumerate(c.groupby(["condition", "arrangement"])):
    g = g.sort_values("K_total")
    col = palette[i % len(palette)]
    ax.errorbar(g.K_total, g.adjusted_delta, yerr=g.adjusted_sem, marker="o",
                ms=5.5, lw=2, color=col, capsize=2.5,
                label=short(cond) + ("" if arr == "single" else f" ({arr[:5]}.)"))
    last = g.iloc[-1]
    ax.annotate(short(cond), (last.K_total, last.adjusted_delta),
                textcoords="offset points", xytext=(7, -1), fontsize=7.8, color=col)
ax.axhline(REQ, color="#d03b3b", lw=1.6, ls="--")
ax.text(6.1, REQ + .012, f"required  +{REQ:.3f}", color="#d03b3b", fontsize=8.5)
ax.set_xlabel("total inserted motif copies in the 200 bp cassette")
ax.set_ylabel("adjusted Δ (real − scrambled)")
ax.set_title("Gate 1b  ·  the lever saturates and then reverses\n"
             "every condition peaks at K = 12 and declines as the cassette fills",
             loc="left")
ax.set_xticks([6, 12, 18, 20]); ax.grid(axis="y", alpha=.9); ax.set_axisbelow(True)
ax.set_xlim(5, 23); ax.legend(frameon=False, fontsize=8, loc="lower left", ncol=2)
fig.tight_layout(); fig.savefig(FIG / "fig2_gate1b_dose_ceiling.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------- Fig 3: H1b
ab = pd.read_csv(RES / "gate2_h1b_ablation.csv")
fig, ax = plt.subplots(figsize=(7.0, 4.3))
groups, vals, errs, cols = [], [], [], []
for motif, st in (("AP1_FOSJUN", "high"), ("AP1_FOSJUN", "low"),
                  ("NFKB_RELA", "high")):
    m = ab[(ab.motif == motif) & (ab.stratum == st)]
    piv = m.pivot_table(index="region_index", columns="arm", values="effect")
    spec = (piv["ablated"] - piv["position_control"]).dropna()
    groups.append(f"{NICE[motif]}\n{st}-response\nregions (n={len(spec)})")
    vals.append(spec.mean()); errs.append(spec.std(ddof=1) / np.sqrt(len(spec)))
    cols.append(MCOL[motif] if st == "high" else "#f2a989")
x = np.arange(len(groups))
ax.bar(x, vals, yerr=errs, color=cols, width=.6,
       error_kw=dict(ecolor=INK2, lw=1.1, capsize=3.5))
ax.axhline(0, color=INK, lw=1.0)
ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8.5)
ax.set_ylabel("motif-specific ablation effect\n(ablate motif − ablate motif-free control)")
ax.set_title("Gate 2 / H1  ·  ablating the native motif: the sign depends on context\n"
             "negative = the motif was contributing POSITIVELY in that real region",
             loc="left")
hi_t = max(v + e for v, e in zip(vals, errs))
lo_t = min(v - e for v, e in zip(vals, errs))
pad = (hi_t - lo_t) * 0.07
ax.set_ylim(lo_t - 3.2 * pad, hi_t + 4.4 * pad)
for xi, (v, e) in enumerate(zip(vals, errs)):
    yy = (v + e + pad * .45) if v >= 0 else (v - e - pad * .45)
    ax.text(xi, yy, f"{v:+.4f}", ha="center",
            va="bottom" if v >= 0 else "top", fontsize=8.5, color=INK2)
br = hi_t + 2.2 * pad
ax.plot([0, 0, 1, 1], [br - pad * .35, br, br, br - pad * .35], color="#d03b3b", lw=1.1)
ax.text(.5, br + pad * .3, "AP-1 high vs low:  interaction p = 0.009", ha="center",
        fontsize=8.5, color="#d03b3b")
ax.grid(axis="y", alpha=.9); ax.set_axisbelow(True)
fig.tight_layout(); fig.savefig(FIG / "fig3_h1b_ablation.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------- Fig 4: H1c + H2
h1c = pd.read_csv(RES / "gate2_h1c_contexts.csv")
h2 = pd.read_csv(RES / "gate2_h2_variants.csv")
fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4))
a = axes[0]
g = h1c[h1c.k == 8].sort_values("adjusted")
a.barh(np.arange(len(g)), g.adjusted, xerr=g["sem"], color=ORANGE, height=.6,
       error_kw=dict(ecolor=INK2, lw=1.0, capsize=2.5))
a.axvline(0, color=INK, lw=1.0)
a.set_yticks(np.arange(len(g))); a.set_yticklabels(g.context, fontsize=8.5)
a.set_xlabel("adjusted Δ at 8 copies")
a.set_title("H1c  ·  AP-1 inserted into four neutral contexts\n"
            "negative in three of four — not one odd background", loc="left")
for i, (_, r) in enumerate(g.iterrows()):
    off = r["sem"] + 0.045
    a.text(r.adjusted + (off if r.adjusted >= 0 else -off), i, f"{r.adjusted:+.3f}",
           va="center", ha="left" if r.adjusted >= 0 else "right", fontsize=8,
           color=INK2)
a.grid(axis="x", alpha=.9); a.set_axisbelow(True)
a.set_xlim(g.adjusted.min() - g["sem"].max() - .26,
           max(g.adjusted.max() + g["sem"].max() + .20, .24))

b = axes[1]
canon = h2["seq"].isin(["TGACTCA", "TGAGTCA"])
b.errorbar(h2.ap1_best_rel[~canon], h2.adjusted_delta[~canon],
           yerr=h2["sem"][~canon], fmt="o", ms=6, color=BLUE, lw=0, elinewidth=1,
           ecolor="#9ec5f4", label="other variants")
b.errorbar(h2.ap1_best_rel[canon], h2.adjusted_delta[canon], yerr=h2["sem"][canon],
           fmt="o", ms=9, color=ORANGE, lw=0, elinewidth=1.2, ecolor="#f2a989",
           label="canonical TRE  TGA(C/G)TCA")
for _, r in h2[canon | (h2["variant"] == "CRE8")].iterrows():
    b.annotate(r["seq"], (r.ap1_best_rel, r.adjusted_delta), fontsize=8,
               textcoords="offset points", xytext=(8, -3), color=INK2)
b.axhline(0, color=INK, lw=1.0)
b.set_xlabel("AP-1 PWM match of the built cassette (relative score)")
b.set_ylabel("adjusted Δ at 8 copies")
b.set_title("H2  ·  hyper-specific to the exact TRE heptamer\n"
            "it does not grade with PWM affinity (mutants ρ = +0.13)",
            loc="left")
b.legend(frameon=False, fontsize=8, loc="lower left")
b.grid(alpha=.9); b.set_axisbelow(True)
fig.tight_layout(); fig.savefig(FIG / "fig4_h1c_h2.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------- Fig 5: H3
freq = pd.read_csv(RES / "gate2_h3_frequency_by_decile.csv")
tab = pd.read_csv(RES / "gate2_h3_motif_vs_assay.csv")
fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4))
a = axes[0]
for resp, col, mk in (("GM21_H3K27ac", ORANGE, "o"), ("GM21_ATAC", BLUE, "s")):
    g = freq[(freq.response == resp) & (freq.decile >= 0)].sort_values("decile")
    a.plot(g.decile, g.AP1_FOSJUN, marker=mk, ms=5.5, lw=2, color=col,
           label="H3K27ac response" if "K27" in resp else "ATAC response")
    a.annotate("H3K27ac" if "K27" in resp else "ATAC",
               (g.decile.iloc[-1], g.AP1_FOSJUN.iloc[-1]), fontsize=8.5, color=col,
               textcoords="offset points", xytext=(6, -2))
a.set_xlabel("decile of measured GM21 response (0 = lowest)")
a.set_ylabel("fraction of regions with a strong AP-1 site")
a.set_title("H3  ·  AP-1 sites are ENRICHED where senescence H3K27ac is high\n"
            "and they track acetylation more steeply than accessibility", loc="left")
a.legend(frameon=False, fontsize=8.5, loc="lower right")
a.grid(alpha=.9); a.set_axisbelow(True); a.set_xlim(-.4, 10.2)

b = axes[1]
t = tab.sort_values("z_gap_difference")
x = np.arange(len(t))
b.bar(x, t.z_gap_difference, color=[MCOL[m] for m in t.motif], width=.6)
b.axhline(0, color=INK, lw=1.0)
b.set_xticks(x); b.set_xticklabels([NICE[m] for m in t.motif], fontsize=9)
b.set_ylabel("z(ATAC) − z(H3K27ac),  motif-present − absent")
b.set_title("H3 predicts AP-1 shifted toward ACCESSIBILITY (positive)\n"
            "measured: AP-1 is the most ACETYLATION-shifted of the four", loc="left")
for xi, v in zip(x, t.z_gap_difference):
    b.text(xi, v + (.006 if v >= 0 else -.006), f"{v:+.3f}", ha="center",
           va="bottom" if v >= 0 else "top", fontsize=8.5, color=INK2)
b.grid(axis="y", alpha=.9); b.set_axisbelow(True)
_sp = t.z_gap_difference
b.set_ylim(_sp.min() - abs(_sp.min()) * .24, _sp.max() + abs(_sp.max()) * .46)
fig.tight_layout(); fig.savefig(FIG / "fig5_h3_atac_vs_k27ac.png", dpi=200)
plt.close(fig)

print("FIGURES_DONE:", sorted(p.name for p in FIG.glob("*.png")))
