#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 98_fig_style.py
#
# Shared publication style, validated palette, size constants. Imported, not run.
#
# CONSUMES: none
# PRODUCES: none (module)
# ---------------------------------------------------------------------------
"""
Shared publication style for the senescence enhancer manuscript figures.

Spec: 300 dpi PNG + vector PDF; 180 mm full width, 85 mm single column;
sans-serif, >= 7 pt at final size; colourblind-safe categorical palette;
no gradients, no 3D, no chartjunk.

Palette is the dataviz default categorical ramp, validated with
scripts/validate_palette.js against a white surface:
  3-slot : all checks pass, worst adjacent CVD dE 9.2 (deutan)
  6-slot : all checks pass, worst adjacent CVD dE 9.1 (protan)
The validator raises a contrast WARN for the lighter slots against white, which
obliges "relief": every series carrying one of those colours also gets a direct
label or a distinct marker, never colour alone.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

MM = 1.0 / 25.4
FULL = 180 * MM          # 7.087 in
SINGLE = 85 * MM         # 3.346 in

OUT = Path("/Users/gabeykim/Downloads/Senescence/output/figures_manuscript")

# categorical slots, fixed order, never cycled
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, VIOLET = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7")
RED = "#d03b3b"          # status: reserved for reference lines / thresholds
GREY = "#8a8a85"
INK, INK2, MUTED = "#1a1a19", "#4a4a47", "#8a8a85"
GRID = "#e6e6e3"
SHADE = "#d9d9d6"

# semantic assignments, used consistently across every figure
C_OIS = BLUE
C_REP = ORANGE
C_CROSS = GREY
C_IR = AQUA
MOTIF_COLOUR = {"NFKB_RELA": BLUE, "AP1_FOSJUN": ORANGE,
                "CEBPB": AQUA, "ETS1": YELLOW}
MOTIF_MARKER = {"NFKB_RELA": "o", "AP1_FOSJUN": "s",
                "CEBPB": "^", "ETS1": "D"}      # secondary encoding, per the relief rule
MOTIF_LABEL = {"NFKB_RELA": "NF-κB", "AP1_FOSJUN": "AP-1",
               "CEBPB": "C/EBPβ", "ETS1": "ETS"}


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "grid.color": GRID, "grid.linewidth": 0.5,
        "legend.frameon": False, "legend.handlelength": 1.4,
        "legend.borderpad": 0.2, "legend.labelspacing": 0.3,
        "lines.linewidth": 1.2, "lines.markersize": 4,
        "pdf.fonttype": 42, "ps.fonttype": 42,   # embed TrueType, keep text editable
        "figure.dpi": 300,
    })


def panel_tag(ax, letter, dx=-0.085, dy=1.045):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left", color=INK)


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    png, pdf = OUT / f"{stem}.png", OUT / f"{stem}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf
