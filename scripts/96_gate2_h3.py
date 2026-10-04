#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 96_gate2_h3.py
#
# Gate 2 H3: matched ATAC-vs-H3K27ac contrast and motif frequency by response decile. CPU.
#
# CONSUMES: results/gate2_motif_calls.csv
# PRODUCES: results/gate2_h3_{motif_vs_assay,dose,frequency_by_decile,conditional}.csv, gate2_h3_results.json
# ---------------------------------------------------------------------------
"""
GATE 2, H3 -- ASSAY-SPECIFIC BIOLOGY (no GPU needed).

Claim under test: AP-1 is a pioneer factor that establishes ACCESSIBILITY. This
model is trained on H3K27ac, which marks the nucleosomes FLANKING an element and
is depleted over the nucleosome-free region itself. If AP-1 occupancy displaces
nucleosomes, regions with strong AP-1 could show an elevated ATAC response while
their H3K27ac response stays flat or goes negative.

The test is a matched within-study contrast. GM21 (GSE205898) has both H3K27ac
(RAS vs EV, n=2 vs 2) and ATAC (RAS vs EV, n=4 vs 4) over the SAME regions and the
SAME perturbation, so the two readouts differ only in assay.

The control that decides whether anything here is AP-1 SPECIFIC: the identical
contrast is run for NF-kB, C/EBP-beta and ETS. If every motif shows ATAC above
H3K27ac then what is being measured is a general property of the two assays, not
something about AP-1. That control is the difference between a real finding and
an artifact of comparing two assays with different dynamic ranges, so the
ATAC-minus-K27ac gap is reported for all four motifs and AP-1 is only called
special if its gap stands apart.

Both readouts are z-scored across regions before differencing, because an
untransformed ATAC log2FC and an untransformed H3K27ac log2FC are not on a common
scale and their raw difference would be meaningless.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/workspace/ois_gates/scripts")
from gate_common import OUT

MOTIFS = ["AP1_FOSJUN", "NFKB_RELA", "CEBPB", "ETS1"]
THRESH = 95          # use the count_95 columns: the most stringent calls
LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def z(x):
    x = np.asarray(x, float)
    m = np.isfinite(x)
    out = np.full_like(x, np.nan)
    out[m] = (x[m] - x[m].mean()) / x[m].std(ddof=1)
    return out


def main():
    outdir = Path(OUT) / "results"
    d = pd.read_csv(outdir / "gate2_motif_calls.csv")
    R(f"regions: {len(d)}")

    K27 = "GM21_RAS_vs_EV"
    ATAC = "GM21_ATAC_RAS_vs_EV"
    ok = np.isfinite(d[K27]) & np.isfinite(d[ATAC])
    d = d[ok].copy()
    R(f"regions with both GM21 H3K27ac and matched GM21 ATAC: {len(d)}")
    R(f"  H3K27ac n=2 vs 2 replicates; ATAC n=4 vs 4 replicates (different n per assay)")
    R(f"  corr(H3K27ac, ATAC) across all regions: "
      f"Pearson {stats.pearsonr(d[K27], d[ATAC]).statistic:+.4f}, "
      f"Spearman {stats.spearmanr(d[K27], d[ATAC]).statistic:+.4f}")
    d["z_k27"] = z(d[K27]); d["z_atac"] = z(d[ATAC])
    d["z_gap"] = d["z_atac"] - d["z_k27"]

    res = {"n_regions": int(len(d)),
           "k27_vs_atac_pearson": float(stats.pearsonr(d[K27], d[ATAC]).statistic),
           "k27_vs_atac_spearman": float(stats.spearmanr(d[K27], d[ATAC]).statistic)}

    # ---------------------------------------------------------------- main table
    R("\n" + "=" * 112)
    R(f"H3 -- PRESENCE OF EACH MOTIF (count_{THRESH} >= 1) vs THE TWO ASSAYS, same cells, "
      f"same contrast")
    R("=" * 112)
    R(f"  {'motif':12s} {'n_with':>7s} {'n_without':>9s} "
      f"{'K27 with':>9s} {'K27 w/o':>9s} {'dK27':>8s} {'p':>9s} | "
      f"{'ATAC with':>9s} {'ATAC w/o':>9s} {'dATAC':>8s} {'p':>9s} | {'z-gap':>7s}")
    rows = []
    for m in MOTIFS:
        has = d[f"{m}_count_{THRESH}"].to_numpy() >= 1
        a, b = d[has], d[~has]
        dk = a[K27].mean() - b[K27].mean()
        da = a[ATAC].mean() - b[ATAC].mean()
        pk = stats.mannwhitneyu(a[K27], b[K27]).pvalue
        pa = stats.mannwhitneyu(a[ATAC], b[ATAC]).pvalue
        gap = a["z_gap"].mean() - b["z_gap"].mean()
        pg = stats.mannwhitneyu(a["z_gap"], b["z_gap"]).pvalue
        R(f"  {m:12s} {has.sum():7d} {(~has).sum():9d} "
          f"{a[K27].mean():+9.4f} {b[K27].mean():+9.4f} {dk:+8.4f} {pk:9.2e} | "
          f"{a[ATAC].mean():+9.4f} {b[ATAC].mean():+9.4f} {da:+8.4f} {pa:9.2e} | "
          f"{gap:+7.4f}")
        rows.append(dict(motif=m, n_with=int(has.sum()), n_without=int((~has).sum()),
                         k27_with=float(a[K27].mean()), k27_without=float(b[K27].mean()),
                         delta_k27=float(dk), p_k27=float(pk),
                         atac_with=float(a[ATAC].mean()), atac_without=float(b[ATAC].mean()),
                         delta_atac=float(da), p_atac=float(pa),
                         z_gap_difference=float(gap), p_z_gap=float(pg)))
    tab = pd.DataFrame(rows)
    tab.to_csv(outdir / "gate2_h3_motif_vs_assay.csv", index=False)
    res["presence"] = rows
    R("\n  z-gap = (z ATAC) - (z H3K27ac); a POSITIVE difference means the motif's")
    R("  regions are shifted toward accessibility relative to acetylation.")
    R("  H3 requires AP-1's z-gap to be positive AND clearly larger than the other three.")

    # ---------------------------------------------------------------- dose
    R("\n" + "=" * 112)
    R("DOSE -- motif copy number vs each assay (Spearman over all regions)")
    R("=" * 112)
    R(f"  {'motif':12s} {'rho vs K27':>11s} {'p':>10s} {'rho vs ATAC':>12s} {'p':>10s} "
      f"{'rho vs z-gap':>13s} {'p':>10s}")
    dose = []
    for m in MOTIFS:
        c = d[f"{m}_count_{THRESH}"].to_numpy()
        r1 = stats.spearmanr(c, d[K27]); r2 = stats.spearmanr(c, d[ATAC])
        r3 = stats.spearmanr(c, d["z_gap"])
        R(f"  {m:12s} {r1.statistic:+11.4f} {r1.pvalue:10.2e} {r2.statistic:+12.4f} "
          f"{r2.pvalue:10.2e} {r3.statistic:+13.4f} {r3.pvalue:10.2e}")
        dose.append(dict(motif=m, rho_k27=float(r1.statistic), p_k27=float(r1.pvalue),
                         rho_atac=float(r2.statistic), p_atac=float(r2.pvalue),
                         rho_zgap=float(r3.statistic), p_zgap=float(r3.pvalue)))
    pd.DataFrame(dose).to_csv(outdir / "gate2_h3_dose.csv", index=False)
    res["dose"] = dose

    # ------------------------------------------- motif frequency by response decile
    R("\n" + "=" * 112)
    R("MOTIF FREQUENCY IN HIGH- vs LOW-RESPONSE REAL REGIONS")
    R("  (the direct check on whether AP-1 sites sit where the model predicts high activity)")
    R("=" * 112)
    freq_out = []
    for resp_name, resp in (("GM21_H3K27ac", d[K27]), ("GM21_ATAC", d[ATAC]),
                            ("IMR90_H3K27ac", d["IMR90_SEN_vs_PRO"])):
        rv = np.asarray(resp, float)
        m = np.isfinite(rv)
        dec = pd.qcut(pd.Series(rv[m]), 10, labels=False, duplicates="drop")
        sub = d[m]
        R(f"\n  {resp_name}: fraction of regions with >=1 call (count_{THRESH}), by decile")
        R(f"  {'decile':>7s} {'n':>7s} {'mean resp':>10s} " +
          " ".join(f"{x:>12s}" for x in MOTIFS))
        for q in range(int(dec.max()) + 1):
            sel = (dec == q).to_numpy()
            fr = [float((sub[f"{x}_count_{THRESH}"].to_numpy()[sel] >= 1).mean())
                  for x in MOTIFS]
            R(f"  {q:7d} {sel.sum():7d} {rv[m][sel].mean():+10.4f} " +
              " ".join(f"{v:12.4f}" for v in fr))
            freq_out.append(dict(response=resp_name, decile=int(q), n=int(sel.sum()),
                                 mean_response=float(rv[m][sel].mean()),
                                 **{x: fr[i] for i, x in enumerate(MOTIFS)}))
        # top vs bottom decile, explicit
        hi = (dec == dec.max()).to_numpy(); lo = (dec == 0).to_numpy()
        R(f"  top vs bottom decile (n={hi.sum()} vs {lo.sum()}):")
        for x in MOTIFS:
            cx = sub[f"{x}_count_{THRESH}"].to_numpy()
            fh, fl = float((cx[hi] >= 1).mean()), float((cx[lo] >= 1).mean())
            od = stats.fisher_exact([[int((cx[hi] >= 1).sum()), int((cx[hi] < 1).sum())],
                                     [int((cx[lo] >= 1).sum()), int((cx[lo] < 1).sum())]])
            R(f"    {x:12s} top {fh:.4f}  bottom {fl:.4f}  ratio {fh/fl if fl else np.nan:6.3f}  "
              f"OR {od.statistic:6.3f}  p {od.pvalue:.3e}")
            freq_out.append(dict(response=resp_name, decile=-1, n=int(hi.sum()),
                                 mean_response=np.nan, motif=x, top_frac=fh,
                                 bottom_frac=fl, odds_ratio=float(od.statistic),
                                 p=float(od.pvalue)))
    pd.DataFrame(freq_out).to_csv(outdir / "gate2_h3_frequency_by_decile.csv", index=False)

    # ---------------------------------------- conditional: within accessible regions
    R("\n" + "=" * 112)
    R("CONDITIONAL -- among regions whose ATAC response is in the TOP quartile,")
    R("  does AP-1 presence still predict lower H3K27ac? (the sharpest form of H3)")
    R("=" * 112)
    qa = d[ATAC].quantile(0.75)
    sub = d[d[ATAC] >= qa]
    R(f"  n = {len(sub)} regions with GM21 ATAC response >= {qa:+.4f} (top quartile)")
    cond = []
    for m in MOTIFS:
        has = sub[f"{m}_count_{THRESH}"].to_numpy() >= 1
        if has.sum() < 20 or (~has).sum() < 20:
            continue
        dk = sub[K27][has].mean() - sub[K27][~has].mean()
        p = stats.mannwhitneyu(sub[K27][has], sub[K27][~has]).pvalue
        R(f"  {m:12s} n_with {has.sum():6d}  K27 with {sub[K27][has].mean():+.4f}  "
          f"without {sub[K27][~has].mean():+.4f}  diff {dk:+.4f}  p {p:.3e}")
        cond.append(dict(motif=m, n_with=int(has.sum()), n_without=int((~has).sum()),
                         delta_k27=float(dk), p=float(p)))
    pd.DataFrame(cond).to_csv(outdir / "gate2_h3_conditional.csv", index=False)
    res["conditional_top_atac_quartile"] = cond

    json.dump(res, open(outdir / "gate2_h3_results.json", "w"), indent=2, default=str)
    Path(outdir, "GATE2_H3_LOG.txt").write_text("\n".join(LOG))
    R("\nGATE2_H3_DONE")


if __name__ == "__main__":
    main()
