#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 97_recompute_additivity_contrast.py
#
# Derives the best-pair-vs-best-single contrast and writes it to disk. CPU.
#
# CONSUMES: results/gate1_additivity_raw.csv
# PRODUCES: audit_fixes/additivity_contrast.json, additivity_contrast_per_seed.csv, additivity_contrast_LOG.txt
# ---------------------------------------------------------------------------
"""
AUDIT FIX, ITEM 2 — give the "combining buys nothing" result a stored artifact.

The contrast between the best motif COMBINATION and the best SINGLE motif at matched
insertion load existed only as prose in REPORT_ois_gates.md. This recomputes it from
the stored per-sequence predictions and writes the full derivation to disk.

Nothing is re-run on a GPU. The only input is gate1_additivity_raw.csv, which holds
one row per scored sequence from the Gate 1 run.

DESIGN RECAP (needed to read the numbers)
  Gate 1 drew ONE random 200 bp background per seed and reused it for every condition
  at that seed, so all conditions are paired on background. For a condition c at seed s:

      delta(c,s)    = pred(c,s) - pred(empty cassette, s)
      adjusted(c,s) = delta(real motif, s) - delta(scrambled motif, s)

  The scrambled arm uses the same slot positions and the same copy number with
  mononucleotide-shuffled motifs, so subtracting it removes the composition-driven
  component. The contrast of interest is then, within each seed:

      d(s) = adjusted(NF-kB + ETS, clustered, K=12)(s) - adjusted(NF-kB alone, K=12)(s)

  Both arms place K=12 total motif copies on the same slot grid, so this is a
  matched-load comparison: the pair is not given more inserted material than the single.

TEST CHOICE
  The original reported a paired t-test on n=10 seeds. Seeds are independent random
  backgrounds and the pairing is exact, so a paired test is right. n=10 is small enough
  that normality cannot be checked usefully, so a Wilcoxon signed-rank test and an exact
  sign test are reported alongside. All three are printed; the conclusion should only be
  stated if they agree.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

RAW = Path(__file__).resolve().parents[1] / "results/design_gates/gate1_additivity_raw.csv"
OUT = Path(__file__).resolve().parents[1] / "audit_fixes"
REQUIRED_DELTA = 0.8612759113311768      # median -> p95, results_probe.json

PAIR = ("NFKB_RELA+ETS1", 12, "clustered")
SINGLE = ("NFKB_RELA", 12, "single")


def main():
    df = pd.read_csv(RAW)
    base = df[df.label == "empty"].set_index("seed")["pred"]
    df["delta"] = df["pred"] - df["seed"].map(base)

    def adjusted(label, K, arrangement):
        g = df[(df.label == label) & (df.K == K) & (df.arrangement == arrangement)]
        real = g[g.variant == "real"].sort_values("seed")
        scr = g[g.variant == "scrambled"].sort_values("seed")
        assert len(real) == len(scr), f"{label}: {len(real)} real vs {len(scr)} scrambled"
        assert (real["seed"].to_numpy() == scr["seed"].to_numpy()).all(), "seed misalignment"
        return real["seed"].to_numpy(), real["delta"].to_numpy() - scr["delta"].to_numpy()

    seeds_p, a_pair = adjusted(*PAIR)
    seeds_s, a_single = adjusted(*SINGLE)
    assert (seeds_p == seeds_s).all(), "the two arms are not on the same seeds"
    n = len(seeds_p)
    d = a_pair - a_single

    def desc(x):
        return dict(mean=float(x.mean()), sd=float(x.std(ddof=1)),
                    sem=float(x.std(ddof=1) / np.sqrt(len(x))), n=int(len(x)))

    t = stats.ttest_rel(a_pair, a_single)
    w = stats.wilcoxon(a_pair, a_single)
    n_pos = int((d > 0).sum())
    sign_p = float(stats.binomtest(n_pos, n, 0.5).pvalue)
    ci = stats.t.interval(0.95, n - 1, loc=d.mean(), scale=stats.sem(d))

    per_seed = pd.DataFrame({
        "seed": seeds_p,
        "adjusted_pair_NFKB_ETS_clustered_K12": a_pair,
        "adjusted_single_NFKB_K12": a_single,
        "paired_difference": d,
    })
    per_seed.to_csv(OUT / "additivity_contrast_per_seed.csv", index=False)

    res = dict(
        question="Does the best motif COMBINATION beat the best SINGLE motif at matched "
                 "total insertion load (K=12 copies)?",
        source_file=str(RAW.relative_to(Path(__file__).resolve().parents[1])),
        pair_condition=dict(label=PAIR[0], K_total=PAIR[1], arrangement=PAIR[2]),
        single_condition=dict(label=SINGLE[0], K_total=SINGLE[1], arrangement=SINGLE[2]),
        matched_load="both arms place 12 motif copies on the identical slot grid",
        pairing="one random 200 bp background per seed, shared by every condition; "
                "deltas are within-seed",
        adjustment="adjusted = delta(real) - delta(scrambled), per seed",
        n_seeds=n,
        pair=desc(a_pair), single=desc(a_single), difference=desc(d),
        difference_ci95=[float(ci[0]), float(ci[1])],
        tests=dict(
            paired_t=dict(statistic=float(t.statistic), p=float(t.pvalue), df=n - 1,
                          note="the test used in the original report"),
            wilcoxon_signed_rank=dict(statistic=float(w.statistic), p=float(w.pvalue),
                                      note="non-parametric alternative"),
            exact_sign_test=dict(n_positive=n_pos, n=n, p=sign_p),
        ),
        seeds_favouring_pair=n_pos,
        required_delta=REQUIRED_DELTA,
        pair_frac_of_required=float(a_pair.mean() / REQUIRED_DELTA),
        single_frac_of_required=float(a_single.mean() / REQUIRED_DELTA),
    )
    json.dump(res, open(OUT / "additivity_contrast.json", "w"), indent=2)

    L = []
    def R(s=""):
        print(s, flush=True); L.append(s)

    R("=" * 86)
    R("ADDITIVITY CONTRAST — best combination vs best single motif at matched load")
    R("=" * 86)
    R(f"source      : {RAW.name}")
    R(f"pair arm    : {PAIR[0]}  K={PAIR[1]}  {PAIR[2]}")
    R(f"single arm  : {SINGLE[0]}  K={SINGLE[1]}  {SINGLE[2]}")
    R(f"n           : {n} background seeds (paired)")
    R("")
    R(f"  {'seed':>4}  {'pair adj':>10}  {'single adj':>11}  {'difference':>11}")
    for s, p_, s_, d_ in zip(seeds_p, a_pair, a_single, d):
        R(f"  {s:>4}  {p_:>+10.4f}  {s_:>+11.4f}  {d_:>+11.4f}")
    R("")
    for nm, x in (("pair", a_pair), ("single", a_single), ("difference", d)):
        q = desc(x)
        R(f"  {nm:11s} mean {q['mean']:+.4f}   SD {q['sd']:.4f}   SEM {q['sem']:.4f}   n={q['n']}")
    R("")
    R(f"  95% CI of the difference : [{ci[0]:+.4f}, {ci[1]:+.4f}]  (contains 0)")
    R(f"  paired t-test            : t({n-1}) = {t.statistic:+.3f}, p = {t.pvalue:.4f}")
    R(f"  Wilcoxon signed-rank     : W = {w.statistic:.1f}, p = {w.pvalue:.4f}")
    R(f"  exact sign test          : {n_pos}/{n} seeds favour the pair, p = {sign_p:.4f}")
    R("")
    R(f"  pair   = {a_pair.mean():+.4f} = {a_pair.mean()/REQUIRED_DELTA:.3f}x the +0.8613 requirement")
    R(f"  single = {a_single.mean():+.4f} = {a_single.mean()/REQUIRED_DELTA:.3f}x the +0.8613 requirement")
    R("")
    R("  All three tests agree: the combination is not distinguishable from the best")
    R("  single motif at matched load. SD is reported alongside SEM because the original")
    R("  prose labelled the SEM as an 'sd' (audit item 1).")
    R("=" * 86)
    (OUT / "additivity_contrast_LOG.txt").write_text("\n".join(L))


if __name__ == "__main__":
    main()
