#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 92_gate1b_ceiling.py
#
# Gate 1b: dose ceiling at K=18 and K=20; merges Gate 1 into one dose curve. GPU STEP.
#
# CONSUMES: gate1_additivity_summary.csv, via 90_gate_common.py
# PRODUCES: results/gate1b_{ceiling_raw,dose_curve,saturation}.csv, gate1b_results.json, GATE1B_LOG.txt
# ---------------------------------------------------------------------------
"""
GATE 1b -- DOSE CEILING.

Gate 1 found the best condition at K=12 reaches 0.773x of the +0.861 required
delta, and that going from K=6 to K=12 was already sub-linear. This asks whether
pushing the insertion load higher closes the gap or saturates.

K is extended to 18 and 20. 20 is the hard ceiling for the 200 bp cassette: the
longest motif here is 10 bp, so 20 copies tile the cassette end to end with no
room left. If the adjusted delta is still short of +0.861 at K=20, no amount of
additional motif insertion into a 200 bp element closes the gap under this model.

Backgrounds are regenerated with the same RNG seed and draw order as Gate 1, so
the K=6 and K=12 points from that run are directly comparable and are merged in
to form one dose curve per condition.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/workspace/ois_gates/scripts")
from gate_common import (BG, CASSETTE, MOTIF_SEQ, OUT, SEQ_LEN, load_model,
                         load_sizes, make_scorer, pick_batch, scramble,
                         window_on)
from gate1_additivity import (NEUTRAL, N_SEEDS, REQUIRED_DELTA, assignment,
                              build_cassette)

LOADS_B = [18, 20]
ARMS = [("NFKB_RELA", ["NFKB_RELA"], "single"),
        ("CEBPB", ["CEBPB"], "single"),
        ("ETS1", ["ETS1"], "single"),
        ("NFKB_RELA+ETS1", ["NFKB_RELA", "ETS1"], "clustered"),
        ("NFKB_RELA+CEBPB+ETS1", ["NFKB_RELA", "CEBPB", "ETS1"], "clustered")]
LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    t0 = time.time()
    outdir = Path(OUT) / "results"
    head, core = load_model()
    batch = pick_batch(core)
    score = make_scorer(head, core, batch)
    sizes = load_sizes()
    bgseq, ws, off = window_on(NEUTRAL[0], NEUTRAL[1], sizes)
    assert off == SEQ_LEN // 2
    lo = SEQ_LEN // 2 - CASSETTE // 2
    R(f"head + trunk loaded (batch {batch}); context {NEUTRAL[0]}:{NEUTRAL[1]:,}")

    # identical draw order to Gate 1 -> identical backgrounds
    rng = np.random.default_rng(20261004)
    backgrounds = ["".join(rng.choice(list("ACGT"), size=CASSETTE,
                                      p=[BG["A"], BG["C"], BG["G"], BG["T"]]))
                   for _ in range(N_SEEDS)]

    jobs, meta = [], []
    for seed in range(N_SEEDS):
        bg = backgrounds[seed]
        jobs.append(bgseq[:lo] + bg + bgseq[lo + CASSETTE:])
        meta.append(dict(label="empty", K=0, arrangement="none", variant="real",
                         seed=seed, bases_inserted=0))
        for K in LOADS_B:
            for label, members, arr in ARMS:
                for variant in ("real", "scrambled"):
                    assign = assignment(K, members, arr)
                    cas, nb = build_cassette(bg, assign, variant, seed)
                    jobs.append(bgseq[:lo] + cas + bgseq[lo + CASSETTE:])
                    meta.append(dict(label=label, K=K, arrangement=arr,
                                     variant=variant, seed=seed, bases_inserted=nb))
    R(f"{len(jobs)} sequences ({N_SEEDS} seeds x {len(jobs)//N_SEEDS} conditions)")

    preds = []
    CH = 50
    for i in range(0, len(jobs), CH):
        preds.append(score(jobs[i:i + CH]))
        done = min(i + CH, len(jobs))
        el = time.time() - t0
        R(f"  scored {done}/{len(jobs)}  ({el:.0f}s, {el/done:.2f}s/seq)")
    df = pd.DataFrame(meta)
    df["pred"] = np.concatenate(preds)
    base = df[df.label == "empty"].set_index("seed")["pred"]
    df["delta"] = df["pred"] - df["seed"].map(base)
    df.to_csv(outdir / "gate1b_ceiling_raw.csv", index=False)

    rows = []
    for (label, K, arr), g in df[df.label != "empty"].groupby(
            ["label", "K", "arrangement"], sort=False):
        raw = g[g.variant == "real"].sort_values("seed")["delta"].to_numpy()
        sc = g[g.variant == "scrambled"].sort_values("seed")["delta"].to_numpy()
        adj = raw - sc
        t = stats.ttest_1samp(adj, 0.0)
        rows.append(dict(condition=label, K_total=K, arrangement=arr, n_seeds=N_SEEDS,
                         bases_inserted=int(g["bases_inserted"].iloc[0]),
                         raw_delta=float(raw.mean()),
                         raw_sem=float(raw.std(ddof=1) / np.sqrt(N_SEEDS)),
                         scrambled_delta=float(sc.mean()),
                         adjusted_delta=float(adj.mean()),
                         adjusted_sem=float(adj.std(ddof=1) / np.sqrt(N_SEEDS)),
                         adjusted_p=float(t.pvalue),
                         frac_of_required=float(adj.mean() / REQUIRED_DELTA)))
    new = pd.DataFrame(rows)

    # merge Gate 1's K=6/12 points for the same conditions into one dose curve
    g1 = pd.read_csv(outdir / "gate1_additivity_summary.csv")
    keep = [(c, a) for c, _, a in ARMS]
    g1 = g1[[tuple(x) in keep for x in zip(g1.condition, g1.arrangement)]]
    cols = ["condition", "K_total", "arrangement", "n_seeds", "bases_inserted",
            "raw_delta", "raw_sem", "scrambled_delta", "adjusted_delta",
            "adjusted_sem", "adjusted_p", "frac_of_required"]
    curve = pd.concat([g1[cols], new[cols]]).sort_values(
        ["condition", "arrangement", "K_total"])
    curve.to_csv(outdir / "gate1b_dose_curve.csv", index=False)

    R("\n" + "=" * 100)
    R("DOSE CEILING -- adjusted delta vs total inserted copy number")
    R("=" * 100)
    R(f"{'condition':24s} {'arrange':11s} {'K':>3} {'bp':>4} {'adjusted':>10} "
      f"{'+-sem':>7} {'frac of 0.861':>14}")
    for _, r in curve.iterrows():
        R(f"{r.condition:24s} {r.arrangement:11s} {r.K_total:3.0f} {r.bases_inserted:4.0f} "
          f"{r.adjusted_delta:+10.4f} {r.adjusted_sem:7.4f} {r.frac_of_required:14.3f}")

    R("\n" + "=" * 100)
    R("SATURATION -- marginal gain per added copy, between consecutive loads")
    R("=" * 100)
    sat = []
    for (c, a), g in curve.groupby(["condition", "arrangement"]):
        g = g.sort_values("K_total")
        ks = g.K_total.to_numpy(); ds = g.adjusted_delta.to_numpy()
        for i in range(1, len(ks)):
            per = (ds[i] - ds[i - 1]) / (ks[i] - ks[i - 1])
            sat.append(dict(condition=c, arrangement=a, K_from=int(ks[i - 1]),
                            K_to=int(ks[i]), delta_from=float(ds[i - 1]),
                            delta_to=float(ds[i]), gain_per_copy=float(per)))
            R(f"  {c:24s} {a:11s} K {ks[i-1]:2.0f}->{ks[i]:2.0f}  "
              f"{ds[i-1]:+.4f} -> {ds[i]:+.4f}   gain/copy {per:+.5f}")
    pd.DataFrame(sat).to_csv(outdir / "gate1b_saturation.csv", index=False)

    best = curve.loc[curve.adjusted_delta.idxmax()]
    R("\n" + "=" * 100)
    R(f"  best across all loads : {best.condition} ({best.arrangement}, K={best.K_total:.0f})")
    R(f"  adjusted delta        : {best.adjusted_delta:+.4f} +- {best.adjusted_sem:.4f} (n={N_SEEDS})")
    R(f"  required              : {REQUIRED_DELTA:+.4f}")
    R(f"  fraction of required  : {best.adjusted_delta/REQUIRED_DELTA:.3f}x")
    R(f"  shortfall             : {REQUIRED_DELTA - best.adjusted_delta:+.4f}")
    R("=" * 100)

    json.dump(dict(loads_b=LOADS_B, n_seeds=N_SEEDS, required_delta=REQUIRED_DELTA,
                   seconds=time.time() - t0, curve=curve.to_dict("records"),
                   saturation=sat,
                   best=dict(condition=best.condition, arrangement=best.arrangement,
                             K_total=int(best.K_total),
                             adjusted_delta=float(best.adjusted_delta),
                             adjusted_sem=float(best.adjusted_sem),
                             frac_of_required=float(best.adjusted_delta / REQUIRED_DELTA))),
              open(outdir / "gate1b_results.json", "w"), indent=2)
    Path(outdir, "GATE1B_LOG.txt").write_text("\n".join(LOG))
    R(f"\nGATE1B_DONE in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
