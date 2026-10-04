#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 91_gate1_additivity.py
#
# Gate 1: motif additivity at matched total copy number, plus spacing. GPU STEP.
#
# CONSUMES: via 90_gate_common.py
# PRODUCES: results/gate1_additivity_{raw,summary,vs_expected}.csv, gate1_spacing.csv, gate1_results.json, GATE1_LOG.txt
# ---------------------------------------------------------------------------
"""
GATE 1 -- MOTIF ADDITIVITY.

Question: do NF-kB, C/EBP-beta and ETS combine additively toward the +0.861 delta
required to move a sequence from the median to the 95th percentile of predicted
activity, or do they saturate below it?

LOAD BALANCING (the thing that would otherwise make combinations look good for
trivial reasons):
  * TOTAL copy number is held constant within each load level K in {6, 12}.
    A single-motif arm places K copies of one motif; a pair places K/2 of each;
    the triple places K/3 of each. Every arm at a given K therefore places
    exactly K motif instances.
  * The SLOT GRID is identical across every arm at a given K: slot j occupies
    span = 200 // K bases starting at j*span, and the motif is centred in its
    slot. Only the IDENTITY of the motif in each slot changes between arms.
  * Motif lengths differ (NF-kB 10 bp, C/EBPb 8 bp, ETS 7 bp), so equal copy
    number is NOT equal inserted bases. Mixtures containing the shorter motifs
    insert FEWER bases than an NF-kB-only arm at the same copy number. That is
    conservative: it biases against combinations, not toward them. Inserted
    bases per arm are reported alongside every result.

PAIRED DESIGN: one random 200 bp background is drawn per seed and reused by
every arm at that seed, so each arm's delta is computed against its own seed's
empty-cassette baseline. This removes background-to-background variance from
the comparison (an improvement on the prior single-motif probe, which drew an
independent background per arm).
"""
import itertools
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

MOTIFS3 = ["NFKB_RELA", "CEBPB", "ETS1"]
LOADS = [6, 12]
N_SEEDS = 10
REQUIRED_DELTA = 0.8612759113311768          # median -> p95, from results_probe.json
NEUTRAL = ("chr1", 265888)                   # design_anchor of the neutral context used before

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def assignment(K, members, arrangement):
    """Which motif sits in each of the K slots."""
    n = len(members)
    if arrangement == "single":
        return [members[0]] * K
    if arrangement == "interleaved":
        return [members[j % n] for j in range(K)]
    if arrangement == "clustered":
        per = K // n
        return [members[min(j // per, n - 1)] for j in range(K)]
    raise ValueError(arrangement)


def build_cassette(bg, assign, variant, seed):
    """Place the assigned motifs on the shared slot grid. Returns (seq, bases_inserted)."""
    K = len(assign)
    cas = list(bg)
    if K == 0:
        return "".join(cas), 0
    span = CASSETTE // K
    nb = 0
    for j, m in enumerate(assign):
        use = MOTIF_SEQ[m] if variant == "real" else scramble(MOTIF_SEQ[m],
                                                              seed * 1009 + j * 31 + len(m))
        L = len(use)
        st = max(0, min(CASSETTE - L, j * span + (span - L) // 2))
        cas[st:st + L] = list(use)
        nb += L
    return "".join(cas), nb


def main():
    t0 = time.time()
    outdir = Path(OUT) / "results"
    outdir.mkdir(parents=True, exist_ok=True)

    head, core = load_model()
    batch = pick_batch(core)
    score = make_scorer(head, core, batch)
    R(f"head + trunk loaded; trunk batch size {batch}")

    sizes = load_sizes()
    bgseq, ws, off = window_on(NEUTRAL[0], NEUTRAL[1], sizes)
    assert off == SEQ_LEN // 2, f"anchor not centred (offset {off}) -- window was clipped"
    lo = SEQ_LEN // 2 - CASSETTE // 2
    R(f"neutral context {NEUTRAL[0]}:{NEUTRAL[1]:,}  window_start={ws:,}  cassette at [{lo}, {lo+CASSETTE})")

    # ---- arm table ------------------------------------------------------------
    arms = [dict(label="empty", members=[], K=0, arrangement="none")]
    for K in LOADS:
        for m in MOTIFS3:
            arms.append(dict(label=f"{m}", members=[m], K=K, arrangement="single"))
        for pair in itertools.combinations(MOTIFS3, 2):
            for arr in ("interleaved", "clustered"):
                arms.append(dict(label="+".join(pair), members=list(pair), K=K,
                                 arrangement=arr))
        for arr in ("interleaved", "clustered"):
            arms.append(dict(label="+".join(MOTIFS3), members=MOTIFS3, K=K,
                             arrangement=arr))
    R(f"{len(arms)} arms ({len(arms)-1} motif arms x 2 variants + 1 empty baseline)")

    # ---- one background per seed, shared by every arm -------------------------
    rng = np.random.default_rng(20261004)
    backgrounds = ["".join(rng.choice(list("ACGT"), size=CASSETTE,
                                      p=[BG["A"], BG["C"], BG["G"], BG["T"]]))
                   for _ in range(N_SEEDS)]

    jobs, meta = [], []
    for seed in range(N_SEEDS):
        bg = backgrounds[seed]
        for a in arms:
            variants = ("real",) if a["K"] == 0 else ("real", "scrambled")
            for variant in variants:
                assign = assignment(a["K"], a["members"], a["arrangement"]) if a["K"] else []
                cas, nb = build_cassette(bg, assign, variant, seed)
                jobs.append(bgseq[:lo] + cas + bgseq[lo + CASSETTE:])
                meta.append(dict(label=a["label"], K=a["K"], arrangement=a["arrangement"],
                                 variant=variant, seed=seed, bases_inserted=nb,
                                 n_members=len(a["members"]), cassette=cas))
    R(f"{len(jobs)} sequences to score "
      f"({N_SEEDS} seeds x {len(jobs)//N_SEEDS} conditions)")

    # ---- score in chunks, checkpointing as we go ------------------------------
    preds = []
    CH = 50
    for i in range(0, len(jobs), CH):
        preds.append(score(jobs[i:i + CH]))
        done = min(i + CH, len(jobs))
        el = time.time() - t0
        R(f"  scored {done}/{len(jobs)}  ({el:.0f}s elapsed, {el/done:.2f}s/seq, "
          f"eta {(len(jobs)-done)*el/done/60:.1f} min)")
    preds = np.concatenate(preds)

    df = pd.DataFrame(meta)
    df["pred"] = preds
    df.drop(columns=["cassette"]).to_csv(outdir / "gate1_additivity_raw.csv", index=False)
    pd.DataFrame([m for m in meta if m["seed"] == 0]).to_csv(
        outdir / "gate1_cassettes_seed0.csv", index=False)

    # ---- paired deltas vs the same seed's empty cassette ----------------------
    base = df[df.label == "empty"].set_index("seed")["pred"]
    df["delta"] = df["pred"] - df["seed"].map(base)

    rows = []
    for (label, K, arr), g in df[df.label != "empty"].groupby(
            ["label", "K", "arrangement"], sort=False):
        real = g[g.variant == "real"].sort_values("seed")
        scr = g[g.variant == "scrambled"].sort_values("seed")
        assert len(real) == len(scr) == N_SEEDS
        raw = real["delta"].to_numpy()
        sc = scr["delta"].to_numpy()
        adj = raw - sc                                   # paired by seed
        t = stats.ttest_1samp(adj, 0.0)
        rows.append(dict(
            condition=label, K_total=K, arrangement=arr, n_members=int(real["n_members"].iloc[0]),
            bases_inserted=int(real["bases_inserted"].iloc[0]), n_seeds=N_SEEDS,
            raw_delta=float(raw.mean()), raw_sem=float(raw.std(ddof=1) / np.sqrt(N_SEEDS)),
            scrambled_delta=float(sc.mean()), scrambled_sem=float(sc.std(ddof=1) / np.sqrt(N_SEEDS)),
            adjusted_delta=float(adj.mean()), adjusted_sem=float(adj.std(ddof=1) / np.sqrt(N_SEEDS)),
            adjusted_t=float(t.statistic), adjusted_p=float(t.pvalue),
            frac_of_required=float(adj.mean() / REQUIRED_DELTA)))
    res = pd.DataFrame(rows).sort_values(["K_total", "n_members", "condition", "arrangement"])
    res.to_csv(outdir / "gate1_additivity_summary.csv", index=False)

    R("\n" + "=" * 110)
    R("GATE 1 -- ADJUSTED DELTA BY CONDITION (adjusted = real - scrambled, paired within seed)")
    R("=" * 110)
    R(f"{'K':>3} {'condition':28s} {'arrange':12s} {'bp':>4} {'raw':>9} {'scram':>9} "
      f"{'adj':>9} {'+-sem':>7} {'p':>9} {'frac of 0.861':>14}")
    for _, r in res.iterrows():
        R(f"{r.K_total:3.0f} {r.condition:28s} {r.arrangement:12s} {r.bases_inserted:4.0f} "
          f"{r.raw_delta:+9.4f} {r.scrambled_delta:+9.4f} {r.adjusted_delta:+9.4f} "
          f"{r.adjusted_sem:7.4f} {r.adjusted_p:9.2e} {r.frac_of_required:14.3f}")

    # ---- additivity: observed combination vs sum of its singles ---------------
    R("\n" + "=" * 110)
    R("ADDITIVITY TEST -- observed combination vs the sum of its constituent singles")
    R("(expected = sum over members of single-motif adjusted delta scaled to that member's "
      "copy share)")
    R("=" * 110)
    add_rows = []
    for K in LOADS:
        singles = {r.condition: r.adjusted_delta
                   for _, r in res[(res.K_total == K) & (res.n_members == 1)].iterrows()}
        sing_sem = {r.condition: r.adjusted_sem
                    for _, r in res[(res.K_total == K) & (res.n_members == 1)].iterrows()}
        for _, r in res[(res.K_total == K) & (res.n_members > 1)].iterrows():
            members = r.condition.split("+")
            share = 1.0 / len(members)
            exp = sum(singles[m] * share for m in members)
            exp_sem = np.sqrt(sum((sing_sem[m] * share) ** 2 for m in members))
            ratio = r.adjusted_delta / exp if exp != 0 else np.nan
            add_rows.append(dict(K_total=K, condition=r.condition, arrangement=r.arrangement,
                                 observed=r.adjusted_delta, observed_sem=r.adjusted_sem,
                                 expected_additive=exp, expected_sem=float(exp_sem),
                                 obs_minus_exp=r.adjusted_delta - exp,
                                 obs_over_exp=float(ratio),
                                 frac_of_required=r.frac_of_required))
    add = pd.DataFrame(add_rows)
    add.to_csv(outdir / "gate1_additivity_vs_expected.csv", index=False)
    R(f"{'K':>3} {'condition':28s} {'arrange':12s} {'observed':>10} {'expected':>10} "
      f"{'obs-exp':>9} {'obs/exp':>8} {'frac req':>9}")
    for _, r in add.iterrows():
        R(f"{r.K_total:3.0f} {r.condition:28s} {r.arrangement:12s} {r.observed:+10.4f} "
          f"{r.expected_additive:+10.4f} {r.obs_minus_exp:+9.4f} {r.obs_over_exp:8.3f} "
          f"{r.frac_of_required:9.3f}")

    # ---- spacing ---------------------------------------------------------------
    R("\n" + "=" * 110)
    R("SPACING -- interleaved vs clustered, paired within seed")
    R("=" * 110)
    sp_rows = []
    for (label, K), g in df[(df.arrangement.isin(["interleaved", "clustered"]))
                            & (df.variant == "real")].groupby(["label", "K"], sort=False):
        a = g[g.arrangement == "interleaved"].sort_values("seed")["delta"].to_numpy()
        b = g[g.arrangement == "clustered"].sort_values("seed")["delta"].to_numpy()
        gs = df[(df.label == label) & (df.K == K) & (df.variant == "scrambled")]
        sa = gs[gs.arrangement == "interleaved"].sort_values("seed")["delta"].to_numpy()
        sb = gs[gs.arrangement == "clustered"].sort_values("seed")["delta"].to_numpy()
        d = (a - sa) - (b - sb)
        t = stats.ttest_rel(a - sa, b - sb)
        sp_rows.append(dict(condition=label, K_total=K, n_seeds=N_SEEDS,
                            interleaved_adj=float((a - sa).mean()),
                            clustered_adj=float((b - sb).mean()),
                            difference=float(d.mean()),
                            difference_sem=float(d.std(ddof=1) / np.sqrt(N_SEEDS)),
                            t=float(t.statistic), p=float(t.pvalue)))
    sp = pd.DataFrame(sp_rows).sort_values(["K_total", "condition"])
    sp.to_csv(outdir / "gate1_spacing.csv", index=False)
    R(f"{'K':>3} {'condition':28s} {'interleaved':>12} {'clustered':>11} {'diff':>9} "
      f"{'+-sem':>7} {'p':>9}")
    for _, r in sp.iterrows():
        R(f"{r.K_total:3.0f} {r.condition:28s} {r.interleaved_adj:+12.4f} "
          f"{r.clustered_adj:+11.4f} {r.difference:+9.4f} {r.difference_sem:7.4f} {r.p:9.3f}")

    best = res.loc[res.adjusted_delta.idxmax()]
    R("\n" + "=" * 110)
    R(f"  best condition overall : {best.condition} (K={best.K_total:.0f}, {best.arrangement})")
    R(f"  adjusted delta         : {best.adjusted_delta:+.4f} +- {best.adjusted_sem:.4f} (n={N_SEEDS} seeds)")
    R(f"  required (median->p95) : {REQUIRED_DELTA:+.4f}")
    R(f"  fraction of required   : {best.adjusted_delta/REQUIRED_DELTA:.3f}x")
    R("=" * 110)

    json.dump(dict(n_seeds=N_SEEDS, loads=LOADS, required_delta=REQUIRED_DELTA,
                   neutral_context=f"{NEUTRAL[0]}:{NEUTRAL[1]}", batch=batch,
                   n_sequences=len(jobs), seconds=time.time() - t0,
                   best=dict(condition=best.condition, K_total=int(best.K_total),
                             arrangement=best.arrangement,
                             adjusted_delta=float(best.adjusted_delta),
                             adjusted_sem=float(best.adjusted_sem),
                             frac_of_required=float(best.adjusted_delta / REQUIRED_DELTA)),
                   summary=res.to_dict("records"),
                   additivity=add.to_dict("records"),
                   spacing=sp.to_dict("records")),
              open(outdir / "gate1_results.json", "w"), indent=2)
    Path(OUT, "results", "GATE1_LOG.txt").write_text("\n".join(LOG))
    R(f"\nGATE1_DONE in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
