# Cross-mechanism transfer experiment — expanded multi-study set

**Date:** 2026-08-24 · pod `<pod-address-redacted>` (RTX 4090, 24 GB)
**Verdict: cross-mechanism transfer FAILS in both directions.**
**Second shuffled-label control failure. Stopped and reported rather than interpreted.**

---

## STEP 2 precondition — verified, not assumed

The corrected targets are signed (range −14.83 … +16.18). Verified on a real 512-gene batch:

```
act_func=None head:  pred range [-0.3773, 0.2832]
                     negatives 1402/2048 (68.5%)
                     head activation attribute: None
VERDICT: PASS -- head emits negative values
```

**Correction to the recorded contract:** `uns['contract']['head']` says "ConvHead + softplus
(non-negative outputs)", but the previous run's `make_head` already passed `act_func=None`.
That contract string appears to have described Borzoi's *native* head, not the fine-tune head.
No config change was needed — but it was checked empirically, and the script was written to
abort if the head could not emit negatives.

---

## STEP 1 — cache pass

**Caching identity verified as a hard preflight** (aborts the run on failure, so no GPU hours
are spent on an invalid cache):

| gene | max\|full − cached\| / max\|full\| |
|---|---|
| A1BG | 6.51e-08 |
| GON7 | 1.13e-07 |
| PTP4A3 | 1.85e-07 |
| ZZEF1 | 1.28e-07 |

Worst 1.85e-07 vs 1e-04 tolerance → **PASSED**. Confirms the 1×1 conv and average pool
commute, so the pooled cache is algebraically exact, not a lossy summary.

| | |
|---|---|
| Genes | 13,211 (from 19,652) |
| Wall clock | **1.66 h** (projected 1.7 h) |
| Peak GPU | 3.77 GB |
| Cache size | 101.5 MB forward + 101.5 MB RC |
| All-zero rows | 0 fwd, 0 rc |
| NaNs | 0 |
| IUPAC sanitised | **31 windows, 70 characters** of 6,926,368,768 scanned |
| Codes seen | Y 27, R 18, M 14, W 6, K 3, B 2 |

(Previous run: 49 windows / 133 chars over 19,652 genes — same rate, fewer genes.)

---

## STEP 3 — results

All AUCs use the **pre-registered threshold of 1.0**, never tuned. Evaluation is on held-out
**chr9 + chr6** throughout. Ambiguous-gene counts are reported on every AUC, as required.

### Headline table (test = chr9 + chr6)

| Condition | test AUC | n_used / 1164 | ambiguous | pos / neg | shuffled mean ± sd | # beat real | empirical p | untrained floor |
|---|---|---|---|---|---|---|---|---|
| **Ceiling: within-OIS** | **0.9079** | 23 (2.0%) | 1141 | 19 / 4 | 0.4697 ± 0.2653 | 1/20 | 0.0952 | 0.5395 |
| **Ceiling: within-replicative** | **0.7826** | 82 (7.0%) | 1082 | 37 / 45 | 0.4796 ± 0.2050 | 0/20 | **0.0476** | 0.6511 |
| **A: OIS → replicative** | **0.4979** | 82 (7.0%) | 1082 | 37 / 45 | 0.5041 ± 0.0782 | **11/20** | **0.5714** | 0.6192 |
| **B: replicative → OIS** | **0.5789** | 23 (2.0%) | 1141 | 19 / 4 | 0.4447 ± 0.1665 | **4/20** | **0.2381** | 0.5132 |
| A matched (9 ON / 12 OFF) | 0.5249 | 82 (7.0%) | 1082 | 37 / 45 | — | — | — | — |
| B matched (9 ON / 12 OFF) | 0.5658 | 23 (2.0%) | 1141 | 19 / 4 | — | — | — | — |

### All-chromosome supplementary (larger n, more stable)

| Condition | AUC | n_used / 13,211 |
|---|---|---|
| within-OIS | 0.8668 | 316 |
| within-replicative | 0.8571 | 1,040 |
| A: OIS → replicative | 0.6014 | 1,040 |
| B: replicative → OIS | 0.5816 | 316 |
| A matched | 0.6049 | 1,040 |
| B matched | 0.5821 | 316 |

Within-mechanism all-chromosome numbers are **inflated by training-set contamination** (the
head saw those genes' targets). The transfer numbers are cleaner — the head never saw the
evaluation mechanism's labels at any locus — but they are still not a clean held-out split.

---

## Answers to the questions asked

**Does cross-mechanism transfer work?** **No. It fails cleanly in both directions.**

**Does either direction reach 0.85?** No. Transfer reaches 0.4979 (A) and 0.5789 (B) on
held-out chromosomes; 0.6014 and 0.5816 on all chromosomes. The bar is 0.85.

**Is either distinguishable from its shuffled control?** **No.** Direction A: p = 0.5714, with
11 of 20 random relabelings beating the real model — and its AUC of 0.4979 sits *below* its own
untrained-head floor of 0.6192. Direction B: p = 0.2381, 4 of 20 beat it.

**How much worse is transfer than within-mechanism?** On the like-for-like all-chromosome
comparison, **≈0.26–0.27 AUC**: 0.8668 → 0.6014 (OIS→rep) and 0.8571 → 0.5816 (rep→OIS).

**Is the asymmetry between A and B biology?** **No — and it was tested rather than assumed.**
Matching training volume to 9 ON / 12 OFF in both directions moved A from 0.4979 to 0.5249 and
B from 0.5789 to 0.5658. Both remain at chance. The A-vs-B gap is not explained by OIS
contributing 15 of 24 ON tasks, and neither direction becomes meaningful when volume is equalised.

---

## 🔴 Second shuffled-label control failure — reported, not explained away

**Three of the four conditions fail the shuffled-label control at p < 0.05:**

| Condition | p | verdict |
|---|---|---|
| within-replicative | 0.0476 | marginally passes (0/20 beat) |
| within-OIS | 0.0952 | **fails** (1/20 beat) |
| B: replicative → OIS | 0.2381 | **fails** (4/20 beat) |
| A: OIS → replicative | 0.5714 | **fails badly** (11/20 beat) |

Per the pre-registered stop rule, I am not interpreting these as evidence that the model
learned senescence biology. The previous single-dataset run failed here at p = 0.143; this is
the second failure, across a different data build, a different target definition, and a
corrected batch structure.

The permutation distributions are also very wide (sd 0.17–0.27) — a mechanical consequence of
the tiny evaluable gene sets below, not a property of the model.

---

## 🔴 The power problem that limits every OIS-evaluated number

Only **23 of 1,164** test genes are evaluable against the OIS contrast, and just **4 of them
are negatives**. An AUC on 4 negatives is not a stable quantity. Diagnosed directly:

| Contrast | genes with \|log2FC\| ≥ 1 | on test chroms | sd of log2FC |
|---|---|---|---|
| OIS pooled (15 v 12) | 316 / 13,211 | 23 / 1,164 | 0.384 |
| — GSE205692 alone (12 v 9) | 178 | 14 | **0.293** |
| — GSE74324 alone (3 v 3) | **1,979** | 182 | **0.902** |
| replicative (9 v 15) | 1,040 | 82 | 0.602 |

**The GSE205692 microarray compresses the OIS contrast and starves the evaluation.** It
contributes 12 of 15 OIS ON tasks and 9 of 12 OFF, so it dominates the pooled OIS mean, and its
dynamic range is ~3× narrower than the RNA-seq study it is pooled with. This was flagged as a
risk in the build report ("largest ON contributor and among the weakest by effect size"); it has
now materialised as a concrete measurement problem.

**Consequences for reading the table:**
- The within-OIS "ceiling" of 0.9079 rests on 23 genes and should not be trusted as a ceiling.
- Direction B (0.5789) is measured on the same 23 genes and is equally fragile.
- **Direction A is the trustworthy result**: evaluated on 82 genes with a balanced 37/45
  split, it came out at **0.4979 — indistinguishable from chance, below its own untrained
  floor, with 11 of 20 shuffled permutations beating it.** That single number is the most
  defensible evidence in the experiment, and it says transfer does not happen.

---

## Study-discrimination control — the one clear pass

"Permute study assignment holding senescence status fixed" was implemented as the question it
exists to answer: can the head reproduce a **study** contrast as well as a **senescence**
contrast, with senescence status balanced on both sides so anything detected is batch?

| | AUC | n_used |
|---|---|---|
| within-OIS senescence contrast | 0.9079 | 23 |
| within-OIS **study** contrast (GSE205692 vs GSE74324, status-balanced) | **0.5803** | 1,042 |

**Residual batch structure is much weaker than the senescence signal.** The per-study centering
is doing its job — this is not a repeat of the "model learns which study generated the sample"
failure mode.

**Not computable for the replicative direction:** GSE175533 is a single study, so no study
contrast exists within it. Reported as not computable rather than substituted with something else.

---

## CDKN1A (chr6, held out) — in every condition

| Trained on | predicted contrast | measured (replicative) | measured (OIS) | direction |
|---|---|---|---|---|
| OIS | **+0.4883** | +0.5146 | +0.2765 | **CORRECT** both ways |
| replicative | **−0.1441** | +0.5146 | +0.2765 | **WRONG** both ways |
| OIS (matched) | +0.4968 | +0.5146 | +0.2765 | CORRECT |
| replicative (matched) | −0.1478 | +0.2765 | +0.2765 | WRONG |

The predicted value depends only on the trained head, so the OIS-trained head gives the same
+0.4883 whether it is being read as "within-OIS" or "transfer to replicative".

**The replicative-trained head still gets CDKN1A wrong** — −0.1441 here, against −0.019 in the
previous run. Same failure, slightly larger magnitude, now with the corrected multi-study
targets. **The OIS-trained head gets it right, and lands remarkably close to the measured
replicative value (+0.4883 vs +0.5146).** That is a genuine asymmetry, but given that
Direction A is at chance overall, one gene agreeing is not evidence of transfer — it is one
gene.

---

## Untrained-head floors are high, and that matters

| Condition | untrained floor |
|---|---|
| within-OIS | 0.5395 |
| within-replicative | 0.6511 |
| A: OIS → replicative | 0.6192 |
| B: replicative → OIS | 0.5132 |

A randomly initialised head scores **0.62–0.65** on the replicative contrast. Any AUC in the
0.5–0.65 band on this metric is therefore not evidence of anything. Direction A's 0.4979 is
below its own floor.

---

## What this means for the project

The build report measured cross-study Spearman of the per-gene ON−OFF response at **+0.531**
between the two OIS studies but only **+0.142** and **+0.093** between replicative and either
OIS study. **That finding reproduces at the model level.** A head trained on one induction
mechanism does not predict the other's response — in the one direction with adequate evaluation
power, it performs at chance.

Stated plainly, as the brief asked: **cross-mechanism transfer fails.** "Senescence" as
currently constructed is not one learnable target state. If that holds up, the project's design
target, beachhead, and validation plan need narrowing to a specific induction mechanism rather
than senescence in general.

Two caveats that keep this from being the final word:

1. **The shuffled-label control failed again**, in 3 of 4 conditions. Under the pre-registered
   rule this means the within-mechanism numbers are *also* not established — the failure is not
   confined to the transfer directions. What is established is the *negative*: nothing here is
   distinguishable from chance except, marginally, within-replicative.
2. **The OIS evaluation is underpowered by construction** because the microarray dominates that
   arm. The clean fix is not more modelling — it is either dropping GSE205692 from the OIS
   contrast (leaving GSE74324's 3 v 3, which alone yields 1,979 evaluable genes) or finding an
   RNA-seq OIS study with proliferating controls to replace it.

**Recommended next step, cheapest first:** re-run the transfer with the OIS arm defined by
**GSE74324 only** (3 ON / 3 OFF, RNA-seq, sd of log2FC 0.902, 1,979 evaluable genes). That trades
task count for evaluation power and would make both directions measurable on a comparable gene
set. It needs no new cache — the existing one covers all 13,211 genes — so it costs minutes,
not hours.

---

## Cost

| | |
|---|---|
| Projected | 3.0 h ≈ $2.25 |
| Actual | **~2.3 h ≈ $1.72** |
| Breakdown | setup 0.42 h · cache 1.66 h · transfer 0.005 h (0.3 min) · retrieval 0.1 h |

Transfer training was far cheaper than projected: each head trains in 1–2 s on cached
embeddings, and the shuffled-label control needed only one training run rather than 20 (see below).

**Control-1 efficiency note:** shuffling which tasks count as ON/OFF does not change the
training targets — the head predicts every task — so the retrained head is bit-identical. This
was *verified* (`retrained_head_bit_identical = True` in all four conditions, matching the
previous run's observation) and the trained head then reused across permutations. Same control,
1/20th the cost.

---

## Artifacts (all verified byte-for-byte against remote)

`~/Downloads/Senescence/infra_borzoi_grelu/execution_results_transfer/`

| File | Bytes |
|---|---|
| `cache_expanded/cache_fwd.npy` | 101,460,608 |
| `cache_expanded/cache_rc.npy` | 101,460,608 |
| `cache_expanded/caching_identity_check.json` | 222 |
| `cache_expanded/chunks_done.txt` | 70 |
| `cache_expanded/sanitize_stats.json` | 325 |
| `transfer_results/results_transfer.json` | 20,220 |
| `cache_expanded.log` | 2,749 |
| `transfer.log` | 8,305 |
| `cache_pass_expanded.py` | 8,384 |
| `transfer_experiment.py` | 21,371 |

No trained head checkpoints were written — heads train in 1–2 s from the cache, so the cache is
the expensive artifact and it is secured. Total local: 225 MB.
