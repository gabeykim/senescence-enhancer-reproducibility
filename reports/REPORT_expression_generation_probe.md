# RNA-seq-only cross-mechanism analysis + generation-readiness probe

**Date:** 2026-08-24 · pod `<pod-address-redacted>` (RTX 4090) · **cost ~0.77 h ≈ $0.58**
**Cache reused, not rebuilt** — verified (13211, 1920) float32 ×2, all finite, 0 all-zero rows, fwd ≠ rc.

**Verdict in one line: pool the data, keep mechanism-specific contrasts, and treat the
model as classification-capable but NOT generation-ready.**

---

## 0. Preconditions

| Check | Result |
|---|---|
| Cache shape / dtype | (13211, 1920) float32, both fwd and rc |
| Finite | yes, both |
| All-zero rows | 0 fwd, 0 rc |
| fwd ≠ rc | yes |
| Prior caching-identity record | worst 1.85e-07 vs 1e-04 tolerance (travelled with the cache) |
| Head emits negatives (signed targets) | **PASS** — range [−0.3773, +0.2832], 1402/2048 (68.5%) negative, `act_func=None` |
| TF32 | disabled on cudnn and matmul |

Setup ~25 min, genome 194 contigs renamed/reindexed, chr1 = 248,956,422.

---

## 1. Dropping the microarray fixed the power problem

That was the point of the rerun, and it worked:

| | previous run (pooled w/ GSE205692) | this run (GSE74324 only) |
|---|---|---|
| OIS evaluable test genes | **23** (4 negatives) | **182** (57 pos / 125 neg) |
| replicative evaluable test genes | 82 | 82 |

**No condition is underpowered** by the ~50-gene threshold. Every OIS-side number in the
previous report rested on 23 genes; those are now superseded.

---

## 2. All eight conditions (test = chr9 + chr6, threshold 1.0, 20 permutations)

| Condition | AUC | evaluable | pos/neg | shuffled mean ± sd | range | beat | p | control | untrained floor |
|---|---|---|---|---|---|---|---|---|---|
| within_OIS (3v3) | 0.7854 | 182 | 57/125 | 0.5034 ± 0.2782 | [0.209, 0.802] | 5/20 | 0.2857 | **FAIL** | 0.4418 |
| **pooled, OIS-specific contrast** | **0.7843** | **182** | 57/125 | 0.4866 ± 0.1498 | [0.212, 0.729] | 0/20 | **0.0476** | **PASS** | 0.5485 |
| within_replicative (9v15) | 0.7826 | 82 | 37/45 | 0.4796 ± 0.2050 | [0.205, 0.779] | 0/20 | **0.0476** | **PASS** | 0.6511 |
| pooled, rep-specific contrast | 0.7808 | 82 | 37/45 | 0.5077 ± 0.2086 | [0.208, 0.770] | 0/20 | **0.0476** | **PASS** | 0.5664 |
| pooled, unified contrast → rep | 0.7592 | 82 | 37/45 | 0.5065 ± 0.1845 | [0.208, 0.757] | 0/20 | **0.0476** | **PASS** | 0.6180 |
| pooled, unified contrast → OIS | 0.6895 | 182 | 57/125 | 0.4939 ± 0.1485 | [0.231, 0.744] | 2/20 | 0.1429 | **FAIL** | 0.4967 |
| B: replicative → OIS | 0.6042 | 182 | 57/125 | 0.4865 ± 0.0968 | [0.293, 0.618] | 3/20 | 0.1905 | **FAIL** | 0.5708 |
| A: OIS → replicative | 0.5135 | 82 | 37/45 | 0.4992 ± 0.0191 | [0.458, 0.542] | 3/20 | 0.1905 | **FAIL** | 0.5784 |

**⚠️ Resolution limit, stated up front:** with 20 permutations the smallest achievable
empirical p is 1/21 = **0.0476**. Every "PASS" above is exactly at that floor — it means
*zero of twenty* permutations beat the real model, which is the strongest statement 20
permutations can make and no stronger. These are not p = 0.001 results. Confirming them
properly needs ~200 permutations, which costs minutes.

---

### Q1 — Which mechanism produces a head that beats its shuffled control?

**Replicative does, on its own.** OIS on its own does not — **but that failure is not
informative**, and it would be a mistake to read it as "no OIS signal."

With a 3v3 design there are only C(6,3) = 20 distinct ON/OFF splits, so 20 permutations
essentially *enumerate the entire null*, and several necessarily land on or near the true
split and its complement. That is exactly what the numbers show: within_OIS has the widest
shuffled distribution of any condition (sd 0.2782, max 0.8024). **A 3v3 design cannot
support a meaningful permutation test at all.** Its AUC of 0.7854 is the highest of the
eight; its control result is uninformative in both directions.

The way to get a testable OIS contrast is to borrow task count from the pooled set, which
is what the next answer is about.

### Q2 — Does pooling help or hurt?

**It depends entirely on whether you pool the DATA or pool the CONTRAST, and the two go in
opposite directions.**

| | OIS | replicative |
|---|---|---|
| trained on that mechanism alone | 0.7854 (control FAILS) | 0.7826 (PASS) |
| **pooled data, mechanism-specific contrast** | **0.7843 (PASS)** | **0.7808 (PASS)** |
| pooled data, unified ON-vs-OFF contrast | 0.6895 (**FAIL**) | 0.7592 (PASS) |

- **Pooling the data is neutral-to-helpful.** AUC is essentially unchanged (0.7843 vs
  0.7854; 0.7808 vs 0.7826), but the extra tasks widen the permutation space enough that
  OIS's control becomes testable and passes.
- **Pooling the contrast hurts, and it hurts precisely where the data build predicted.**
  The *same pooled head* scores 0.7843 on OIS with an OIS-specific contrast and drops to
  **0.6895, failing its control**, when forced through a single unified ON-vs-OFF axis.
  Replicative barely moves (0.7808 → 0.7592) — consistent with replicative supplying 9 of
  the 12 pooled ON tasks, so the unified axis is mostly the replicative axis.

**This is the ρ ≈ 0.1 finding reproducing at the model level, and it is now actionable:
train on both mechanisms, but do not collapse them to one senescence score.** A shared
trunk representation serves both; a shared output contrast does not.

**Transfer still fails in both directions** — A: 0.5135 (below its own 0.5784 untrained
floor), B: 0.6042, both p = 0.1905. This is now a well-powered negative rather than the
previous run's broken one.

### CDKN1A — third consecutive repeat of the same pattern

| Head / contrast | predicted | measured | direction |
|---|---|---|---|
| within_OIS | +0.5694 | +1.5039 (OIS) | **CORRECT** |
| pooled, OIS-specific | +0.5856 | +1.5039 (OIS) | **CORRECT** |
| within_replicative | −0.1441 | +0.5146 (rep) | **WRONG** |
| pooled, rep-specific | −0.1522 | +0.5146 (rep) | **WRONG** |
| pooled, unified → OIS | +0.0006 | +1.5039 (OIS) | correct sign, ~zero magnitude |

Every replicative-contrast head predicts CDKN1A in the wrong direction — −0.019, then
−0.1441, now −0.1522 across three runs. Every OIS-contrast head gets it right. Whatever
this is, it is systematic and reproducible, not noise.

### Study-discrimination control

Computable only for the pooled conditions (the single-study ones have no study contrast):
**study AUC 0.5343** (263 evaluable, status-balanced) against senescence AUC 0.7843 /
0.7808. Residual batch structure is well below the senescence signal — the per-study
centering is holding.

---

## 3. GENERATION-READINESS PROBE

Run against **pooled_mechspecific_OIS** — the best condition that passes its control
(AUC 0.7843, p = 0.0476, 182 evaluable).

> *Selection correction:* the rerun script's own selector examined only 4 of the 8
> conditions and chose within_replicative. Re-running the selection over all 8 picked
> pooled_mechspecific_OIS instead — marginally higher AUC, same p, and 2.2× the evaluable
> genes. The probe used the corrected choice.

### (c) Dynamic range — adequate

| | min | 5% | 50% | 95% | max | sd | IQR |
|---|---|---|---|---|---|---|---|
| predicted contrast | −2.822 | −0.994 | −0.135 | +0.510 | +1.306 | **0.441** | 0.517 |
| measured log2FC | −6.045 | — | — | — | +7.498 | **0.868** | — |

Predicted sd is **0.508×** the measured sd — the model compresses the range about 2×, but
it does span a usable spread. **This is not a degenerate narrow-band predictor.** Candidates
can be ranked.

### (a) In-silico saturation mutagenesis — concentrated, but the motif claim is unsupported

200 bp around each TSS, every position × 3 alternative bases (601 forward passes per
region), on the 20 most-responsive held-out genes. **13 of 20 regions completed** within
the 1,500 s budget (trunk batch size 4, peak 22.29 GB) — reported as run, not as planned.

| Metric | Value | Interpretation |
|---|---|---|
| Top-10% concentration | **0.291** (diffuse null = 0.10) | ~2.9× more concentrated than uniform — **not diffuse** |
| Mean \|Δ\| per single mutation | 0.00733 | 1.7% of the prediction sd |
| Max \|Δ\| per single mutation | 0.194 | 44% of the prediction sd — some positions matter a lot |
| Motif enrichment (in-motif / out-motif importance) | 2.151 ± 1.501 | **n = 2 regions, p = 0.474 — NOT SUPPORTED** |

**🔴 The headline motif question could not be answered.** Only **2 of 13** windows contained
any consensus match for NF-κB / C/EBPβ / AP-1 / ETS at all, so the enrichment statistic
rests on two regions and is not significant. That is a limitation of scanning strict IUPAC
consensus patterns in a 200 bp window — not evidence that attributions ignore motifs. What
*can* be said is that attributions are spatially concentrated rather than diffuse, which is
the necessary precondition; the specific attribution-to-motif link remains untested.

**Per-gene sign agreement on these 13 most-responsive genes: 10/13 (76.9%), Spearman 0.533.**
Three of the most strongly-responding held-out genes get the wrong sign (COL15A1 measured
−6.04 → predicted +0.11; FOXQ1 +5.08 → −0.004; COL10A1 +3.73 → −0.26).

### (b) Motif-insertion response — a real, motif-specific, dose-dependent gradient

k copies of each motif written into a random 200 bp cassette at the centre of a neutral
held-out context (ACTL7A, chr9, measured log2FC +0.0000). k ∈ {0,1,2,3,4,6,8}, 5 random
backgrounds each. Scrambled motifs (same composition, shuffled order) as the control.

| Motif | variant | ρ(pred, k) | slope/copy | Δ (k=0→8) |
|---|---|---|---|---|
| C/EBPβ | **real** | **+0.818** | 0.01701 | **+0.13181** |
| C/EBPβ | scrambled | −0.034 | −0.00002 | −0.00173 |
| NF-κB p65 | **real** | **+0.717** | 0.01283 | **+0.08855** |
| NF-κB p65 | scrambled | +0.386 | 0.00047 | +0.00357 |
| AP-1 | **real** | **+0.600** | 0.01714 | **+0.12693** |
| AP-1 | scrambled | +0.229 | 0.00015 | +0.00075 |
| ETS | **real** | +0.530 | 0.00251 | +0.01833 |
| ETS | scrambled | +0.180 | 0.00028 | +0.00153 |

**This is the clearest positive result in the probe.** All four real motifs produce a
monotonic increase in predicted senescence activity; all four scrambled controls produce
essentially nothing (|Δ| ≤ 0.0036). The response is **sequence-specific, not composition
artifact**. Directed evolution would have a real gradient to climb.

**But the lever is weak.** The largest achievable effect — 8 copies of C/EBPβ, i.e. a
saturated 200 bp cassette — is **+0.132, only 29.9% of the prediction sd** across held-out
sequences (0.441). Moving a designed sequence from the median (−0.135) to the 95th
percentile (+0.510) requires ≈ +0.65; maximal motif engineering delivers ≈ +0.13, about a
fifth of that.

### Verdict: **classification-only — not generation-ready**

Stated plainly, against the three options:

- **Not "neither."** The signal is real: predictions span a usable range, attributions are
  ~2.9× concentrated rather than diffuse, and motif insertion produces a specific,
  monotonic, dose-dependent response that scrambled controls do not.
- **Not "generation-ready."** The maximum achievable design effect is ~0.3 sd. A directed-
  evolution loop would climb a genuine but shallow gradient and plateau well short of the
  range the model already assigns to natural sequences. Per-gene sign agreement is 77% on
  the *most responsive* held-out genes, so the objective being climbed is itself unreliable
  at the level of individual sequences. And the attribution-to-motif link — the thing that
  would let you reason about *what* to build rather than just hill-climb — is untested
  (n = 2).

**Classification-capable at the resolution limit of a 20-permutation test; usable as a
ranker; not yet usable as a design objective.**

---

## 4. What I would do next, cheapest first

1. **Re-run the controls with 200 permutations** (~minutes, no new cache). Every current
   PASS sits at p = 0.0476, the floor of a 20-permutation test. This is the single cheapest
   way to convert "passes at the resolution limit" into a real confidence statement, and it
   should be done before any decision rests on these results.
2. **Re-run the ISM motif analysis with PWM scanning rather than strict consensus**, over a
   wider window (1–2 kb) — the motif question failed on scanning stringency, not on the
   model. Also cheap; the ISM importances are already saved for the 13 regions.
3. **If both hold up, the design target is OIS with a mechanism-specific contrast**, trained
   on pooled data. That is the only configuration here that combines the highest AUC, a
   passing control, adequate evaluation power, and correct CDKN1A direction.
4. **Do not build a unified senescence score.** Two independent lines of evidence now say
   the mechanisms do not share an axis: ρ ≈ 0.1 in the data, and 0.7843 → 0.6895 with a
   control failure when the contrast is unified.

The accessibility/enhancer route is not yet forced. Transfer failed, but a per-mechanism
head does clear its control, and the motif-insertion gradient is real. What is not
established is that the gradient is strong enough to design against.

---

## 5. Artifacts

All at **`/Users/gabeykim/Downloads/Senescence/infra_borzoi_grelu/execution_results_generation_probe/`**
— every file verified byte-for-byte against the remote before the pod was released.

| File | Bytes |
|---|---|
| `results_rerun.json` | 26,459 |
| `results_generation_probe.json` | 10,803 |
| `heads.pt` (within_OIS, within_replicative, pooled) | 465,213 |
| `motif_insertion.csv` | 10,570 |
| `ism_per_region.csv` | 1,653 |
| `ism_importance.npz` (13 regions × 200 positions) | 15,225 |
| `task_indices.npz` | 1,964 |
| `rerun.log` | 10,149 |
| `probe.log` | 4,654 |
| `rerun_transfer_rnaseq.py` | 17,893 |
| `generation_probe.py` | 19,439 |
| `setup_pod.sh` | 2,955 |

The trunk cache was **not** re-copied — it is unchanged and already held at
`infra_borzoi_grelu/execution_results_transfer/cache_expanded/` (2 × 101,460,608 bytes).

**Cost:** ~0.77 h × $0.75 = **$0.58**. Setup 25 min, rerun 0.5 min, probe ~28 min
(ISM dominated).
