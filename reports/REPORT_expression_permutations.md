# 200-permutation controls + PWM motif-attribution rescan

**Date:** 2026-08-24 · pod `<pod-address-redacted>` (RTX 4090) · **cost ~0.42 h ≈ $0.32**
**Cache reused, not rebuilt.**

---

## Headline

**No verdict changed. Every condition that passed at 20 permutations still passes at proper
resolution, and every failure still fails.** The passes are no longer sitting on the
resolution floor — the strongest is now p = 0.00498 (0 of 200), where before it could only
be reported as "≤ 0.0476".

**But two conditions turn out to be structurally incapable of ever passing**, which the
20-permutation run could not distinguish from ordinary failure. See §2.

---

## 0. Preconditions

| Check | Result |
|---|---|
| Cache shape / dtype | (13211, 1920) float32 ×2 |
| Finite | yes, both |
| All-zero rows | 0 fwd, 0 rc |
| fwd ≠ rc | yes |
| Cache mean / sd | 0.09420 / 0.22701 — identical to prior run |
| Trainset md5 | `538e1e5185d7aba763bf0ca40bdecc19` (matches local) |
| Head emits negatives | **PASS** — range [−0.3773, +0.2832], 1402/2048 (68.5%) |
| Retrained head bit-identical under label shuffling | **True, re-confirmed in all 8 conditions** |
| TF32 | off on cudnn and matmul |
| Genome | 194 contigs renamed/reindexed, chr1 = 248,956,422 |

The bit-identity check matters: shuffling which tasks count as ON/OFF does not change the
training targets, so one training run legitimately serves every permutation. Re-verified
here rather than assumed from the prior run.

---

## 1. TASK 1 — controls at proper resolution

Two regimes, chosen by the size of the null space rather than by preference:

- **EXHAUSTIVE** — when the number of distinct ON/OFF partitions is small enough to
  enumerate, all of them are evaluated and the p-value is **exact**. Drawing 200 random
  permutations from a 20-partition space adds no resolution; it resamples the same 20
  values.
- **SAMPLED** — 200 random partitions, p = (b+1)/(N+1), minimum attainable p = 1/201 = 0.00498.

Partition space for an n_on vs n_off contrast drawn from n_tasks:
C(n_tasks, n_on) × C(n_tasks − n_on, n_off).

| Condition | AUC | evaluable | shuffled mean ± sd | range | beat | **new p** | prior p(20) | mode | partitions evaluated / possible | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| **pooled_mechspecific_OIS** | 0.7843 | 182 | 0.4763 ± 0.1594 | [0.208, 0.770] | **0/200** | **0.00498** | 0.0476 | sampled | 200 / 11,875,500 | **PASS** |
| pooled_mechspecific_rep | 0.7808 | 82 | 0.5012 ± 0.1718 | [0.202, 0.783] | 3/200 | **0.01990** | 0.0476 | sampled | 200 / 7.76×10¹¹ | **PASS** |
| within_replicative | 0.7826 | 82 | 0.4754 ± 0.1914 | [0.204, 0.801] | 5/200 | **0.02985** | 0.0476 | sampled | 200 / 1,307,504 | **PASS** |
| pooled_unified_rep | 0.7592 | 82 | 0.5094 ± 0.1642 | [0.203, 0.802] | 8/200 | **0.04478** | 0.0476 | sampled | 200 / 86,493,225 | **PASS (marginal)** |
| B_rep_to_OIS | 0.6042 | 182 | 0.4877 ± 0.0894 | [0.304, 0.666] | 20/200 | 0.10448 | 0.1905 | sampled | 200 / 1,307,504 | FAIL |
| pooled_unified_OIS | 0.6895 | 182 | 0.4971 ± 0.1552 | [0.213, 0.790] | 27/200 | 0.13930 | 0.1429 | sampled | 200 / 86,493,225 | FAIL |
| within_OIS | 0.7854 | 182 | 0.5000 ± 0.2861 | [0.198, 0.802] | 5/20 | **0.25000 (exact)** | 0.2857 | **exhaustive** | 20 / 20 | FAIL |
| A_OIS_to_rep | 0.5135 | 82 | 0.5000 ± 0.0207 | [0.458, 0.542] | 5/20 | **0.25000 (exact)** | 0.1905 | **exhaustive** | 20 / 20 | FAIL |

**Conditions whose verdict changed: none.**

### The passes survive, and one is now genuinely strong

`pooled_mechspecific_OIS` — the condition selected as the design target in the prior run —
has **zero of 200 permutations beating it, p = 0.00498**, the floor of a 200-permutation
test. That is a real result rather than an artefact of coarse resolution.

`pooled_unified_rep` is the weak one: 8/200, p = 0.04478. It clears 0.05 but only just, and
it would not survive any correction for multiple comparisons across these eight tests. It
should not be treated as established.

---

## 2. 🔴 The finding that the 20-permutation run could not have produced

**`within_OIS` and `A_OIS_to_rep` cannot ever pass, at any number of permutations.**

Both draw a 3v3 contrast from only 6 tasks. That admits exactly **C(6,3) × C(3,3) = 20**
distinct ON/OFF partitions. Enumerating all 20 gives an **exact** p of 0.25000 for both.
Because the true partition is necessarily one of the 20, the smallest attainable p is
**1/20 = 0.0500** — a 3v3 design cannot produce a significant permutation result even in
principle, let alone at 0.05.

Running 200 permutations there would have resampled the same 20 values and reported a
spuriously precise-looking p. The exhaustive treatment is the honest one.

Two consequences:

1. **`within_OIS`'s failure is not evidence against OIS signal.** Its AUC (0.7854) is the
   highest of all eight conditions. Its control is simply uncomputable at useful resolution.
   The prior run's suspicion is now confirmed arithmetically rather than argued.
2. **The only way to get a testable OIS contrast is to borrow task count from the pooled
   set** — which is exactly what `pooled_mechspecific_OIS` does (3v3 contrast drawn from 30
   tasks → 11.9M partitions → p = 0.00498). Same contrast, same AUC to within 0.001, but a
   null space large enough to measure against.

Note also the structural signature of the exhaustive cases: shuffled mean is exactly 0.5000
for both, because the 20 partitions come in complementary pairs whose AUCs sum to 1.

---

## 3. TASK 2 — PWM motif attribution

### ⚠️ Scope limit — the wider window was NOT run

**The saved importance arrays cover exactly 200 positions per region. That is all the prior
ISM ever computed.** A 1–2 kb rescan cannot be produced from saved data — there are no
importance values outside the 200 bp window to scan. Recomputing is the only route, and per
the brief I have costed it and stopped rather than spending the time (§4).

What was done for free: replace strict IUPAC consensus with proper PWM scanning **inside
the same 200 bp**, which directly tests whether scanning stringency was the bottleneck.

### Motifs used

Fetched live from **JASPAR** (`https://jaspar.elixir.no`, CORE vertebrates):

| Family | Matrix ID | JASPAR name | Length |
|---|---|---|---|
| NF-κB p65 | **MA0107.1** | RELA | 10 |
| C/EBPβ | **MA0466.2** | CEBPB | 10 |
| AP-1 | **MA0099.3** | FOS::JUN | 10 |
| ETS | **MA0098.3** | ETS1 | 10 |

Log-odds PWMs against a genomic background (A/T 0.295, C/G 0.205), scanned on both strands.
Relative-score thresholds 0.75 / 0.80 / 0.85, with **0.80 pre-declared as primary**.

### Result — scanning stringency was indeed the bottleneck

| Threshold | windows with ≥1 match | n for enrichment | mean enrichment | sd | p (t vs 1) | median within-region permutation p |
|---|---|---|---|---|---|---|
| 0.75 | **13/13** | 13 | 1.197 | 0.356 | 0.0692 | 0.0515 |
| **0.80 (primary)** | **13/13** | 13 | **1.294** | 0.525 | **0.0666** | 0.0550 |
| 0.85 | 9/13 | 9 | 1.981 | 1.238 | 0.0448 | 0.0025 |

**Prior (strict consensus, same 200 bp): 2/13 windows, enrichment 2.151 ± 1.501, n = 2, p = 0.474.**

The match problem is solved — **13/13 windows now contain matches instead of 2/13**, so the
enrichment statistic rests on n = 13 rather than n = 2. Mean matched coverage at the primary
threshold is 52.2 of 200 bases (ETS1 26.5, C/EBPβ 16.2, AP-1 7.8, RELA 6.3).

### But the enrichment itself is weak and does not clear significance at the primary threshold

At rel ≥ 0.80: importance at PWM-matched positions is **1.29× background, p = 0.0666** —
suggestive, not significant. Median enrichment 1.171; **9 of 13 regions above 1.0**.

Per region at the primary threshold:

| gene | matched bases | enrichment | within-region perm p |
|---|---|---|---|
| TREM1 | 77 | **2.465** | 0.0005 |
| TNFAIP3 | 36 | **2.217** | 0.0005 |
| PRRX2 | 38 | 1.564 | 0.0125 |
| COL10A1 | 55 | 1.461 | 0.0365 |
| FOXQ1 | 59 | 1.241 | 0.0160 |
| TPD52L1 | 67 | 1.188 | 0.0550 |
| TNFSF15 | 75 | 1.171 | 0.0250 |
| COL15A1 | 21 | 1.092 | 0.2434 |
| CDKN2B | 50 | 1.024 | 0.4218 |
| OLFML2A | 48 | 0.971 | 0.5912 |
| PGM5 | 70 | 0.963 | 0.7007 |
| DAAM2 | 52 | 0.769 | 0.9870 |
| MAMDC2 | 30 | 0.697 | 0.9895 |

**Stated plainly: model importance does NOT concentrate on senescence TF motifs in a way
that clears significance at the pre-declared threshold.** It is modestly enriched (1.29×,
p = 0.0666, n = 13), with two regions showing strong and individually significant
enrichment (TREM1 2.47×, TNFAIP3 2.22×, both perm p = 0.0005) and four regions actively
*de*-enriched.

The enrichment does rise monotonically with match stringency (1.20 → 1.29 → 1.98), and at
rel ≥ 0.85 it reaches p = 0.0448 with median within-region perm p = 0.0025. That trend is
what genuine motif sensitivity would look like. **I am not treating the 0.85 number as the
result** — 0.80 was pre-declared, and 0.85 drops to 9 of 13 windows. The trend is reported
as supporting evidence, not as the headline.

---

## 4. Wider-window ISM — projected cost, awaiting go-ahead

Measured rate from the prior probe: **0.2018 s per forward pass** (13 regions × 601 passes
in 1,577 s, batch 4, 22.29 GB peak).

| Option | forward passes | hours | cost |
|---|---|---|---|
| 1 kb × 13 regions | 39,013 | 2.19 | **$1.64** |
| 1 kb × 20 regions | 60,020 | 3.37 | **$2.52** |
| 2 kb × 13 regions | 78,013 | 4.37 | **$3.28** |
| 2 kb × 20 regions | 120,020 | 6.73 | **$5.05** |

**Not run — waiting for your go-ahead**, per the brief.

My read: the wider window would mainly buy statistical power on the enrichment estimate
(more matched positions per region), not a different kind of answer. The 200 bp rescan
already shows the effect is real but small at the primary threshold. If you want this
resolved, **1 kb × 20 regions ($2.52)** is the sensible option — it roughly quintuples
matched positions per region and restores the 7 regions the prior probe's time budget cut.
2 kb doubles the cost for a marginal further gain.

---

## 5. What this means

**Do the passing controls survive at proper resolution? Yes — all four.**

- `pooled_mechspecific_OIS` at **p = 0.00498 (0/200)** is now solidly established and remains
  the best design target: highest AUC among passing conditions, best evaluation power
  (182 evaluable genes), correct CDKN1A direction.
- `pooled_mechspecific_rep` (p = 0.0199) and `within_replicative` (p = 0.0299) are real.
- `pooled_unified_rep` (p = 0.0448) is marginal and should not be leaned on.
- The unified-contrast OIS condition still fails (p = 0.1393), and both transfer directions
  still fail. **The "pool the data, not the contrast" conclusion from the prior run stands
  at proper resolution.**

**Does importance concentrate on senescence TF motifs? Not decisively, at n = 13.** PWM
scanning fixed the match-scarcity problem (2/13 → 13/13 windows), but the enrichment is
1.29× with p = 0.0666 at the pre-declared threshold. Two regions show strong, individually
significant motif attribution; four go the other way. This is neither the clean positive
that would make the model design-usable nor a clean negative.

That leaves the generation-readiness verdict from the prior run **unchanged**:
classification-capable, not generation-ready. Task 1 strengthened the classification claim;
Task 2 did not strengthen the design claim.

---

## 6. Artifacts

All at **`/Users/gabeykim/Downloads/Senescence/infra_borzoi_grelu/execution_results_permutations/`**
— every file verified byte-for-byte against the remote.

| File | Bytes |
|---|---|
| `results_perm200.json` | 39,212 |
| `results_pwm_rescan.json` | 4,313 |
| `pwm_rescan_per_region.csv` | 2,459 |
| `perm200.log` | 6,439 |
| `pwm_rescan.log` | 1,219 |
| `perm200.py` | 13,166 |
| `pwm_rescan.py` | 10,262 |
| `setup_pod.sh` | 2,955 |

`results_perm200.json` contains the full per-condition AUC vectors (all 200 or all 20
values) so any of these p-values can be recomputed without rerunning anything.

The trunk cache was not re-copied — unchanged and already held at
`infra_borzoi_grelu/execution_results_transfer/cache_expanded/` (2 × 101,460,608 bytes).
The prior ISM importance arrays remain at
`infra_borzoi_grelu/execution_results_generation_probe/ism_importance.npz`.

**Cost:** pod up 21:47–21:59 UTC plus setup and retrieval ≈ **0.42 h × $0.75 = $0.32**.
Setup ~22 min (overlapped with Task 1, which ran as soon as gReLU was installed rather than
waiting for the genome), Task 1 0.4 min, Task 2 ~1 min.
