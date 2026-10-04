# ERRATA for the working reports in this directory

**The reports in `reports/` are a working record of how the analysis actually
proceeded. They are superseded by the manuscript. They have not been rewritten,
because the record of what was concluded when is itself worth keeping — but
every discrepancy found between a report and the stored data is listed here.**

Five discrepancies are known. All five were found by auditing the reports
against the machine-written artifacts in `results/`; none changes a conclusion,
and all five are corrected in the manuscript. Where a report and a stored file
disagree, **the stored file in `results/` is authoritative** — it was written by
the script that computed the number, whereas the report was typed afterwards.

How these were found: `audit/NUMBER_AUDIT.md` checked 72 claimed numbers against
their source files (60 matched); `audit_fixes/AUDIT_FIXES.md` then closed the
open items. The fifth was found later still, while building Figure 2.

---

## 1. Directional agreement between accessibility and H3K27ac

| | |
|---|---|
| **Report says** | `REPORT_enhancer_cross_mechanism.md:225` — "12,432 (60.6% of all; **93.3%** of the co-differential)" |
| **Stored value** | `results/enhancer_cross_mechanism/step2_stats.json` → `n_same_direction = 12429`, i.e. **93.24%** of the 13,330 co-differential anchors |
| **Correct** | **12,429 / 13,330 = 93.24%** |
| **Manuscript** | Figure S1 panel B and its caption |

Resolved by recomputing from the stored per-anchor table: the recomputation
returns 12,429, matching the JSON exactly. It is **not** a tie or threshold
difference — there are zero exact-zero signs, and `>` versus `>=` changes
nothing. 12,432 is not reproducible by any variant tested and the generating
script never computes 93.3%. Full workings:
`audit_fixes/provenance_directional_agreement.md`.

## 2. IR versus OIS correlation

| | |
|---|---|
| **Report says** | `REPORT_replicative_ceiling.md:234` and `:301` — "+0.084" |
| **Stored value** | `results/replicative_ceiling/replicative_correlations.csv`, comparison I2 at threshold 0.5 → `0.08349433131616606` |
| **Correct** | **+0.083** |
| **Manuscript** | Figure S3 and the shared-control section |

The paired framing "+0.115 vs +0.084" becomes "+0.115 vs +0.083".

## 3. CDKN2A distal signal in BJ replicative

| | |
|---|---|
| **Report says** | `REPORT_replicative_ceiling.md:268` — "−0.17, **+0.22**, +0.31" |
| **Stored value** | `results/replicative_ceiling/locus_sanity_checks.csv` → `0.21488386363151507` |
| **Correct** | **+0.215** (3 dp) or **+0.21** (2 dp) |
| **Manuscript** | Figure 4, which plots the stored value |

A double-rounding artifact: the same report's own table at line 252 correctly
prints +0.215, and the prose at line 268 rounds that a second time to +0.22.

## 4. Where the motif dose response peaks

| | |
|---|---|
| **Report says** | `REPORT_ois_gates.md` §1.4 — "Every condition peaks at K = 12 and declines beyond it" |
| **Stored value** | `results/design_gates/gate1b_dose_curve.csv` — ETS peaks at **K = 6** (+0.2277, versus +0.2147 at K = 12); the three-motif combination also *rises* from K = 18 (+0.1455) to K = 20 (+0.2753) |
| **Correct** | **Four of five conditions peak at K = 12. ETS peaks at K = 6. All five are below their own peak at K = 20.** |
| **Manuscript** | Figure 6A, which plots and annotates the stored values |

The report is internally inconsistent: its own table in the same section bolds
+0.2277 at K = 6 for ETS. The conclusion — saturation with an interior optimum,
best achievable 0.773× of the requirement — is unaffected.

## 5. Whether all within-replicative values fall inside the cross-mechanism range

| | |
|---|---|
| **Report says** | `REPORT_replicative_ceiling.md:24` — "**The three within-replicative values (+0.054, +0.144, +0.268) fall inside the range of the four cross-mechanism values (−0.164 to +0.164).**" |
| **Stored value** | `results/replicative_ceiling/replicative_correlations.csv` at threshold 0.5: R2 = +0.0536 (inside), R1 = +0.1439 (inside), **R3 = +0.2682 (outside; the upper bound is +0.1639)** |
| **Correct** | **Two of three, not three of three.** |
| **Manuscript** | Figure 2 and its caption |

Found last, while building Figure 2, and missed by the earlier number audit
because that audit verified R3 = +0.268 and the span −0.164 to +0.164 *separately*
and never compared them to each other. It is a check on a relational claim, not
on a value.

The argument is unaffected: no within-replicative comparison approaches the
within-OIS ceiling of +0.653, and the cleanest comparison — R2, same lab, same
pipeline, same genome build, no technical confound — is +0.054. Accurate
wording: *two of the three within-replicative comparisons (+0.054, +0.144) fall
inside the cross-mechanism range; the third (+0.268) exceeds it only slightly
and remains far below the within-OIS ceiling of +0.653.*

---

## Redactions

`reports/REDACTIONS.txt` lists 18 lines across 4 reports where an ephemeral
rented-GPU host address (`host:port`) was replaced with
`<pod-address-redacted>`. No scientific content was altered. Those instances are
the only modifications made to any report in this directory.

## Not an erratum, but worth stating

`REPORT_borzoi_grelu_infrastructure.md:1181` gives "cumulative approximate spend
across all four pod sessions: ~$2.72". That is correct **for the first four pods
only** and is not a project total; it predates five later GPU sessions. The
full accounting is in `audit_fixes/compute_cost_table.csv` (≥ $10.53 across ten
sessions, of which one is a reconstructed lower bound). Do not quote the $2.72
figure as a project total.
