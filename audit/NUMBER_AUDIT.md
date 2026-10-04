# Number audit — claim inventory vs saved artifacts

Audit date 2026-10-04. Local only, no GPU, nothing recomputed that was already computed.
Paths are relative to `~/Downloads/Senescence/`.

**Result: of 72 audited line items — 60 MATCH, 4 MISMATCH, 1 NOT FOUND, 3 existence-only
checks FOUND, and 4 that verify numerically but are misattributed or mislabelled.** One
accession (`GSE254358`) is confirmed non-existent against live GEO, reproducing the project's
own finding. All internal n-arithmetic is consistent.

Counts are per line item in `audit_table.csv`; several items bundle more than one figure (e.g.
B1 covers three values), so the number of individual figures checked is larger than 72.

Where a claim was verified by reading a stored float and rounding it, the stored float is given
in full so the rounding can be re-checked.

---

## Group A — within-replicative reproducibility

Primary source: `output/replicative_ceiling_test/replicative_correlations.csv`, rows at the
primary threshold 0.5 unless noted.

| # | Claimed | Source file | Stored value | Verdict |
|---|---|---|---|---|
| A1 | R2 = +0.054, n = 134,336 | `output/replicative_ceiling_test/replicative_correlations.csv` (R2, thr 0.5) | `0.05361793271687495`, n `134336` | **MATCH** |
| A2 | R1 = +0.144, n = 144,754 | same (R1, thr 0.5) | `0.14393701486681684`, n `144754` | **MATCH** |
| A3 | R3 = +0.268, n = 124,962 | same (R3, thr 0.5) | `0.268230379137125`, n `124962` | **MATCH** |
| A4 | Within-OIS ceiling = +0.653, n = 88,530 | `output/enhancer_test/cross_mechanism_correlations.csv` (b, thr 0.5) | `0.6527365495816069`, n `88530` | **MATCH** |
| A5 | Cross-mechanism range −0.164 to +0.164 | `replicative_correlations.csv` (X1–X4, thr 0.5) | min X3 `−0.1640636031375235`, max X4 `+0.1638694927961774` | **MATCH** |
| A6 | R1 threshold sensitivity +0.063 to +0.327 | same (R1, thr 0.0 and 1.0) | `0.06329813662889718` … `0.3267534819931631` | **MATCH** |
| A7 | OIS ceiling threshold sensitivity +0.640 to +0.678 | `cross_mechanism_correlations.csv` (b, thr 0.0 and 1.0) | `0.639847571824792` … `0.6778941781755103` | **MATCH** |
| A8 | X1 = +0.114 | `replicative_correlations.csv` (X1, thr 0.5) | `0.11423358278017126`, n `123121` | **MATCH** |
| A9 | R2 descriptor: "same lab/pipeline/aligner/hg18 build" | same, `note` column | "within-replicative, SAME lab/pipeline/build, **different cell line** (upper bound)" | **MISATTRIBUTED DESCRIPTOR** — see D2 below |

---

## Group B — locus checks

| # | Claimed | Source file | Stored value | Verdict |
|---|---|---|---|---|
| B1 | CDKN2A promoter, 3 replicative: −0.85, −0.63, −0.43 | `output/replicative_ceiling_test/locus_sanity_checks.csv` | `−0.8482555751480071`, `−0.6318056309802202`, `−0.4306104011673332` | **MATCH** |
| B2 | CDKN2A promoter, 2 OIS: −1.03, −1.07 | same | GM21 `−1.0291839495298947`, IMR90_prior `−1.0675939010597921` | **MATCH** (but see D5) |
| B3 | CDKN2A distal, replicative: −0.17, **+0.22**, +0.31 | same | `−0.16624652949671184`, **`+0.21488386363151507`**, `+0.3113982819661036` | **MISMATCH** on the middle value |
| B4 | CDKN2A distal, OIS: +0.24, +0.49 | same | `+0.23961104732871571`, `+0.4899655221683543` | **MATCH** |
| B5 | Model chr9 promoter −0.594 vs measured −1.341 | `output/ois_enhancer_run/results/locus_checks.csv` | predicted `−0.5936548057943583`, measured_IMR90 `−1.3411686420440674` | **MATCH** |
| B6 | Model chr9 distal +0.198 vs measured +0.421 | same | predicted `0.19800283095123722`, measured_IMR90 `0.42054158449172974` | **MATCH** |

---

## Group C — expression-phase cross-study correlations

| # | Claimed | Source file | Stored value | Verdict |
|---|---|---|---|---|
| C1 | OIS vs OIS = +0.531 | `output/trainset_expanded/VERIFY_LOG.txt:45` (GSE205692 vs GSE74324); also `REPORT_multistudy_trainset.md:28,361` | `+0.531` | **MATCH** |
| C2 | Replicative vs IMR90-OIS = +0.142 | `VERIFY_LOG.txt` (GSE175533 vs GSE74324) | `+0.142` | **MATCH** |
| C3 | Replicative vs GM21-OIS = +0.093 | `VERIFY_LOG.txt` (GSE175533 vs GSE205692) | `+0.093` | **MATCH** |
| C4 | Pooled mech-specific OIS AUC 0.7843, p 0.00498, 0/200 | `infra_borzoi_grelu/execution_results_permutations/REPORT_permutations.md:56`, `results_perm200.json` | `0.7843`, `0/200`, `0.00498` | **MATCH** |
| C5 | Pooled unified on OIS 0.6895, 27/200, p 0.139 | `REPORT_permutations.md:61`, `perm200.log:58` | `0.6895`, `27/200`, `0.13930` | **MATCH** |
| C6 | IR vs replicative, independent controls = +0.115 | `replicative_correlations.csv` (I3, thr 0.5) | `0.11527686999154328` | **MATCH** |
| C7 | IR vs OIS = **+0.084** | `replicative_correlations.csv` (I2, thr 0.5) | **`0.08349433131616606`** | **MISMATCH** |
| C8 | Shared-control: observed +0.564, constructed null +0.523 | `output/replicative_ceiling_test/shared_denominator_null.txt` | `+0.5642`, `+0.5231` | **MATCH** |

---

## Group D — OIS enhancer model

Primary source: `output/ois_enhancer_run/results/results_train.json`.

| # | Claimed | Stored value | Verdict |
|---|---|---|---|
| D1 | Spearman held-out IMR90 = +0.5484, 86.1% of +0.637 | `eval_rest.IMR90_SEN_vs_PRO.spearman 0.5484116314300708`; `frac_of_ceiling 0.8609287777552131`; `reproducibility_ceiling 0.637` | **MATCH** |
| D2 | AUC = 0.9070 | `0.907004556857948` | **MATCH** |
| D3 | Evaluable 26,936 of 100,070, ambiguous 73,134 | `26936` / `100070` / `73134` | **MATCH** |
| D4 | chr9+chr6: +0.5250 (82.4%), AUC 0.8899 | `0.5249507759070756`, `0.8240985493046713`, `0.8899005853011095` | **MATCH** |
| D5 | IMR90 SEN vs QUI: +0.5895, AUC 0.9354 | `0.5895035068695315`, `0.9354146775192059` | **MATCH** |
| D6 | Shuffled-label 0/200 on both, p = 0.00498 | `spearman.n_beating 0`, `auc.n_beating 0`, `empirical_p 0.004975124378109453` | **MATCH** |
| D7 | Untrained head −0.0903, AUC 0.3842 | `−0.09030808019676641`, `0.3842417077772162` | **MATCH** |
| D8 | Shuffled embeddings: +0.548 → +0.025 | `0.5484116314300708` → `0.024962857674289466` | **MATCH** |
| D9 | Prediction sd 0.5137 vs measured 0.9601 | `eval_rest.pred_sd 0.5137194991111755`, `measured_sd_IMR90 0.960066556930542` | **MATCH** (but see D4 below) |

The 200 stored permutation Spearmans were re-counted: 0 of 200 reach the real +0.5484, and the
stored maximum (`0.094796`) matches the list. Consistent.

---

## Group E — training set build

| # | Claimed | Source file | Stored value | Verdict |
|---|---|---|---|---|
| E1 | 112,099 consensus regions, 111,671 surviving | `output/ois_enhancer_trainset/region_stats.json`; `output/ois_enhancer_run/ois_cache/sanitize_stats.json` | `n_consensus_regions 112099`; `n_regions_cached 111671` | **MATCH** |
| E2 | Losses 247 past end, 154 before start, 21 liftover, 6 ambiguous; 0.38% | `output/ois_enhancer_trainset/STEP2_LOG.txt` §5 and §7 | `247`, `154`, `21`, `6`; 428/112,099 = 0.3818% | **MATCH** |
| E3 | Liftover 99.92% midpoint, 99.96% span | `output/enhancer_test/harmonization_stats.json` | `99.91758611410118`, `99.95996797437951` | **MATCH but MISATTRIBUTED** — see D1 |
| E4 | 111,671 regions → **8,459** windows, 13.2 per window | `output/ois_enhancer_run/ois_cache/sanitize_stats.json` | `n_regions_cached 111671`, `n_windows 8459` (ratio 13.201) | **MATCH** — but see **D3**, a second window count exists |
| E5 | Median K27ac→ATAC summit distance 245 bp | `output/ois_enhancer_trainset/STEP1_LOG.txt` §2 | `245 bp` | **MATCH** |
| E6 | Cross-dataset reproducibility +0.637 | `output/ois_enhancer_trainset/step2_stats.json` | `spearman_gm21_vs_imr90_primary 0.6365320223147795` (n = 83,716) | **MATCH** |
| E7a | 35% of diff-ATAC regions show no H3K27ac change | `output/enhancer_test/step2_stats.json` (thr 0.5) | `frac_also_diff 0.6497684621009018` → 35.02% | **MATCH** |
| E7b | **93.0%** directional agreement where both move | same `n_same_direction 12429` / `n_also_diff_k27ac 13330` = **93.24%**; `REPORT_enhancer_cross_mechanism.md:225` states **93.3%** (12,432) | **MISMATCH** — see D6 |
| E8 | Quiescence +0.880, +0.036; sd 0.46 vs 0.97 | `STEP2_LOG.txt` §1 | `+0.8800`, `+0.0355`; `sd 0.4598` vs `sd 0.9652` | **MATCH** |
| E9 | Batch PCA status R² 0.777 vs study 0.448 | `output/ois_enhancer_trainset/batch_pca_stats.csv` (raw pooled, PC1) | `r2_status 0.776726601064073`, `r2_study 0.44766651259202994` | **MATCH** |
| E10 | Identity worst deviation 3.31e-07 vs 1e-04 | `output/ois_enhancer_run/ois_cache/preflight.json` | `identity_worst 3.3053336317829245e-07`, `tolerance 0.0001` | **MATCH** |
| E11 | Bin arithmetic 5,387 pairs, 0 inconsistent | same | `bin_arithmetic_pairs_checked 5387`, `bin_arithmetic_bad 0` | **MATCH** |

---

## Group F — design gates

| # | Claimed | Source file | Stored value | Verdict |
|---|---|---|---|---|
| F1 | Required delta = +0.861 | `output/ois_enhancer_run/results/results_probe.json` | `delta_median_to_p95 0.8612759113311768` | **MATCH** |
| F2 | Best = +0.6657, **sd** 0.0487, n=10, NF-κB+ETS clustered K=12 | `output/ois_gates/results/gate1_results.json` → `best` | `adjusted_delta 0.6656626462936401`, `adjusted_sem 0.048732499060084866`, `K_total 12`, `clustered` | **MATCH on value; MISLABELLED** — see D7 |
| F3 | Shortfall +0.196, ratio 0.773× | same | `frac_of_required 0.7728796748359079`; 0.8612759 − 0.6656626 = `0.1956133` | **MATCH** |
| F4 | Best pair over best single +0.0639, **sd** 0.0414, p = 0.157 | `output/ois_gates/REPORT_ois_gates.md` §1.3 (prose only) | recomputed from `gate1_additivity_raw.csv`: mean `0.0639`, SEM `0.0414`, SD `0.1309`, p `0.1571` | **MATCH on value; NOT IN ANY RESULTS FILE; MISLABELLED** — see D7, D8 |
| F5 | **All five conditions peak at K=12, decline at K=18 and K=20** | `output/ois_gates/results/gate1b_dose_curve.csv` | ETS1 peaks at **K=6** (`+0.2277` vs `+0.2147` at K=12); triple **rises** K18→K20 (`+0.1455` → `+0.2753`) | **MISMATCH** |
| F6 | Clustered beat interleaved 8/8, significant in 3, up to +0.294 | `output/ois_gates/results/gate1_spacing.csv` | 8/8 negative differences; 3 with p<0.05; max abs difference `0.2939751744270324` | **MATCH** |
| F7 | NF-κB +0.665 at 8 copies, scrambled +0.125, ratio 5.3× | `results_probe.json` → `motif_insertion.summary` | `0.6653652131557464`, `0.12504179477691646`; ratio 5.321 | **MATCH** |
| F8 | C/EBPβ scrambled control = +0.068 | same | `0.06829091310501101` | **MATCH** |
| F9 | AP-1 delta −0.464, Spearman −0.375, scrambled −0.016 | same | `−0.46432673633098603`, `−0.37484156984116956`, `−0.015941154956817583` | **MATCH** |

---

## Group G — AP-1 resolution

| # | Claimed | Source file | Stored value | Verdict |
|---|---|---|---|---|
| G1 | Ablation high-response = −0.1055, p = 0.032 | `output/ois_gates/results/gate2_ap1_results.json` → `H1b_interaction.AP1_FOSJUN.by_stratum.high` | `specific −0.10550387452046073`, `p 0.03221350924508854` | **MATCH** |
| G2 | Ablation low-response = +0.0658 | same, `.low` | `0.06576001347275451` | **MATCH** |
| G3 | Interaction p = 0.0091 | same, `.interaction.welch_p` | `0.009063013838304235` (Mann-Whitney `0.00757959203158261`) | **MATCH** |
| G4 | Continuous Spearman −0.364, p = 0.011, n = 48 | same, `.continuous_vs_gm21` | `−0.36354754667824574`, `0.011088281505031718`, `48` | **MATCH** |
| G5 | NF-κB positive control −0.1404, p = 0.0003 | same, `NFKB_RELA.by_stratum.high` | `−0.1404095912973086`, `0.00029890902470616266` | **MATCH** |
| G6 | AP-1 ISM attribution 3.61×, p = 0.0002, n = 16 | `output/ois_gates/results/gate2_h1a_ism_enrichment.csv` | mean `3.6086`, p `0.000216`, n `16` | **MATCH** |
| G7 | AP-1 frequency across deciles 0.086 → 0.179, OR 2.31, p = 5.9e-94 | `output/ois_gates/results/gate2_h3_frequency_by_decile.csv` | `0.086139` → `0.178562`, OR `2.306186`, p `5.912336e-94` | **MATCH** (see D9 on "to 0.179") |
| G8 | AP-1 regions H3K27ac +0.206, p = 2e-283, n = 16,909 | `output/ois_gates/results/gate2_h3_motif_vs_assay.csv` | `delta_k27 0.206123`, p `1.999463e-283`, `n_with 16909` | **MATCH** |
| G9 | Top ATAC quartile, AP-1 still higher H3K27ac, p = 4.9e-05 | `output/ois_gates/results/gate2_h3_conditional.csv` | `delta_k27 +0.048217`, p `4.900534e-05`, n_with `4657` | **MATCH** |
| G10 | Sign agreement 50/50 and 20/20 | `results_probe.json` → `dynamic_range` | `sign_agreement_top50 50`, `sign_agreement_top20 20` | **MATCH** |

---

## Group H — provenance and infrastructure

| # | Claimed | Source file | Stored value | Verdict |
|---|---|---|---|---|
| H1 | Borzoi 185,892,699 parameters | `infra_borzoi_grelu/REPORT.md:155`; `execution_results/task1_output.log:3` | `185,892,699` | **MATCH** |
| H2 | 524,288 bp input | `infra_borzoi_grelu/REPORT.md:156` | `524,288 bp` | **MATCH** |
| H3 | 32 bp bins | `infra_borzoi_grelu/REPORT.md:157` | `32 bp bins` | **MATCH** |
| H4 | 7,611 tracks | `REPORT.md:158`; `execution_results/task1_output_tracks.csv` (7,611 rows) | `7,611` (CHIP 3,886 + RNA 1,543 + CAGE 1,276 + DNASE 674 + ATAC 232 = 7,611) | **MATCH** |
| H5 | hg38 | `REPORT.md:159` | `hg38` (model's declared metadata) | **MATCH** |
| H6 | IUPAC counts, 19,652-gene run | `infra_borzoi_grelu/execution_results_finetune/sanitize_stats.json` | `genes_sanitized 49`, `chars_replaced 133` (Y49 W17 K6 R33 M19 B4 S5); `REPORT.md:1102` "49 of 19,652" | **FOUND** |
| H7 | IUPAC counts, 13,211-gene run | `infra_borzoi_grelu/execution_results_transfer/cache_expanded/sanitize_stats.json` | `genes_sanitized 31`, `chars_replaced 70`; `total_bases_scanned 6,926,368,768` = 13,211 × 524,288 exactly | **FOUND** |
| H8 | IUPAC counts, enhancer run | `output/ois_enhancer_run/ois_cache/sanitize_stats.json` | `windows_sanitized 33`, `chars_replaced 98` (R24 Y41 B4 M10 K6 S4 W9) | **FOUND** |
| H9 | **Total compute spend across all sessions** | — | no artifact contains a project-wide total | **NOT FOUND** — see D10 |

---

# Mismatches, not-founds and flags

## D1 — E3 liftover percentages are from the wrong build *(misattribution)*

The claim "99.92% midpoint success, 99.96% span preservation" is filed under *training set
build*. Those numbers are real but belong to the **cross-mechanism enhancer test**
(237,628/237,824 midpoints, `output/enhancer_test/harmonization_stats.json`, and
`REPORT_enhancer_cross_mechanism.md:83–85`).

The OIS training set's own liftover is **112,078 / 112,099 = 99.98%, 21 failed**
(`output/ois_enhancer_trainset/STEP1_LOG.txt` §3, `REPORT_ois_enhancer_trainset.md:42`,
`region_stats.json: liftover_pct 99.98126655902372`).

**Assessment: misattribution.** Both figures are correct for their own build; pairing the
237,824-anchor liftover with the 112,099-region training set in a Methods section would
misdescribe the training set. Note also that "99.96% span preservation" is *conditional* —
4,994 / 4,996 of the anchors whose **both ends lifted**, not 99.96% of all anchors.

## D2 — A9 R2 descriptor drops a material qualifier *(descriptive error)*

Claim: "same lab/pipeline/aligner/hg18 build". Stored note: "within-replicative, SAME
lab/pipeline/build, **different cell line** (upper bound)". R2 is IMR90 vs **BJ** — two
different cell lines. "aligner" and "hg18" do not appear in the note (hg18 is corroborated
separately by `quantification_stats.json`, which records the hg19→hg18 liftover).

**Assessment: transcription error.** Omitting "different cell line" changes what the +0.054
bounds: as written it reads like a technical-replicate floor, when it is a cross-cell-line
comparison.

## D3 — E4: two different window counts exist *(version drift)*

| Artifact | Regions | Windows | Ratio |
|---|---|---|---|
| `ois_enhancer_trainset/step2_stats.json`, `STEP2_LOG.txt` §6, `REPORT_ois_enhancer_trainset.md:14,279,368`, `window_covering_per_chrom.csv` (sums exactly) | **111,677** | **8,432** | 13.244 |
| `ois_enhancer_run/ois_cache/sanitize_stats.json` | **111,671** | **8,459** | 13.201 |

The covering was computed *before* the 6 segmental-duplication regions were dropped
(`STEP2_LOG.txt` §7), so the executed cache used a re-derived covering. The claim pairs
111,671 with 8,459 — internally consistent and matching the cache that was actually run.

**Assessment: version drift, claim takes the correct (executed) version.** The risk is the
other direction: the trainset report says 8,432 in four places. If the Methods cite the report,
they will disagree with the claim. Decide which number the manuscript uses and make it one.

## D4 — D9: prediction sd exists at two values *(different denominators, easily confused)*

| Artifact | pred sd | measured sd | n |
|---|---|---|---|
| `results_train.json` → `eval_rest` | `0.5137194991111755` | `0.960066556930542` | 100,070 |
| `results_probe.json` → `dynamic_range` | `0.5138489007949829` | `0.9613737463951111` | 111,671 |

The claim uses the `results_train.json` pair (0.5137 / 0.9601). Both are correct for their own
region set; they are not alternative estimates of the same quantity.

**Assessment: not an error, but state n.** "Prediction sd 0.5137 against measured 0.9601" is
only unambiguous if the preprint says it is the 100,070 non-test regions.

## D5 — B2 omits a third OIS dataset *(selective reporting, not a number error)*

`locus_sanity_checks.csv` carries **three** OIS H3K27ac columns at the CDKN2A promoter:
GM21 `−1.029`, IMR90-prior `−1.068`, and **IMR90-OIS-same-lab `−0.0548`**. The claim cites the
first two as "two OIS datasets". The third does not show promoter-down and is n=9 regions.

**Assessment: the two quoted values are correct.** Flagged because "CDKN2A promoter goes down
in OIS" is stated more strongly by two of three than by three of three.

## D6 — E7b directional agreement: three different values *(mismatch + internal inconsistency)*

| Source | Same-direction count | Percentage of co-differential |
|---|---|---|
| Claim | — | **93.0%** |
| `output/enhancer_test/step2_stats.json` | `n_same_direction 12429` | 12,429/13,330 = **93.24%** |
| `output/enhancer_test/REPORT_enhancer_cross_mechanism.md:225` | **12,432** | **93.3%** |

Two problems. (a) The claim's 93.0% matches neither source. (b) The report's count (12,432)
disagrees with the machine-readable JSON (12,429) by 3 regions — the report's own derived
percentage, 93.3%, is computed from the larger count.

**Assessment: transcription error compounded by an unexplained count discrepancy.** 93.0% looks
like a downward round of 93.2–93.3%. The 12,429 vs 12,432 gap is the more serious finding: one
of the two was produced by a different pass and the discrepancy is not documented anywhere.
`step2_stats.json` is the machine-written artifact and should be treated as authoritative unless
the report's provenance can be re-established.

## D7 — F2, F4: "sd" is actually SEM *(mislabelling, factor of 3.16)*

| Claimed as "sd" | Stored field | True SD (n=10) |
|---|---|---|
| 0.0487 (best condition) | `gate1_results.json: best.adjusted_sem` = `0.048732499060084866` | **0.1541** |
| 0.0414 (pair over single) | SEM of the paired difference | **0.1309** |

The source files and `REPORT_ois_gates.md` both label these **SEM** (`±SEM` column header,
`adjusted_sem` field name). The claim inventory relabels them "sd".

**Assessment: transcription error, and the one most likely to mislead a referee.** Reporting
+0.6657 ± 0.0487 as a standard deviation understates the spread of the 10 background seeds by
3.16×. Either keep SEM and say so, or quote SD = 0.154.

## D8 — F4 has no machine-readable source *(never written to a results file)*

The pair-over-single contrast (+0.0639, p = 0.157) appears **only in the prose of
`output/ois_gates/REPORT_ois_gates.md` §1.3**. No CSV or JSON in `output/ois_gates/results/`
stores it; `gate1_additivity_summary.csv` has the per-condition values but not this contrast.

It **is** exactly reproducible from `gate1_additivity_raw.csv` — this audit recomputed it and
got mean 0.0639, SEM 0.0414, SD 0.1309, paired t = +1.543, p = 0.1571, n = 10, matching the
claim. Flagged because it is the single most load-bearing negative result in Gate 1 ("combining
motif families buys nothing") and it currently has no stored artifact.

**Assessment: never computed into a file, though correct.** Worth emitting to a CSV before
submission so it has provenance like every other number.

## D9 — G7 "0.086 to 0.179" is first-to-last decile, not min-to-max *(precision nuance)*

Both endpoints are correct (`0.086139`, `0.178562`), but the AP-1 frequency series is
**non-monotonic** and peaks at decile 6 (`0.1989`), then falls across deciles 7–9. "0.086 to
0.179" invites the reading that 0.179 is the maximum. Full series:

`0.086, 0.093, 0.109, 0.128, 0.155, 0.173, 0.199, 0.197, 0.196, 0.179`

**Assessment: not a number error; a wording risk.** "Rises monotonically" would be wrong; "top
decile 2.07× the bottom (OR 2.31)" is right and is what the OR encodes.

## D10 — H9 total compute spend: NOT FOUND

No artifact records a project-wide total. The only aggregate is
`infra_borzoi_grelu/REPORT.md:1181`: **"Cumulative approximate spend across all four pod
sessions: ~$2.72"** (~$0.75 + ~$0.55 + ~$0.91 + ~$0.51) — explicitly limited to the first four
pods and explicitly "derived from SSH command timestamps, not pulled invoices".

Session costs that exist individually and are **not** in that total:

| Session | Stated cost | Source |
|---|---|---|
| Pods 1–4 (expression phase → calibration) | ~$2.72 | `infra_borzoi_grelu/REPORT.md:1181` |
| Pod 5 — real fine-tune | **not recorded** | REPORT.md:1170 says pod 5 was still running at write time |
| Transfer / expanded cache | ~$1.72 (2.3 h) | `execution_results_transfer/REPORT_transfer.md:254` |
| Generation probe | $0.58 (0.77 h) | `execution_results_generation_probe/REPORT_generation_probe.md:266` |
| Permutations (200) | $0.32 (0.42 h) | `execution_results_permutations/REPORT_permutations.md:256` |
| OIS enhancer run | $2.65 (3.53 h) | `output/ois_enhancer_run/REPORT_ois_enhancer_run.md:4` |
| OIS gates | $0.64 (0.857 h) | `output/ois_gates/REPORT_ois_gates.md` §Cost |
| Trainset builds, ceiling test, enhancer test | $0.00 (local) | `REPORT.md:1159–1161` |

Documented sessions sum to **~$8.63**, with pod 5 unknown. That sum is this audit's arithmetic,
not a stored figure. **Do not quote a total until pod 5 is priced and the figure is written to
an artifact.**

## D11 — B3: CDKN2A distal BJ value *(double-rounding)*

Stored `0.21488386363151507`. Claim says **+0.22**. Correct rounding is +0.215 (3 dp) or +0.21
(2 dp). The source report is itself inconsistent: `REPORT_replicative_ceiling.md:252` gives
**+0.215** in the table, while line 268 writes **+0.22** in the prose. The claim inherited the
prose value, which is +0.215 rounded a second time.

**Assessment: transcription error via double-rounding.** Small, but it is the number used to
argue the replicative distal signal is "inconsistent in sign", so it will be read closely.

## D12 — C7: IR vs OIS *(rounding error, inherited from the report)*

Stored `0.08349433131616606` → **+0.083**. Claim says **+0.084**, and so does
`REPORT_replicative_ceiling.md:234` and `:301`.

**Assessment: transcription error originating in the source report, not in the claim
inventory.** Correcting it means correcting the report too. Note the paired framing "+0.115 vs
+0.084" becomes "+0.115 vs +0.083".

## D13 — F5: "all five conditions peak at K=12" is false *(mismatch)*

From `output/ois_gates/results/gate1b_dose_curve.csv`:

| Condition | K=6 | K=12 | K=18 | K=20 | Peak |
|---|---|---|---|---|---|
| NF-κB+ETS (clustered) | +0.3919 | **+0.6657** | +0.5444 | +0.3620 | K=12 |
| NF-κB (single) | +0.4196 | **+0.6018** | +0.5686 | +0.3884 | K=12 |
| C/EBPβ (single) | +0.1877 | **+0.5837** | +0.3585 | +0.2576 | K=12 |
| NF-κB+C/EBPβ+ETS (clustered) | +0.3033 | **+0.3889** | +0.1455 | **+0.2753** | K=12, but **rises** K18→K20 |
| **ETS (single)** | **+0.2277** | +0.2147 | +0.1541 | +0.1022 | **K=6** |

Two defects: ETS peaks at K=6, not K=12; and the triple does not decline monotonically — it
rises from +0.1455 at K=18 to +0.2753 at K=20.

**Assessment: transcription error, and it originates in `REPORT_ois_gates.md` itself**, whose
prose says "Every condition peaks at K = 12 and declines beyond it" while its own table
(§1.4) bolds +0.2277 at K=6 for ETS. The report is internally inconsistent.

Accurate statement: **four of five conditions peak at K=12; ETS peaks at K=6; all five are
lower at K=20 than at their peak.** The headline conclusion (saturation with an interior
optimum, no headroom above 0.773×) is unaffected.

---

# Accession verification (re-run live, 2026-10-04)

Queried `eutils.ncbi.nlm.nih.gov` directly (`esearch db=gds …[ACCN] AND gse[ETYP]`, then
`esummary`), and `encodeproject.org` JSON — not read from `output/tables/accession_verification_log.txt`.

| Accession | Live status | Title / type | BioProject |
|---|---|---|---|
| GSE175533 | **LIVE** | Revisiting the Hayflick Limit… (n=147) | PRJNA732700 |
| GSE74324 | **LIVE** | BRD4 connects enhancer remodeling… (RNA-seq, n=72) | PRJNA299683 |
| GSE74238 | **LIVE** | BRD4 connects enhancer remodeling… (ChIP-Seq, n=22) | PRJNA299593 |
| GSE206402 | **LIVE** | Escape From OIS… [ATAC-seq] (n=25) | PRJNA850441 |
| GSE106146 | **LIVE** | Histone acetyltransferase p300 induces de novo super-enhancers… (n=82) | PRJNA415730 |
| GSE205692 | **LIVE** | Escape From OIS… (n=34) | PRJNA847199 |
| GSE210285 | **LIVE** | Integrated multi-omics… senescence landscape (n=4) | PRJNA865031 |
| GSE205898 | **LIVE** | Escape From OIS… [ChIP-seq] (n=30) | PRJNA848201 |
| GSE220545 | **LIVE** | Transcriptome of senescence-committed vs uncommitted… (n=6) | PRJNA910178 |
| GSE146585 | **LIVE** | Functional genomics assays… (n=223) | PRJNA610876 |
| GSE99028 | **LIVE** | Expression profiling of etoposide-induced senescent cells (n=6) | PRJNA387040 |
| GSE206493 | **LIVE** | Escape From OIS… [RNA-seq] (n=12) | PRJNA851072 |
| GSE206496 | **LIVE** | Escape From OIS… (n=101) | PRJNA851070 |
| **GSE254358** | **DOES NOT EXIST** | `esearch` count = 0 | — |
| ENCSR200OML | **LIVE** | Experiment, ATAC-seq, released | — |
| ENCSR978WIX | **LIVE** | Annotation (ChromBPNet-model), released | — |
| ENCSR169RNU | **LIVE** | Annotation, released | — |

BioProjects re-verified live, all **LIVE**: PRJNA732700, PRJNA299683, PRJNA299593, PRJNA850441,
PRJNA847199, PRJNA865031, PRJNA910178, PRJNA851070, PRJNA851072, PRJNA848201, PRJNA439280,
PRJNA387040.

**GSE254358.** Confirmed non-existent, reproducing the project's own finding. `REPORT.md:11`
already states "GSE254358 does not exist", and lines 16–23 record the same `count: 0` plus an
HTTP 404 on the FTP path and six near-miss accessions that were checked. Two residual exposures:

1. `data/GSE254358/` still exists and contains three files that look like data but are not:
   `GSE254358_family.soft.gz` and `GSE254358_series_matrix.txt.gz` are 990 bytes each and
   decompress to `mailto:webadmin@ncbi.nlm.nih.gov` — NCBI error payloads, not records.
2. The accession is still referenced in 10+ files including `REPORT.md`,
   `REPORT_GSE175533_MULTI.md`, `output/encode_imr90_section.md`, `output/gse206402_section.md`,
   `output/tables/design_adequacy_log.txt`, and five scripts under `scripts/`.

Those references are to the *audit of* the non-existent accession, which is legitimate. But
nothing in the tree marks the `data/GSE254358/` directory as holding error payloads, and a
reader or a script globbing `data/GSE*` would treat it as a dataset.

**`output/tables/accession_verification_log.txt` does not cover** GSE106146, GSE205692,
GSE220545, GSE146585, GSE99028 or GSE254358 — six accessions that are used in the work. This
audit verified them directly; all are live except GSE254358.

---

# Internal n-consistency

Every AUC report in `results_train.json` was re-checked programmatically:

- `n_evaluable + n_ambiguous == n_regions_total` — **8/8 pass** (3 tasks × 2 region sets, plus
  control2 and control3).
- `n_positive + n_negative == n_evaluable` — **8/8 pass**.
- `train 88,282 + early_stop 11,788 == non_test 100,070` — **pass**.
- `non_test 100,070 + test_chr9_chr6 11,601 == 111,671` — **pass**, matches the h5ad shape and
  `n_regions_cached`.
- Permutation list length == 200, and recounting `n_beating` from the 200 stored Spearmans
  gives 0, matching the stored `n_beating` and the stored max (0.094796) — **pass**.
- Region losses: `247 + 154 + 21 + 6 = 428 = 112,099 − 111,671`, and 428/112,099 = 0.3818% —
  **pass**, the claimed 0.38% is correct.

**No claimed n is internally inconsistent.**

---

# Claimed precision exceeding the source

| Claim | Issue |
|---|---|
| F2 "sd 0.0487", F4 "sd 0.0414" | Not precision but **wrong statistic** — these are SEM (D7) |
| E7b "93.0%" | Source supports 93.24% (JSON) or 93.3% (report); 93.0% is supported by neither (D6) |
| B3 "+0.22" | Source `0.21488` supports +0.215 or +0.21 (D11) |
| C7 "+0.084" | Source `0.08349` supports +0.083 (D12) |
| G7 "0.086 to 0.179" | Endpoints exact; "to" implies a maximum the series does not have (D9) |
| D1 "86.1%" | `frac_of_ceiling 0.8609287777552131` → 86.09%; quoting 86.1% is fine, but it is a ratio of a 3-dp ceiling (0.637) to a 4-dp Spearman, so the third significant figure is not meaningful |
| E6 "+0.637" | Stored `0.6365320223147795`; correct to 3 dp, but it is the thr-0.5 value (n = 83,716) — the same quantity is +0.620/+0.624/+0.657 at the other three thresholds, so the threshold must be stated |

---

# Files consulted

`output/replicative_ceiling_test/`: `replicative_correlations.csv`, `locus_sanity_checks.csv`,
`batch_pca_stats.csv`, `step2_stats.json`, `quantification_stats.json`,
`shared_denominator_null.txt`, `REPORT_replicative_ceiling.md`.
`output/enhancer_test/`: `cross_mechanism_correlations.csv`, `harmonization_stats.json`,
`step2_stats.json`, `locus_sanity_checks.csv`, `REPORT_enhancer_cross_mechanism.md`.
`output/ois_enhancer_trainset/`: `region_stats.json`, `step2_stats.json`, `batch_pca_stats.csv`,
`locus_sanity_checks.csv`, `window_covering_per_chrom.csv`, `STEP1_LOG.txt`, `STEP2_LOG.txt`,
`REPORT_ois_enhancer_trainset.md`.
`output/ois_enhancer_run/`: `results/results_train.json`, `results/results_probe.json`,
`results/locus_checks.csv`, `ois_cache/preflight.json`, `ois_cache/sanitize_stats.json`,
`REPORT_ois_enhancer_run.md`.
`output/ois_gates/results/`: `gate1_results.json`, `gate1_additivity_raw.csv`,
`gate1_additivity_summary.csv`, `gate1_spacing.csv`, `gate1b_dose_curve.csv`,
`gate2_ap1_results.json`, `gate2_h1a_ism_enrichment.csv`, `gate2_h3_motif_vs_assay.csv`,
`gate2_h3_conditional.csv`, `gate2_h3_frequency_by_decile.csv`; `REPORT_ois_gates.md`.
`output/trainset_expanded/`: `VERIFY_LOG.txt`, `REPORT_multistudy_trainset.md`.
`output/tables/accession_verification_log.txt`.
`infra_borzoi_grelu/`: `REPORT.md`, `execution_results/task1_output.log`,
`execution_results_finetune/sanitize_stats.json`,
`execution_results_transfer/cache_expanded/sanitize_stats.json`,
`execution_results_transfer/REPORT_transfer.md`,
`execution_results_generation_probe/REPORT_generation_probe.md`,
`execution_results_permutations/REPORT_permutations.md`, `.../perm200.log`,
`.../results_perm200.json`.
`REPORT.md` (top level), `data/GSE254358/`.
