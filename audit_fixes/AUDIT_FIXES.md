# Audit fixes — eight items closed

Date 2026-10-04. Local only, no GPU. No model, cache or training result was recomputed; the one
derived quantity (item 2) comes from a stored per-sequence CSV. Paths are relative to
`~/Downloads/Senescence/`.

**Two of the eight changed the finding rather than just documenting it:**

- **Item 4** — the number audit's explanation for the 8,432 vs 8,459 window gap was **wrong**.
  It is not the 6 dropped regions. Both numbers reproduce exactly from the *same* 111,671
  regions; the entire difference is the covering criterion.
- **Item 3** — the 12,429 vs 12,432 discrepancy is **not** a filtering or tie-handling
  difference. 12,429 is reproducible from the stored per-anchor table; **12,432 is not
  reproducible by any variant tested** and has no computational basis.

One correction to the audit's own wording: the GSE254358 files do not "decompress to" anything
(item 7) — they are not compressed.

---

## ITEM 1 — SEM versus SD

Swept every report, CSV column and figure in the project for dispersion values: 20 locations,
in `dispersion_audit.csv`. SEM→SD conversions for every stored gate value are in
`dispersion_sem_to_sd.csv`.

**Verdict: 7 need a fix, 3 need clarifying, 10 are already correct.** The project is in better
shape than the two caught instances suggested — every CSV column is explicitly named `*_sem`,
and every permutation-null dispersion in every report is correctly labelled "sd".

### The only mislabelled dispersion in any report

| File | Line | Stated | Actual | Correct |
|---|---|---|---|---|
| `output/ois_gates/REPORT_ois_gates.md` | 86 (§1.3) | `± 0.0414` (unlabelled) | SEM over 10 seeds | **SEM 0.0414; SD 0.1309** |

Everything else flagged "YES" is in the **claim inventory / manuscript draft**, not in a report:

| Claim | Stated | Actual | Correct |
|---|---|---|---|
| Gate 1 best condition, +0.6657 | "sd 0.0487" | SEM (`adjusted_sem`) | **SEM 0.0487; SD 0.1541** |
| Pair over single, +0.0639 | "sd 0.0414" | SEM | **SEM 0.0414; SD 0.1309** |

### Figures: error bars are SEM and say so nowhere

All five error-bar usages in the project are in `output/ois_gates/scripts/make_figures.py`, all
are SEM, and **no figure states it**. Four figures need a caption line:

| Figure | Error bar | n |
|---|---|---|
| `fig1_gate1_additivity.png` | `xerr = adjusted_sem` | 10 seeds |
| `fig2_gate1b_dose_ceiling.png` | `yerr = adjusted_sem` | 10 seeds |
| `fig3_h1b_ablation.png` | `yerr` = SEM of the per-region specific effect | 24 regions per stratum |
| `fig4_h1c_h2.png` | `xerr`/`yerr` = `sem` columns | 10 seeds |

`fig5_h3_atac_vs_k27ac.png` has no error bars (point estimates over all 111,671 regions).

### Already correct — do not change

| File | What | Why it is right |
|---|---|---|
| `REPORT_ois_enhancer_run.md:158-161` | "Shuffled mean ± sd" −0.0445 ± 0.0607 (Spearman), 0.4673 ± 0.0428 (AUC) | genuine SD of the 200-permutation null; fields `control1_shuffled_labels.{spearman,auc}.sd`; **n = 200** |
| `infra_borzoi_grelu/REPORT.md:43, 1012` | "mean ± sd (20 permutations)" 0.5605 ± 0.0940 | field `control1_shuffled_labels.std`; **n = 20** |
| `REPORT_generation_probe.md:43-52` | "shuffled mean ± sd" | SD of the 20-permutation null; **n = 20** |
| `REPORT_permutations.md:54-61` | "shuffled mean ± sd" | SD of the 200-permutation null; **n = 200** |
| all gate CSVs | `*_sem` / `sem` columns | field names state the quantity; `n_seeds` column carries n |
| `REPORT_ois_gates.md:65,175,194,213` | `±SEM` column headers | explicit |

### Clarify (value correct, label implicit)

| File | What | n |
|---|---|---|
| `REPORT_ois_enhancer_run.md:271` | "1.995 ± 0.908" — is SD (`enrichment_sd`); SEM would be 0.214 | 18 regions |
| `infra_borzoi_grelu/REPORT.md:1077-1079` | "+0.022±0.146 / +0.244±0.059 / +0.080±0.022" — SD (`specificity_probe_quiescent.*.std`) | 15 / 9 / 9 samples |
| `REPORT_generation_probe.md:160` | "2.151 ± 1.501" — SD (`ism.motif_enrichment_sd`) | **n = 2 regions**; the report already flags this NOT SUPPORTED — never quote without that caveat |

### n for every reported mean with a dispersion

| Quantity | n is over |
|---|---|
| Gate 1 / 1b adjusted Δ, spacing, H1c contexts, H2 variants | **10** random 200 bp background seeds |
| H1b ablation effects | **24** regions per response stratum (48 AP-1 total, 24 NF-κB) |
| H1a ISM enrichment | **18** regions (16 with an AP-1 match) |
| Shuffled-label nulls, OIS enhancer run | **200** permutations |
| Shuffled-label nulls, expression phase | **20** permutations (later re-run at 200) |
| ISM PWM enrichment, OIS run | **18** regions |
| ISM motif enrichment, generation probe | **2** regions |
| Quiescence specificity probe | **15 / 9 / 9** samples |
| H3 motif-vs-assay, decile frequencies | **111,671** regions (16,909 AP-1-positive) |

---

## ITEM 2 — the additivity result now has a stored artifact

Recomputed from `output/ois_gates/results/gate1_additivity_raw.csv` by
`recompute_additivity_contrast.py`. Outputs: `additivity_contrast.json`,
`additivity_contrast_per_seed.csv`, `additivity_contrast_LOG.txt`.

**What was compared.** Best combination `NF-κB + ETS, clustered, K = 12` against best single
`NF-κB alone, K = 12`. Both place **12 motif copies on the identical slot grid**, so this is a
matched-load contrast. Per seed, `adjusted = Δ(real) − Δ(scrambled)`, each Δ taken against that
seed's own empty cassette — the backgrounds are shared across conditions, so the pairing is exact.

| | mean | SD | SEM | n |
|---|---|---|---|---|
| pair (NF-κB+ETS clustered K12) | +0.6657 | 0.1541 | 0.0487 | 10 |
| single (NF-κB K12) | +0.6018 | 0.1549 | 0.0490 | 10 |
| **paired difference** | **+0.0639** | **0.1309** | **0.0414** | **10** |

95% CI of the difference **[−0.0297, +0.1575]** — contains zero.

**Test used, and whether it is the right one.** The original used a **paired t-test on 10
seeds**: `t(9) = +1.543, p = 0.1571`. Pairing is correct — seeds are independent random
backgrounds and every condition sees the same background. n = 10 is too small to check
normality usefully, so two distribution-free alternatives were run:

| Test | Statistic | p |
|---|---|---|
| Paired t (original) | t(9) = +1.543 | **0.1571** |
| Wilcoxon signed-rank | W = 14.0 | **0.1934** |
| Exact sign test | 7/10 seeds favour the pair | **0.3438** |

**No materially different answer.** All three fail to reject at 0.05, and the CI contains zero.
The conclusion — the combination is not distinguishable from the best single motif at matched
load — does not depend on the parametric assumption. The non-parametric tests are, if anything,
weaker, so the t-test is the *most* favourable of the three to the combination.

---

## ITEM 3 — 12,429 versus 12,432: resolved, and not as expected

**The correct value is 12,429. The authoritative file is
`output/enhancer_test/step2_stats.json`. The directional-agreement figure is 93.24%, not 93.0%
and not 93.3%.** Provenance note: `provenance_directional_agreement.md`.

**Method.** Rather than guess at a filtering difference, the quantity was recomputed from the
stored per-anchor table `output/enhancer_test/region_responses.csv` using the criterion in
`scripts/52_enhancer_cross_mechanism.py:223`.

| Source | Same-direction count | % of co-differential |
|---|---|---|
| Recomputation from `region_responses.csv` | **12,429** | **93.2408%** |
| `step2_stats.json` (`n_same_direction`) | **12,429** | 93.2408% |
| `STEP2_LOG.txt:81` (same run) | 60.6% of all diff-ATAC → **12,429** | — |
| `REPORT_enhancer_cross_mechanism.md:225` | **12,432** | stated 93.3% |
| Claim inventory | — | stated 93.0% |

**It is not a tie or boundary difference.** Both were tested explicitly:

- Exact-zero signs among the differential-ATAC rows: **0**. There are no ties for `np.sign` to
  resolve differently.
- Switching the threshold from `>=` to strict `>`: **no change at all** (n_diff 20,515,
  co-differential 13,330, same-direction 12,429 either way).

So the three plausible mechanisms — tie handling, boundary inclusion, a different filter — are
all ruled out. **12,432 is not reproducible from the stored data by any variant tested.** The
script never computes 93.3% either; that percentage exists only in the report, derived from the
larger count.

**Assessment:** a report-side transcription error (12,429 → 12,432), with the 93.3% derived
from it. The machine-written JSON and the run log agree with each other and with a fresh
recomputation. Use **12,429 / 13,330 = 93.24%**.

Secondary observation worth keeping: at the stricter `|log2FC| ≥ 1.0` threshold the agreement is
**3,373 / 3,376 = 99.91%**, which is a stronger statement than the 0.5 threshold gives.

---

## ITEM 4 — 8,432 versus 8,459: the audit's explanation was wrong

**The manuscript should use 8,459 windows for 111,671 regions = 13.201 regions per window.**
Provenance note: `provenance_window_count.md`.

The number audit attributed the gap to "the covering being computed before 6
segmental-duplication regions were dropped, with the cache re-deriving it". **That cannot be
right on its face**: dropping regions can only hold the required window count equal or lower
under the same algorithm, yet the count *rose* by 27 while the region count *fell* by 6.

**What was actually done.** Both covering algorithms were re-implemented from their source and
run on the **same** 111,671 cached regions (centres from
`ois_enhancer_run/ois_cache/region_window_map.csv`):

| Algorithm | Source | Windows on the same 111,671 regions |
|---|---|---|
| Point-covering — region treated as dimensionless | `scripts/71_ois_trainset.py:233-238` | **8,432** |
| Width-aware — full 1 kb region must fit inside the predicted span | `ois_enhancer_run/scripts/ois_cache.py:93` | **8,459** |

Both reproduce their respective reported numbers **exactly**, from an identical region set.
The 6 dropped regions contribute nothing.

**The real cause** is the covering criterion:

- The planner (`71_ois_trainset.py`) advanced by `searchsorted(pos, pos[i] + PRED_SPAN)` — it
  counted a region as covered if its *point* position fell in the 196,608 bp predicted span.
- The executed cache (`ois_cache.py`) requires `c − 500 >= s_lo and c + 500 <= s_hi` — the
  **whole 1 kb region** must lie inside the span, which shrinks the effective span to
  196,608 − 1,000 = **195,608 bp** and costs 27 extra windows genome-wide.

Chromosome-end clamping (`ois_cache.py:86`) contributes **nothing**: the width-aware run above
was unclamped and already gives exactly 8,459.

**Status of the two numbers.** 8,432 was a *planning estimate* that under-counted because it
ignored region width; it appears in `REPORT_ois_enhancer_trainset.md` at lines 14, 20, 279 and
368, in `step2_stats.json`, and in `window_covering_per_chrom.csv` (which sums to exactly
111,677 / 8,432). 8,459 is what was executed, recorded in `ois_cache/sanitize_stats.json` and
confirmed by `windows.csv` (8,459 rows) and `region_window_map.csv` (111,671 rows).

The "13.2× compression" headline survives either way (13.244 vs 13.201), but the two are not
the same quantity and should not be mixed.

---

## ITEM 5 — liftover attribution

Correct attribution: `liftover_attribution.md`. **Three builds, three different liftovers.**
The 99.92% / 99.96% pair belongs to the enhancer test, not the OIS training set.

| Build | Liftover | Midpoint success | Failures | Span integrity | Source |
|---|---|---|---|---|---|
| **Enhancer cross-mechanism test** (237,824 anchors) | hg19 → hg38 | **237,628 / 237,824 = 99.92%** | 196 (0.08%) | both ends lift 4,996/5,000 = 99.92%; of those, span within 50 bp of 1 kb **4,994 / 4,996 = 99.96%** | `enhancer_test/harmonization_stats.json`; `STEP1_LOG.txt:30-31` |
| **Replicative ceiling test** (same 237,824 anchors, reused) | hg19 → **hg18** | **237,512 / 237,824 = 99.87%** | 312 (0.13%) | both ends lift 4,995/5,000; span within 50 bp **100.00%** | `replicative_ceiling_test/quantification_stats.json`; `STEP1_LOG.txt` |
| **OIS enhancer training set** (112,099 summits) | hg19 → hg38 | **112,078 / 112,099 = 99.98%** | 21 (0.02%) | **no span-integrity check was run** | `ois_enhancer_trainset/region_stats.json`; `STEP1_LOG.txt:[3]` |

Two consequences for Methods:

1. The training set's liftover is **99.98%, 21 failed** — and those 21 are the same 21 that
   appear in the loss breakdown.
2. **There is no 99.96%-equivalent for the OIS build.** A span-integrity check was never
   performed on it. Do not carry the enhancer-test span figure across; either state that no
   span check was run, or run one.

Also note the 99.96% is conditional — it is 4,994 of the 4,996 anchors whose *both ends*
lifted, not 99.96% of all anchors.

### Other Methods figures checked for the same error

| Figure | Correct build | Cross-contamination? |
|---|---|---|
| 112,099 consensus → 111,671 surviving | OIS trainset | **No** |
| Loss breakdown 247 / 154 / 21 / 6 (0.38%) | OIS trainset | **No** — all four appear in `ois_enhancer_trainset/STEP2_LOG.txt` §5 and §7 |
| 237,824 anchors | Enhancer test; **reused verbatim** by the replicative ceiling test | **No**, but say "the same anchor set" when it appears under the replicative results |
| Cross-dataset reproducibility +0.637 | OIS trainset (thr 0.5, n = 83,716) | **No** |
| Within-OIS ceiling +0.653 (n = 88,530) | Enhancer test | **No** |
| **"35% show no H3K27ac change; 93.x% directional agreement"** | **Enhancer test** | **YES — filed under "training set build" in the claim inventory.** The OIS trainset's own ATAC concordance is a *different, inverted* statistic: 25,082/57,169 = 43.9% of differential-**H3K27ac** regions are also differential in ATAC, 40.8% agreeing in direction (`ois_enhancer_trainset/STEP2_LOG.txt` §4). Do not interchange them — they condition on opposite assays. |
| **Consensus-region criterion** | differs by build | **Worth stating explicitly.** Enhancer test: 1 kb **midpoint**-centred anchors from a **combined ATAC + K27ac** consensus, ≥2 sample support, 500 bp dedup. OIS trainset: 1 kb **summit**-centred regions from **K27ac peaks only**, summit supported by ≥2 of 4 samples, 250 bp summit cluster. This difference is already documented in `ois_enhancer_trainset/STEP2_LOG.txt` §2 as the reason +0.653 and +0.637 are not expected to be equal. |

---

## ITEM 6 — documented compute cost

Full table: `compute_cost_table.csv` (15 rows: 10 GPU sessions, 5 local zero-cost).
Rate is $0.75/hour throughout, user-stated, on every pod.

| Session | What ran | Elapsed | Cost | Source |
|---|---|---|---|---|
| Pod 1 | Task 1 contract smoke test + Task 2 sweep (118 genes) | ~60 min | $0.75 | stored |
| Pod 2 | venv + gReLU, genome setup, corrected Task 2 | ~44 min | $0.55 | stored |
| Pod 3 | gReLU install, genome staging, Tasks 3 + 4 | ~73 min | $0.91 | stored |
| Pod 4 | fine-tune calibration pass (nothing trained) | ~41 min | $0.51 | stored |
| **Pod 5** | full trunk cache, 19,652 genes + head training + setup | **cache 2.53 h measured; setup and training untimed** | **≥ $1.90** | **reconstructed** |
| Transfer | expanded cache 13,211 genes (1.66 h) + transfer heads | ~2.3 h | $1.72 | stored |
| Generation probe | setup + ISM and motif-insertion probe | ~0.77 h | $0.58 | stored |
| Permutations | 200-permutation null on cached embeddings | ~0.42 h | $0.32 | stored |
| OIS enhancer run | setup + cache 1.16 h + training/200 controls + probe 1.7 h | 3.53 h | $2.65 | stored |
| OIS gates | setup + Gate 1/1b + Gate 2 + figures | 0.857 h | $0.64 | stored |
| 5 local sessions | trainset builds, enhancer test, replicative ceiling test, audits | — | $0.00 | stated local-only |

- **Nine sessions carry a stored cost figure, totalling $8.63.**
- **Pod 5 has no stored cost.** It is reconstructed from one stored elapsed figure —
  `CACHE COMPLETE 2.53 h` in `execution_results_finetune/cache_pass_clean.log`. Setup and head
  training on that pod were never timed, so **$1.90 is a lower bound**, not an estimate of the
  session.
- **Documented total: ≥ $10.53.**

Every per-session figure in the project is itself derived from SSH command timestamps, not from
a pulled RunPod invoice — each source report says so.

The stale aggregate in `infra_borzoi_grelu/REPORT.md:1181` ("Cumulative approximate spend across
all four pod sessions: ~$2.72") is correct for pods 1–4 only and should not be quoted as a
project total.

**Recommendation: omit the figure from the manuscript body.** It is not a scientific result, it
cannot be stated precisely (pod 5 is a lower bound, and no figure anywhere is invoice-backed),
and a referee gains nothing from it. If a reproducibility or data-availability statement wants a
cost signal, the defensible wording is: *"approximately $10.50 of GPU time at $0.75/h across ten
sessions; see `output/audit_fixes/compute_cost_table.csv`"* — with "approximately" and the
artifact reference both doing real work. Do **not** write "$8.63" (it silently omits the
single largest compute session) and do **not** write a bare "$10.53" (it implies a precision the
sources do not support).

---

## ITEM 7 — GSE254358 error payloads labelled

`data/GSE254358/README.md` is in place. Files kept, nothing deleted.

Live re-verification on 2026-10-04 reproduces the original finding and adds one control:

| Check | Result |
|---|---|
| `esearch db=gds term=GSE254358[ACCN]` | **count = 0**, empty idlist |
| FTP `…/geo/series/GSE254nnn/GSE254358/` | **HTTP 404** |
| Parent range dir `…/geo/series/GSE254nnn/` | **HTTP 200** |
| **Sibling `…/GSE254357/`** (new control) | **HTTP 200** — an accession one lower in the same folder resolves, so the 404 is specific to this accession, not the path |
| `esearch db=bioproject term=GSE254358` | **count = 0** |
| Six transposition typos | all six resolve to real but unrelated series |

**Correction to the number audit's wording.** The audit said the files "decompress to
`mailto:webadmin@ncbi.nlm.nih.gov`". They are **not gzip at all** — `file` reports
`XML 1.0 document text`, and `gzcat` fails with `not in gzip format`. They are uncompressed
Apache *"Object not found!" Error 404* XHTML pages, saved under `.gz` names. The mailto line is
the page's `webmaster` link, which `suppl_listing.txt` scraped out.

| File | Size | md5 |
|---|---|---|
| `GSE254358_family.soft.gz` | 990 B | `af4cfff178d8056c25c148f96e742185` |
| `GSE254358_series_matrix.txt.gz` | 990 B | `af4cfff178d8056c25c148f96e742185` (byte-identical to the above) |
| `suppl_listing.txt` | 66 B | `58320622d970ec45733291510b51c32a` |

sha256 for all three is in the README.

---

## ITEM 8 — merged accession log

`accession_log_merged.csv` — **44 accessions**: 29 GEO, 12 BioProject, 3 ENCODE. Columns:
accession, repository, role in project, verification method, date verified, outcome, detail.
Every row re-verified live on 2026-10-04 against `eutils.ncbi.nlm.nih.gov` or
`encodeproject.org` — none trusted from the prior log.

| Outcome | n |
|---|---|
| LIVE | 37 |
| LIVE but unrelated (the six rejected GSE254358 typo candidates) | 6 |
| **DOES NOT EXIST** (GSE254358) | 1 |

The six accessions the old `output/tables/accession_verification_log.txt` omitted — GSE106146,
GSE205692, GSE220545, GSE146585, GSE99028, GSE254358 — are all present. So are nine that neither
log had: the **GSE146585 subseries actually quantified by the replicative ceiling test**
(GSE146567 IMR90 replicative, GSE146559 BJ replicative, GSE146563 IMR90 OIS same-lab), three
further subseries named in that report (GSE146566, GSE146568, GSE146582), GSE74328, GSE169767
and GSE105937.

### Coverage cross-check

Every `GSE…`, `PRJNA…` and `ENCSR…` token in every `.md`, `.py`, `.sh`, `.txt`, `.json` and
`.csv` in the project was extracted and diffed against the log.

- Raw scan: 5,316 distinct accessions. **5,267 of those come from one file**,
  `infra_borzoi_grelu/execution_results/task1_output_tracks.csv` — the Borzoi model's own
  7,611-track provenance table. Those are the pretrained trunk's training-track accessions, not
  datasets this project obtained or cites, and they do not belong in a Methods accession log.
- Excluding that one bulk table: **49 distinct accessions on the Methods surface; the log covers
  44; 0 in the log are unused.**
- **The 5 not in the log are all ENCODE internal cross-references inside a downloaded payload**,
  `data/ENCSR200OML/experiment.json`: `ENCSR143LJB` (its `related_annotations` entry) and
  `ENCSR535HHO`, `ENCSR571ISC`, `ENCSR585BFP`, `ENCSR938RZZ` (pipeline `reference_files` on
  individual ENCODE file records). None appears in any report, script or table; they are ENCODE's
  own metadata, not project citations.

**No accession used anywhere in the project is missing from the merged log.**

---

## Files written

| Path | What |
|---|---|
| `output/audit_fixes/AUDIT_FIXES.md` | this report |
| `output/audit_fixes/dispersion_audit.csv` | item 1 — 20-row correction table |
| `output/audit_fixes/dispersion_sem_to_sd.csv` | item 1 — SEM→SD for every stored gate value |
| `output/audit_fixes/recompute_additivity_contrast.py` | item 2 — derivation script |
| `output/audit_fixes/additivity_contrast.json` | item 2 — full result |
| `output/audit_fixes/additivity_contrast_per_seed.csv` | item 2 — the 10 per-seed values |
| `output/audit_fixes/additivity_contrast_LOG.txt` | item 2 — run log |
| `output/audit_fixes/provenance_directional_agreement.md` | item 3 |
| `output/audit_fixes/provenance_window_count.md` | item 4 |
| `output/audit_fixes/liftover_attribution.md` | item 5 |
| `output/audit_fixes/compute_cost_table.csv` | item 6 |
| `data/GSE254358/README.md` | item 7 — in place, in that directory |
| `output/audit_fixes/accession_log_merged.csv` | item 8 |

Source reports were **not** edited; they remain the historical record, as instructed.
