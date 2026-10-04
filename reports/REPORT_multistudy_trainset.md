# Expanded multi-study senescence training set

**Built:** 2026-08-23 · local only, no GPU
**Supersedes:** `output/trainset/senescence_trainset.h5ad` (33 tasks × 19,652 genes, WI-38 only)
**Product:** `output/trainset_expanded/senescence_trainset_multistudy.h5ad` — **85 tasks × 13,211 genes**

**Scripts:** `40_fetch_expansion_series.sh` → `41_parse_expansion_metadata.py` →
`42_build_expanded_trainset.py` → `43_emit_expanded_anndata.py` → `44_verify_expanded_trainset.py`

---

## Headline answers to the three questions asked

| Question | Answer |
|---|---|
| Independent biological lineages in the ON arm | **3 cell lines** (WI-38, IMR90, GM21) / **9 independent culture units** — up from 1 cell line and ~3 units |
| Induction mechanisms in the ON arm | **2** (replicative, oncogene-induced). Not 3 — see below |
| Is the effective n materially better than n≈3? | **Yes, roughly 3×, but it is still small.** 9 independent ON units, not the 24 the task count suggests |

**The honest version:** the ON arm grew from 9 tasks to 24, but from ~3 independent
biological units to **9**. Most of the task-count growth is serial timepoints resampled
from the same cultures, which do not add the kind of replication a shuffled-label control
tests. This is a real improvement — three cell lines from three labs now have to agree
before a gene counts as senescence-responsive — but it is a 3× improvement on a very small
number, not a move into well-powered territory.

**And one finding materially changes the plan:** the two OIS studies agree with each other
on which genes respond (Spearman ρ = **+0.531**), but replicative senescence agrees with
neither (ρ = **+0.093** and **+0.142**). The mechanisms do *not* appear to share much
regulatory logic at the transcriptome level. That is the assumption a single binary
contrast would have smuggled in, and it now has a measurement against it.

---

## 1. What actually made it in — and two datasets that did not

The brief named GSE206402, GSE210285 and GSE74324 as additions. Two of the three could not
be used as specified. Both failures are about **what was deposited**, not about the
accessions being wrong — all five candidates were verified live on GEO FTP and via esummary
before any of this.

| Series | Cell line | Mechanism | Deposited | Verdict |
|---|---|---|---|---|
| **GSE175533** | WI-38 | replicative | TPM | **training** (existing primary, rebuilt through the same path) |
| **GSE74324** | IMR90 | OIS (HRAS-V12) | RPKM | **training** — 9 of 72 samples |
| **GSE205692** | GM21 | OIS (H-RAS-G12V) | log2 RMA microarray | **training** — 21 of 34 samples |
| **GSE99028** | IMR90 | etoposide (TIS) | raw counts | **held out as transfer test** |
| GSE206493 | GM21 | OIS | raw counts | **excluded — no control arm exists** |
| GSE210285 | 2BS | replicative | ATAC bigWigs only | **excluded — no expression data exists** |

### 🔴 GSE206493 is not the matched control-containing RNA-seq the brief describes

All 12 samples are H-RAS-G12V-overexpressing GM21 at day 14 or 27, differing **only** by
shRNA (`pLKO` non-targeting vs `shPOU2F2 #324` vs `#325`). There is no empty-vector, no
proliferating, and no untreated arm anywhere in the series. The single processed file is
named `GSE206493_Counts_POU2F2_KD.txt.gz`, which says the same thing. **No ON-vs-OFF
contrast can be formed from it.**

The empty-vector-vs-RAS contrast in that study exists only in **GSE205692** — the
undisclosed 4th subseries the prior audit found (34 samples, Affymetrix HTA 2.0,
`PRJNA847199`). That is what I used instead. GSE206402 itself is ATAC-seq; peaks cannot
supply per-gene expression targets, so it cannot contribute to an expression-target set
regardless.

### 🔴 GSE210285 has no expression data at all

4 samples, all ATAC-seq (2 Growing, 2 Senescence). The series summary describes RNA-seq and
ChIP-seq, but none was deposited under this accession, and a GEO search for a sibling 2BS
senescence RNA-seq series returned no such record. Only bigWigs — no peak files either.
Excluded, as the brief instructed.

### GSE74324 — 9 of 72 samples used

Confirmed as the 8-arm perturbation panel the prior audit described. Only the shRenilla
no-treatment arm (`P_Ren_1-3`, `PQ_Ren_1-3`, `H_Ren_1-3`) is a clean unperturbed P/Q/S
baseline. The `_DMSO` arm is the vehicle control of the JQ1 sub-experiment — a different
protocol branch — so it is **not** pooled in as a second set of replicates. The other 63
samples are knockdown or drug conditions.

---

## 2. Genome build — no liftover is needed, and that is not luck

| Series | Declared build |
|---|---|
| GSE175533 | hg38 |
| GSE74324 | **hg19 / GRCh37** |
| GSE99028 | **hg19** |
| GSE205692 | **hg19 / GRCh37** |
| *(GSE206402, excluded)* | hg19 |
| *(GSE210285, excluded)* | hg38 |

Three of the four included studies are hg19. **No liftover was performed and none is
required**, because no deposited coordinate enters the training set:

- targets are per-gene expression keyed by **gene symbol**, never by position;
- `.var` intervals are built independently from **GENCODE v44 (hg38)** TSS coordinates;
- the two deposited coordinate columns that exist (`GSE74324.start_position/end_position`,
  and GPL17586's `seqname/start/stop`) are **not read** by the loaders — GSE74324 uses only
  `external_gene_name` + values, GSE205692 only `gene_assignment` symbols.

Verified by spot-check: `.var` gives CDKN1A `chr6:36,676,460` (hg38; the hg19 value is
~32 kb away) and MKI67 `chr10:128,126,423` (hg38; hg19 is ~1.8 Mb away).

This matters because the brief's concern was well-founded — 524,288 bp windows are exactly
the case where the commonly-cited ~99.9% liftover success rate does not apply, since that
figure is for point variants. Symbol-keyed targets sidestep the problem entirely rather than
surviving it. **0 intervals lost to liftover, because there was no liftover.**

The residual annotation cost is 147 symbols with no GENCODE v44 gene record (dropped), plus
53 + 26 whose 524 kb window would run past a chromosome end or before position 0 (dropped,
never clipped — clipping would shift the TSS off-centre and break the contract).

---

## 3. Quantification harmonisation

Four studies, **four different quantification pipelines**. Nothing was merged raw.

| Study | Deposited | Conversion | Exact? |
|---|---|---|---|
| GSE175533 | TPM | `log1p(TPM)` | — |
| GSE74324 | RPKM | `TPM = RPKM / ΣRPKM × 1e6`, then `log1p` | **exact rescale** |
| GSE99028 | raw counts | ÷ GENCODE v44 union-exon length → TPM → `log1p` | approximate |
| GSE205692 | log2 RMA intensity (sva-corrected by submitters) | × ln2 onto a natural-log axis | **not an expression unit** |

The GSE205692 rescale is cosmetic — it makes the magnitudes comparable but does **not** make
microarray intensity the same measurement as TPM. Only the correction in §5 addresses that,
which is why the PCA is reported before *and* after.

Microarray values come from the submitters' own processed series-matrix table, not from
re-running RMA on the CEL files locally — that would have introduced a *fifth* pipeline for
no benefit.

---

## 4. Gene identification across studies

**Different annotation vintages produce different gene sets**, and this cost real genes
before it was fixed. GSE205692's HTA 2.0 platform is ~2013-era and calls the SASP chemokine
`IL8`; the RNA-seq studies call it `CXCL8`. The naive intersection therefore **deleted one of
the five marker genes this build validates on**.

Fixed by mapping every study onto current HGNC-approved symbols (45,044 approved symbols;
56,651 unambiguous alias/previous mappings). Old symbols reused by more than one current gene
(1,441 of them) are **left unmapped**, not guessed — guessing would silently merge two
different genes' measurements.

| Study | Renamed | Collisions left alone | Merged after rename |
|---|---|---|---|
| GSE175533 | 1,077 | 75 | 6 |
| GSE74324 | 1,686 | 53 | 3 |
| GSE99028 | 250 | 27 | 3 |
| GSE205692 | 1,378 | 24 | 3 |

**Intersection: 13,274 → 13,669 symbols** after harmonisation, and CXCL8 is back.

| Study | Symbols | Lost to the intersection |
|---|---|---|
| GSE175533 | 33,481 | 19,812 (59.2%) |
| GSE74324 | 31,502 | 17,833 (56.6%) |
| GSE99028 | 23,549 | 9,880 (42.0%) |
| GSE205692 | 17,784 | **4,115 (23.1%)** |

**GSE205692 is the binding constraint.** The microarray covers only 17,784 genes, so
including it costs the set ~6,400 genes relative to an RNA-seq-only intersection. That is a
real price paid for the third cell line — stated here rather than buried.

Final universe after intervals + expression filter: **13,211 genes** (down from 19,652 in the
single-study build — a 33% reduction, entirely attributable to requiring agreement across
four annotation vintages).

### Excel date corruption

Audited on **every** study, not just the one where it was known. GSE175533's 27 corrupted
symbols (MARCH1-11, SEPT1-14, DEC1) are dropped upstream as before. **Residual date-like
symbols in all four loaded matrices: 0.** No new supplementary file carries the defect.

---

## 5. 🔴 Batch structure — the largest threat, quantified

### Before correction: study identity is almost the entire signal

| Space | PC | Var explained | R²(study) | R²(senescence status) | Dominated by |
|---|---|---|---|---|---|
| raw pooled log | PC1 | **74.7%** | **0.997** | 0.136 | **STUDY** |
| raw pooled log | PC2 | 13.1% | 0.988 | 0.037 | **STUDY** |
| raw pooled log | PC3 | 3.2% | 0.003 | 0.038 | status |
| raw pooled log | PC4 | 2.0% | 0.110 | 0.148 | status |

**PC1 carries 74.7% of all variance and is 99.7% explained by which study a sample came
from.** Senescence status explains 13.6% of it. Pooling these raw would have handed the
model study identity as the single most learnable feature in the data — and it would have
looked like an improvement while failing in exactly the same way as before.

### Correction applied: arm-balanced per-study per-gene centering, then per-study scale equalisation

For each study *s* and gene *g*, the reference subtracted is the **unweighted mean of that
study's ON-arm mean and OFF-arm mean**; the result is then divided by that study's pooled SD.

Why this form:

- Subtracting a per-study per-gene constant removes the study offset **and simultaneously
  every per-gene constant nuisance factor** — gene length, annotation-version mismatch, probe
  affinity, GC/3′ bias. Those are precisely the terms that differ between TPM, RPKM→TPM,
  counts→TPM and microarray intensity, so the four pipelines are reconciled *by
  construction*. It also means GSE99028's approximate union-exon lengths cannot bias the
  corrected target — the error is a per-gene constant within the study and cancels exactly.
- Balancing on **arm means** rather than the grand mean makes the reference independent of
  how many replicates each arm has. ON:OFF ratios here are 9:15, 3:3, 12:9, 2:2 — grand-mean
  centering would shift each study by a composition-dependent amount and reintroduce study
  structure through the back door.
- **Rejected:** centering on the OFF arm alone (forces every OFF sample to ~0, so all OFF
  tasks would carry only replicate noise and the model would be trained to fit it).
- **Rejected:** per-gene z-scoring (dividing by a within-study SD estimated from n=2–3
  amplifies noise for invariant genes).

| Study | Pooled SD of centered values |
|---|---|
| GSE175533 | 0.353 |
| GSE74324 | 0.446 |
| GSE99028 | 0.373 |
| **GSE205692** | **0.164** |

GSE205692's dynamic range is ~2.4× compressed relative to the RNA-seq studies — the
microarray + sva signature — which is what the scale equalisation step addresses.

### After correction

| Space | PC | Var explained | R²(study) | R²(status) | Dominated by |
|---|---|---|---|---|---|
| study-centered | PC1 | 31.0% | 0.214 | 0.055 | study |
| study-centered | PC2 | 20.2% | 0.208 | 0.193 | study |
| **centered + scaled** | PC1 | 23.9% | 0.217 | 0.035 | study |
| **centered + scaled** | **PC2** | **20.8%** | 0.202 | **0.265** | **STATUS** |
| centered + scaled | PC4 | 7.0% | 0.057 | 0.160 | status |

**R²(study) on PC1 falls from 0.997 to 0.217**, and senescence status becomes the dominant
explanation on PC2. Figure: `output/figures/multistudy_batch_pca.png`.

**Stated plainly: batch structure is strongly suppressed but not eliminated.** R² of ~0.20–0.22
for study on the top two PCs is not zero, and no correction of this kind can make it zero when
one of the four studies is a different measurement technology. Any training run on this set
**must** keep a study-shuffled control alongside the label-shuffled one — if a model can be
trained to predict study identity from these targets, the residual 0.22 is doing work.

`.X` carries the corrected target. `layers['log_native']` (uncorrected) and
`layers['study_centered']` (centered, unscaled) are both retained, so nothing is destroyed
and the correction can be revisited without a rebuild.

### ⚠️ Contract consequence: the softplus head must change

The corrected targets are **centered, therefore signed** (range −14.83 to +16.18). The
existing head is `ConvHead + softplus`, which cannot represent negative values. **Training on
`.X` with a softplus output head will silently fail** — it would floor the entire negative
half at zero. The output nonlinearity must be changed to linear. This is recorded in
`uns['contract']['NOTE_softplus']`. It is a required consequence of batch correction, not an
oversight.

---

## 6. Per-study marker sanity checks — all four pass

log2 fold change, ON arm vs **that study's own** OFF arm, computed per study independently.

| Study | CDKN1A ↑ | IL6 ↑ | CXCL8 ↑ | LMNB1 ↓ | MKI67 ↓ | *CDKN2A (report only)* |
|---|---|---|---|---|---|---|
| GSE175533 | +0.51 ✅ | +3.01 ✅ | +2.08 ✅ | −3.33 ✅ | −3.13 ✅ | *+0.20* |
| GSE74324 | +1.50 ✅ | +4.29 ✅ | +9.24 ✅ | −4.11 ✅ | −4.02 ✅ | *+1.43* |
| GSE99028 | +2.53 ✅ | +2.03 ✅ | +4.97 ✅ | −3.34 ✅ | −3.78 ✅ | *+1.29* |
| GSE205692 | **+0.15** ✅ | +1.02 ✅ | +5.43 ✅ | **−0.47** ✅ | **−0.83** ✅ | *+0.12* |

**Studies with ≥1 marker failure: none.** No study was excluded on marker grounds.

CDKN2A behaves as prior work established — it moves least of all six genes in the two
replicative/weak-effect studies (+0.20, +0.12) and was not used as a criterion.

**⚠️ GSE205692 passes on direction but weakly on magnitude.** Its CDKN1A response is +0.15
against +0.51/+1.50/+2.53 elsewhere, and its arrest markers are ~4–8× smaller. Its marker
evidence rests mainly on CXCL8 (+5.43). Partly the compressed microarray range, partly the
submitters' sva correction. **This study contributes a weaker-magnitude contrast than its
21-sample count suggests** — it is the largest contributor by task count and among the
weakest by effect size.

### GSE205692 day window — and a correction to my own reasoning

This study is about **escape** from OIS, so late RAS timepoints may have resumed
proliferating. Rather than assert a window, I measured it (SASP score = mean z of
CDKN1A/IL6/CXCL8; arrest score = mean z of −LMNB1/−MKI67):

| RAS day | n | SASP | Arrest | Composite | Assigned |
|---|---|---|---|---|---|
| 8 | 3 | **−0.18** | **+1.54** | +1.36 | excluded |
| 13–25 | 12 | +0.21 … +1.27 | −0.42 … +1.88 | +0.13 … +2.48 | **ON** |
| 32–56 | 10 | −0.39 … −0.12 | −0.36 … +0.50 | −0.76 … +0.11 | excluded |
| *empty vector* | 9 | −0.64 | −0.76 | −1.40 | OFF |

Per-marker log2FC vs empty vector:

| Group | CDKN1A | IL6 | CXCL8 | LMNB1 | MKI67 |
|---|---|---|---|---|---|
| day 8 | +0.14 | **−1.05** | +1.45 | −0.71 | −0.77 |
| days 32–56 | +0.08 | +0.25 | +0.60 | −0.39 | −0.53 |
| **days 13–25 (ON)** | +0.15 | +1.02 | **+5.43** | −0.47 | −0.83 |

- **Days 32–56 confirmed as escape:** CXCL8 +0.60 vs +5.43, MKI67 −0.53 vs −0.83. Attenuated
  signature, consistent with a partially escaped population. Correctly excluded.
- **Day 8: my stated reason was wrong.** I excluded it a priori as "arrest not yet
  established". Measurement shows day 8 is the *most arrested* timepoint in the series
  (arrest z +1.54, higher than any ON day). What it lacks is the SASP — **IL6 is −1.05, below
  the empty-vector control**. The exclusion stands, but on the opposite ground: day 8 is
  arrest-without-SASP, i.e. early incomplete senescence, and it **fails the pre-registered
  "CDKN1A/IL6/CXCL8 up" criterion outright**. Including it would place a phenotypically
  different arm under the same ON label as three fully SASP-positive arms.

---

## 7. GSE99028 — held out as the transfer test, not trained on

**Decision: hold out, joining the GSE175533 RIS arm.** Reasoning, stated as the brief asked:

- Etoposide and doxorubicin are **both topoisomerase-II poisons**. This is the closest
  available proxy to the project's eventual doxorubicin validation, which makes it worth more
  as a test than as training data.
- It is a **2-vs-2** contrast. As training data it adds 2 ON tasks — negligible against the
  failure mode being fixed, which is inability to separate signal from noise at low n.
- Holding it out keeps **one induction mechanism entirely unseen**. Given the §8 finding that
  replicative and OIS responses barely correlate, an unseen-mechanism test is the only honest
  way to check whether the model has learned senescence rather than two study-specific
  programmes.

**Caveat, not correctable from these data:** the two proliferating controls carry
`treatment: None` with no shRNA, while the ON arm carries a non-targeting shRNA. The contrast
therefore confounds etoposide with shRNA vector + selection. Small next to the etoposide
effect, but not zero.

**Assumption flagged — column-to-GSM mapping.** GEO does not state which matrix column is
which GSM. `PD29`/`PD34` map by literal string. `SC-etoposide-rep{1,2}` is read as the SOFT's
`sh-NTC etoposide-replicate {1,2}` ("SC" = scrambled control): it is the only remaining
pairing, GSM order matches column order left-to-right, and SOFT `mutation:` is `Wild-type` for
the SC pair vs `cGAS shRNA knockdown` for the shcGas pair. The `shcGas` samples are excluded —
cGAS knockdown suppresses the SASP arm, so they are a deliberately partial phenotype.

Held-out contents: **GSE99028** 2 ON + 2 OFF (`role='test_transfer'`), **GSE175533 RIS**
15 senescent + 3 baseline (`role='holdout_ris'`).

---

## 8. 🔴 Binary contrast vs multi-task — recommendation, with a measurement behind it

**Recommendation: multi-task.** Keep one output task per sample (as the existing contract
already does) and define the senescence contrast at *evaluation* time.

The brief correctly noted that a single binary contrast assumes the mechanisms share
regulatory logic, and that the project has not established this. It is now measured. Spearman
correlation of the per-gene ON−OFF response, between studies, in the corrected space:

| Pair | ρ |
|---|---|
| GSE205692 (GM21, OIS) vs GSE74324 (IMR90, OIS) | **+0.531** |
| GSE175533 (WI-38, replicative) vs GSE74324 (IMR90, OIS) | **+0.142** |
| GSE175533 (WI-38, replicative) vs GSE205692 (GM21, OIS) | **+0.093** |

**The two OIS studies — different cell lines, different labs, different platforms — agree
moderately. Replicative senescence agrees with neither.** ρ ≈ 0.1 is very weak agreement.

This is a genuine, useful result: the cross-study consistency check was impossible in the
single-study build, and it says the shared-logic assumption is **not** supported for
replicative-vs-OIS. Collapsing to one binary label would force the head to fit one weight
vector to two largely uncorrelated response programmes, and the most likely outcome is that
it fits the larger/cleaner one (OIS, 15 of 24 ON tasks) and calls it senescence.

**Tradeoff, stated:** multi-task costs more head parameters and yields no single "senescence
score" without choosing a contrast. It also cannot rescue the fact that ρ ≈ 0.1 means there
may be little *common* signal to learn. The upside is that it keeps the mechanisms separable,
lets per-mechanism performance be read off directly, and makes the GSE99028 transfer test
interpretable — under a single binary contrast, failure on TIS would be uninterpretable
because you could not tell which training mechanism it failed to generalise from.

Signal survives the correction, so the target is not gutted:

| Marker | Corrected ON−OFF (SD units) | Native log2FC |
|---|---|---|
| CXCL8 | +13.16 | +5.33 |
| IL6 | +5.35 | +3.36 |
| CDKN1A | +1.00 | +0.27 |
| CDKN2A | +0.66 | +0.23 |
| LMNB1 | −4.64 | −1.63 |
| MKI67 | −5.07 | −1.33 |

Genome-wide: mean ON−OFF −0.055, sd 0.888, **1,878 / 13,211 genes move >1 SD**.

---

## 9. The object

```
output/trainset_expanded/senescence_trainset_multistudy.h5ad
  .X                      (85, 13211) float32   corrected target, range [-14.83, +16.18]
  .layers['log_native']   uncorrected log expression
  .layers['study_centered'] centered, not scaled
  .obs   study, matrix_column, sample, cell_line, mechanism, mechanism_class,
         platform_kind, quantification, cls, cls_binary, role, detail, bio_unit
  .var   chrom, start, end, strand, tss, gene_type, split
  .uns   contract, target_transform, splits, arm_rationale, quantification_by_study,
         symbol_harmonisation, batch_correction, provenance
```

**Preserved from the existing build:** 524,288 bp TSS-centered intervals (verified: all
widths 524288, all TSS offsets 262144); float32 targets; **chr9 + chr6 held out entirely as
test**; disjoint chromosome lists.

| Split | Chromosomes | Genes |
|---|---|---|
| train | 19 chroms (chr1–5, 7, 10–15, 17–19, 21, 22, X, Y) | 10,625 |
| val | chr8, chr16, chr20 | 1,422 |
| test | **chr9, chr6** | 1,164 |

Verified: zero chromosome overlap between any pair of splits; 0 duplicate gene names; 0
duplicate task names; CDKN2A (chr9) and CDKN1A (chr6) both land in **test**, as before.

| Role | Tasks | |
|---|---|---|
| `train_eligible` | 63 | ON 24, OFF 27, quiescent 12 |
| `test_transfer` | 4 | GSE99028: 2 ON, 2 OFF |
| `holdout_ris` | 18 | GSE175533 RIS: 15 senescent, 3 baseline |
| *(dropped)* | *123* | *perturbation arms, escape timepoints, hTERT, transition PDLs* |

Companion CSVs: `multistudy_sample_assignments.csv` (208 rows, all samples with role),
`multistudy_task_table.csv` (85 retained tasks), `per_study_markers.csv`,
`batch_pca_stats.csv`.

---

## 10. Replication accounting — the number that actually matters

The task count is not the number that governs whether a shuffled-label control can be beaten.
Serial timepoints from one culture and replicate wells of one experiment are not independent
draws.

| Study | ON tasks | **Independent culture units** |
|---|---|---|
| GSE175533 (WI-38, replicative) | 9 | **3** |
| GSE205692 (GM21, OIS) | 12 | **3** |
| GSE74324 (IMR90, OIS) | 3 | **3** |
| **Total** | **24** | **9** |

- GSE205692's 12 ON tasks come from only **3** independent cultures — the other 9 are
  additional day-timepoints sampled from those same 3.
- GSE175533's 9 ON tasks come from **one WI-38 serial-passage lineage**; the 3 counted units
  are parallel replicate cultures of that lineage. This is the generous reading — even those
  descend from a single WI-38 stock.

OFF arm: 27 tasks, **9** independent units.

| | Before (failed build) | After |
|---|---|---|
| ON tasks | 9 | 24 |
| Cell lines | 1 (WI-38) | **3** (WI-38, IMR90, GM21) |
| Independent culture units | ~3 | **9** |
| Distinct mechanisms | 1 | **2** (replicative, OIS) |
| Independent studies/labs | 1 | **3** |

**"2 mechanisms", not 3.** GSE74324's `HRASV12` and GSE205692's `HRASG12V` are the same
oncogenic-RAS mechanism under two spellings of the same allele. OIS is represented twice — by
two different cell lines, which is genuine lineage replication but not mechanism diversity.
Counting it as 3 would overstate exactly the quantity this build exists to improve.

---

## 11. What this does and does not fix — read before training

**Improved, materially:**
- 1 → 3 cell lines, 1 → 3 labs, ~3 → 9 independent culture units in the ON arm.
- Cross-study consistency is now *measurable* — it was structurally impossible before.
- A gene must move consistently across three independent studies to look senescence-responsive,
  which is a much harder bar for a spurious grouping to clear than the single-lineage version.

**Not fixed, and honest about it:**
1. **9 independent ON units is still small.** A 3× improvement on ~3. This may well not be
   enough to clear a shuffled-label control at p < 0.05; the expansion is real but modest,
   and the brief was right to anticipate that.
2. **The gene universe shrank 33%** (19,652 → 13,211), almost entirely to accommodate
   GSE205692's 17,784-gene microarray coverage. Fewer training pairs per task.
3. **Residual batch structure R² ≈ 0.22** on the top PCs. A **study-shuffled control is now
   mandatory** alongside the label-shuffled one.
4. **GSE205692 is the largest ON contributor (12 of 24 tasks) and among the weakest by effect
   size** (CDKN1A +0.15, LMNB1 −0.47). It is also the only non-RNA-seq study. If the next run
   improves, an ablation without it is the first thing to check.
5. **ρ ≈ 0.1 between replicative and OIS** is the most consequential finding here. It may mean
   there is little shared regulatory logic to learn — in which case a single senescence head
   is the wrong target regardless of how much data is added. Multi-task keeps that visible
   instead of averaging it away.
6. **The softplus head must be changed to linear** before training on `.X`.

**Suggested first check on the next run:** train on OIS only (15 ON tasks, 2 cell lines) and
test on replicative, then the reverse. If the ρ ≈ 0.1 finding is real, that cross-mechanism
transfer will fail — and it is far cheaper to learn that from a held-out split than from
another full fine-tune plus control suite.
