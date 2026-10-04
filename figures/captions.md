# Manuscript figures — captions and provenance

Generated 2026-10-04, locally, no GPU. Every plotted value is read from a stored artifact;
nothing was re-run on a model. Two quantities are *derived* from stored arrays rather than read
from a summary — both are marked below and both reproduce the stored summary exactly.

**Format.** 300 dpi PNG + vector PDF (TrueType embedded, `pdf.fonttype 42`, text remains
editable; 10 of 12 PDFs are fully vector, Figures 5 and S1 each carry one raster XObject,
the hexbin density layer). Full width 180 mm, single column 85 mm. Sans-serif.
**Minimum type size is 7 pt at final size** - verified: no `fontsize` below 7 anywhere in
`figures_main.py` or `figures_supp.py`, and all rcParams tick/legend/label sizes are 7.
Categorical palette validated colourblind-safe with the dataviz validator against a white
surface (3-slot: worst adjacent dE 9.2 deutan; 6-slot: 9.1 protan); every low-contrast slot
also carries a direct label or a distinct marker.

**Dispersion.** A prior audit found four figures with unlabelled SEM bars. Every error bar
in this set is named in its own caption AND on the panel itself.

| Figure | File | Width x height (mm) | Error bars |
|---|---|---|---|
| 1 | `figure1_data_landscape.png` | 85 x 78 | none (counts) |
| 2 | `figure2_reproducibility_matrix.png` | 176 x 112 | none (point estimates) |
| 3 | `figure3_threshold_sensitivity.png` | 80 x 77 | none (point estimates) |
| 4 | `figure4_cdkn2a_promoter_distal.png` | 81 x 80 | none (means over stated regions) |
| 5 | `figure5_model_performance.png` | 175 x 66 | none (distributions / point estimates) |
| 6 | `figure6_design_ceiling.png` | 180 x 73 | **SEM**, n = 10 backgrounds |
| 7 | `figure7_ap1_context_dependence.png` | 176 x 73 | **SEM**, n = 24 regions per bar (B) |
| S1 | `figureS1_atac_vs_k27ac_concordance.png` | 169 x 75 | none |
| S2 | `figureS2_batch_structure.png` | 167 x 67 | none |
| S3 | `figureS3_shared_control_artifact.png` | 158 x 67 | none |
| S4 | `figureS4_exhaustive_3v3_null.png` | 81 x 72 | none |
| S5 | `figureS5_motif_dose_response_all.png` | 177 x 71 | **SEM**, n = 5 backgrounds |

---

## Figure 1 — Data landscape (single column)

**Caption.** Senescence chromatin datasets used in this study, by induction mechanism. Filled
bars give the number of replicates in the senescent arm, hatched bars the matched non-senescent
control arm; colour encodes mechanism. The shaded region and dashed line mark n = 2. Every
H3K27ac arm in every study is n ≤ 2 per side; only GM21 ATAC reaches n = 4 vs 4. GSE106146
contributes two contrasts (replicative and irradiation) that share one proliferating control
arm (†), which is the source of the shared-denominator artifact quantified in Figure S3. WI-38
replicative ATAC (GSE175533) is deposited as a log2 fold change only and has no replicate-level
data, so it is not plotted. This n = 2 ceiling bounds every correlation reported in this paper.

**Sources.** `output/figures_manuscript/dataset_table.csv`, assembled from
`output/replicative_ceiling_test/step2_stats.json` (`n_per_arm`),
`output/replicative_ceiling_test/scripts/60_fetch_replicative_ceiling.sh` (arm → accession),
`output/ois_enhancer_trainset/ois_enhancer_trainset.h5ad` (`obs`: cell line, mechanism,
control arm, n_replicates), `output/ois_enhancer_trainset/STEP1_LOG.txt`, and
`output/audit_fixes/accession_log_merged.csv`.

---

## Figure 2 — The reproducibility matrix (full width) · CENTREPIECE

**Caption.** Per-region Spearman correlation of the senescence response between every pair of
datasets, at the primary response-magnitude threshold |log2FC| ≥ 0.5, grouped by whether the
two datasets share an induction mechanism. Each point is one comparison; n (regions entering
that correlation) is printed beside it. The dashed blue line marks the within-OIS ceiling,
ρ = +0.653 (GM21 vs IMR90 H3K27ac, n = 88,530). Grey shading spans the four cross-mechanism
H3K27ac comparisons, −0.164 to +0.164. **The ringed point is R2** (IMR-90 vs BJ replicative,
ρ = +0.054, n = 134,336): the single comparison with no technical confound — same lab, same
pipeline, same genome build — and therefore the load-bearing value. Two of the three
within-replicative H3K27ac comparisons (R1 +0.144, R2 +0.054) fall inside the cross-mechanism
span; the third (R3 +0.268) sits just above it. No within-replicative comparison comes close to
the within-OIS ceiling. Two datasets of the *same* mechanism are not reliably more concordant
than two datasets of *different* mechanisms. No error bars: each point is a single correlation
over the stated n.

**Sources.** `output/replicative_ceiling_test/replicative_correlations.csv` (R1, R2, R3, A1, A2,
X1–X4 at threshold 0.5); `output/enhancer_test/cross_mechanism_correlations.csv` (a, b, c, d at
threshold 0.5).

> ⚠ **Mismatch flagged — see "Discrepancies" below.** The source report states that all three
> within-replicative values fall inside the cross-mechanism range. R3 = +0.2682 does not.

---

## Figure 3 — Threshold sensitivity (single column)

**Caption.** Spearman correlation as a function of the response-magnitude threshold applied to
both datasets. Within-OIS agreement (solid, blue) is essentially flat across thresholds, +0.640
at |log2FC| ≥ 0 to +0.678 at ≥ 1.0 — the agreement is a property of the data, not of the
filter. The within-replicative comparisons (dashed, orange) behave differently: R1 climbs
five-fold, +0.063 to +0.327, meaning its apparent agreement is manufactured by discarding
low-magnitude regions. R2 and R3 are shown for comparison. Endpoint values are printed at each
line end. No error bars: each point is a single correlation; n at each threshold is in the
source table and falls as the threshold rises.

**Sources.** `output/enhancer_test/cross_mechanism_correlations.csv` (comparison b, all four
thresholds); `output/replicative_ceiling_test/replicative_correlations.csv` (R1, R2, R3, all
four thresholds).

---

## Figure 4 — CDKN2A promoter versus distal (single column)

**Caption.** Mean H3K27ac log2 fold change at the CDKN2A locus, split into promoter-proximal
(filled, n = 10 anchors) and distal (hatched, n = 67 anchors) regions, for all seven dataset
arms that quantify it. Bars are means over the stated anchors; no error bars. **The promoter
falls in every arm without exception**, including the two n = 1 vs 1 arms — the assay is
working and the locus behaves as expected everywhere. **The distal signal does not follow**: it
is inconsistent in sign across the three replicative datasets (−0.17, +0.22, +0.31) and rises
only in the two n = 2 vs 2 OIS datasets (+0.24, +0.49). The n = 1 vs 1 same-lab OIS arm is flat
at both promoter (−0.05) and distal (−0.05), consistent with its single replicate. Because the
promoter result reproduces everywhere, the distal inconsistency cannot be attributed to
measurement failure — it is a property of the replicative response.

**Sources.** `output/replicative_ceiling_test/locus_sanity_checks.csv` (rows `locus == CDKN2A`;
columns `*_mean` for all seven arms).

---

## Figure 5 — OIS model performance (full width, 4 panels)

**Caption.** Performance of the trained OIS enhancer head on held-out IMR90, the cross-study
validation split.
**(A)** Predicted versus measured per-region response, hexagonal density (log colour scale),
n = 100,070 held-out regions. Spearman ρ = +0.5484, which is 86.1% of the +0.637 cross-dataset
reproducibility ceiling measured on the same regions — the ceiling, not 1.0, is the right
reference.
**(B)** Shuffled-label null. Histogram of the Spearman correlation under 200 label
permutations (all 200 values stored); the real value is marked in red. 0 of 200 permutations
reach it, giving the smallest attainable p for 200 permutations, **p = 0.00498 = 1/201**. The
AUC null is reported as summary statistics (mean ± SD, observed range, 0/200) because the
per-permutation AUCs were not saved — see "Data not available" below.
**(C)** Controls, all evaluated on the same n = 26,936 evaluable regions: the trained head
against an untrained head and against a head fed region-to-embedding-shuffled inputs. Both
controls collapse to chance (dotted line, AUC 0.5).
**(D)** ROC on the evaluable subset, n = 26,936 of 100,070 held-out regions (17,240 up,
9,696 down); the remaining 73,134 are ambiguous at |measured log2FC| < 1.0 and are excluded by
the pre-registered evaluation rule. AUC = 0.9070.
No error bars in this figure.

**Sources.** `output/ois_enhancer_run/results/results_train.json` (all summary statistics, the
200 stored permutation Spearmans, controls); `output/ois_enhancer_run/results/pred_rest.npy`
and `idx_rest.npy` (per-region predictions); `output/ois_enhancer_trainset/ois_enhancer_trainset.h5ad`
(measured `IMR90_SEN_vs_PRO`).
*Derived:* panels A and D are computed from the stored prediction arrays rather than read from
a summary. The ROC AUC so computed is **0.9070**, reproducing the stored `auc` 0.907004557 to
four decimal places.

---

## Figure 6 — The design ceiling (full width, 3 panels)

**Caption.** How far motif insertion can move predicted senescence H3K27ac in a 200 bp cassette.
**(A)** Dose response: adjusted Δ (real motif minus the matched scrambled control, paired within
background) against total inserted motif copy number, one line per condition. Error bars are
**SEM over n = 10 independent random backgrounds**. The dashed red line is the delta required to
move a sequence from the median to the 95th percentile of predicted activity, +0.861. Four of
the five conditions peak at K = 12 and decline at K = 18 and K = 20; **ETS peaks at K = 6**, not
K = 12. At K = 20 the cassette is tiled end to end with motif and the response collapses. The
best condition reaches +0.666, or 0.773× the requirement.
**(B)** Arrangement at matched load: the same conditions plotted interleaved versus clustered,
paired within background. Clustered is higher in 8 of 8 condition×load pairs, significantly in
3 (paired t-test, p < 0.05, blue). The largest arrangement effect, +0.294, exceeds the gain from
adding a second motif family.
**(C)** The best combination against the best single motif at matched load (12 copies each),
paired by background (grey lines, n = 10; blue line = means). Paired difference +0.0639,
95% CI [−0.030, +0.157], paired t(9) = +1.54, p = 0.157. Wilcoxon signed-rank p = 0.193 and an
exact sign test p = 0.344 agree. Combining motif families buys nothing measurable over using
more of the best single family.

**Sources.** `output/ois_gates/results/gate1b_dose_curve.csv` (A);
`output/ois_gates/results/gate1_spacing.csv` (B);
`output/audit_fixes/additivity_contrast_per_seed.csv` (C, the 10 per-background values);
`output/ois_gates/results/gate1_results.json` (`required_delta`).

---

## Figure 7 — AP-1 context dependence (full width, 3 panels)

**Caption.** The AP-1 insertion result reverses when the motif is tested in its native context.
**(A) INSERTION.** Change in predicted response from 0 to 8 inserted copies in random 200 bp
background, for four transcription-factor families, each against a mononucleotide-shuffled
control of the same motif at the same positions (hatched). NF-κB, C/EBPβ and ETS raise the
prediction; **AP-1 lowers it** (−0.464) against a flat scrambled control (−0.016).
n = 5 backgrounds × 7 copy numbers per arm.
**(B) ABLATION.** The same motifs removed from their native sites in real genomic regions.
Plotted is the motif-specific effect: shuffling the matched bases minus shuffling an equally
long motif-free stretch in the same window, so a generic "perturbing 10 bp changes the
prediction" component is subtracted. Error bars are **SEM over n = 24 regions per bar**.
Negative means the motif was contributing *positively* in that context. **AP-1's sign flips
with context** — −0.1055 (p = 0.032) in high-response regions, +0.0658 (p = 0.135) in
low-response regions, interaction p = 0.0091 (Welch t-test; Mann-Whitney p = 0.0076). NF-κB is
the positive control and behaves as expected (−0.1404, p < 0.001), validating the assay.
**(C)** Fraction of regions carrying a strong AP-1 match (JASPAR MA0099.3, relative score
≥ 0.95) across deciles of the measured GM21 H3K27ac response, n = 11,168 regions per decile.
AP-1 sites rise from 0.086 in the bottom decile to 0.179 in the top (OR = 2.31, Fisher exact
p = 5.9 × 10⁻⁹⁴); the series is non-monotonic and peaks at decile 6 (0.199).
Panels A and B point in opposite directions, and panel C shows that the genome-wide
distribution agrees with the ablation, not the insertion: the insertion probe was only ever run
in near-zero-response backgrounds, which is exactly the regime where panel B says AP-1 is
negative.

**Sources.** `output/ois_enhancer_run/results/results_probe.json`
(`motif_insertion.summary`, panel A); `output/ois_gates/results/gate2_ap1_results.json`
(`H1b_interaction`, panel B); `output/ois_gates/results/gate2_h3_frequency_by_decile.csv`
(panel C, itself derived from `gate2_motif_calls.csv`).

---

# Supplementary figures

## Figure S1 — Accessibility versus H3K27ac concordance

**Caption.** **ENHANCER CROSS-MECHANISM BUILD** (237,824 GM21 anchors) — not the OIS training
set. **(A)** Per-anchor GM21 ATAC response against GM21 H3K27ac response, same cells and same
RAS-vs-EV contrast, hexagonal density on a log colour scale, n = 54,498 anchors with both assays
measured. Dashed lines mark |log2FC| = 0.5. **(B)** Of the 20,515 anchors differential in ATAC,
7,185 (35.0%) show no H3K27ac change; of the 13,330 that are co-differential, 12,429 agree in
direction = **93.24%**. No error bars. The OIS training set reports a superficially similar
43.9% figure, but that statistic conditions on differential *H3K27ac* and asks about ATAC — the
opposite direction — and the two are not interchangeable.

**Sources.** `output/enhancer_test/step2_stats.json` (`active_vs_open["0.5"]`);
`output/enhancer_test/region_responses.csv` (panel A density).
*Derived and checked:* the same-direction count was recomputed from the per-anchor table and
returns **12,429**, matching the stored `n_same_direction` exactly (the figure asserts this at
run time). Note the source report states 12,432 / 93.3%; see
`output/audit_fixes/provenance_directional_agreement.md`.

## Figure S2 — Batch structure before and after correction

**Caption.** Principal-component structure of the OIS H3K27ac signal matrix (108,174 regions ×
10 samples, two studies), with the variance explained by each PC printed under its label.
**(A)** Raw pooled coverage: PC1 (43.8% of variance) is dominated by senescence status
(R² = 0.777) but PC2 (25.9%) is dominated by study (R² = 0.483) — pooling raw coverage mixes
batch into the second component. **(B)** The per-study log2 fold-change contrasts actually used
as targets: study R² falls to ~0 on every PC (3 × 10⁻³¹ or smaller) while status R² rises to
0.908 on PC1. The "correction" is the target definition itself — each contrast is computed
inside one study, so the study-level offset cancels by construction. No pooled batch correction
was applied, deliberately: the train/validate split is cross-study, so any correction fitted
across both studies would leak validation information. No error bars.

**Sources.** `output/ois_enhancer_trainset/batch_pca_stats.csv`;
`output/ois_enhancer_trainset/STEP2_LOG.txt` §3 (matrix dimensions and the rationale).

## Figure S3 — The shared-control artifact

**Caption.** Spearman correlations on one axis showing that a shared control arm manufactures
agreement. The IR-versus-replicative comparison appears strong at +0.564 (n = 110,699) — but
both contrasts divide by the *same* GSE106146 proliferating sample. Three constructed nulls
show what that alone buys: two entirely unrelated numerators over the same shared denominator
still correlate at +0.523 (n = 180,422). Recomputing IR versus replicative with independent
control arms (comparison I3) gives **+0.115** (n = 138,174). Almost all of the apparent
agreement was the shared denominator. No error bars; n printed beside each bar.

**Sources.** `output/replicative_ceiling_test/shared_denominator_null.txt` (I1 and the three
nulls); `output/replicative_ceiling_test/replicative_correlations.csv` (I3 at threshold 0.5).

## Figure S4 — Exhaustive enumeration of the 3v3 permutation null

**Caption.** Every one of the 20 distinct label partitions of the within-OIS 3-versus-3 contrast,
sorted by AUC; the red line is the real AUC (0.7854). This is an exhaustive enumeration
(`mode = exact_exhaustive`), not a sample: the partition space contains exactly 20 members, and
the true labelling is one of them. Five partitions reach the real AUC, giving an exact
p = 0.25. Even had none reached it, the smallest attainable p would be **1/20 = 0.05** — a 3v3
design cannot produce significance below that threshold regardless of effect size. n = 182
evaluable genes (57 positive, 125 negative). This is why the pooled contrasts, which have a far
larger partition space, were used for the reported permutation tests. No error bars.

**Sources.** `infra_borzoi_grelu/execution_results_permutations/results_perm200.json`
(`within_OIS`: `all_aucs`, `real_auc`, `partition_space`, `mode`, `n_beating_real`,
`empirical_p`, `n_evaluable`).

## Figure S5 — Motif-insertion dose response, all four families

**Caption.** Full version of Figure 7A. Predicted response against inserted motif copy number
for each of the four transcription-factor families, real motif (solid, filled markers) against a
mononucleotide-shuffled control at the same positions (dashed, open markers). Error bars are
**SEM over n = 5 independent random backgrounds** at each copy number. Printed per panel: the
0→8 copy delta for the real motif and for its scrambled control, and the Spearman correlation
between copy number and predicted response. NF-κB gives the largest positive response
(Δ +0.665, ρ = +0.846); AP-1 is the only family with a negative, monotonic dose response
(Δ −0.464, ρ = −0.375) against an essentially flat scrambled control (−0.016).

**Sources.** `output/ois_enhancer_run/results/motif_insertion.csv` (per-sequence predictions:
motif, variant, seed, k, pred); `output/ois_enhancer_run/results/results_probe.json`
(`motif_insertion.summary` for the printed deltas and ρ).

---

# Discrepancies between stored data and the manuscript draft

A prior audit found four mismatches, all stated to be already corrected in the draft. Building
these figures surfaced **one further mismatch that the audit did not catch**, because the audit
checked individual numbers but not a *relational* claim between two of them.

### NEW — Figure 2: "all three within-replicative values fall inside the cross-mechanism range"

`output/replicative_ceiling_test/REPORT_replicative_ceiling.md:24` states:

> "**The three within-replicative values (+0.054, +0.144, +0.268) fall inside the range of the
> four cross-mechanism values (−0.164 to +0.164).**"

Stored values at threshold 0.5 (`replicative_correlations.csv`, `cross_mechanism_correlations.csv`):

| | ρ | inside [−0.164, +0.164]? |
|---|---|---|
| R2 IMR-90 vs BJ | +0.0536 | yes |
| R1 IMR-90 vs GSE106146 | +0.1439 | yes |
| **R3 BJ vs GSE106146** | **+0.2682** | **no — exceeds the upper bound +0.1639** |

**Two of three, not three of three.** The audit verified R3 = +0.268 (MATCH) and the
cross-mechanism range −0.164 to +0.164 (MATCH) separately, and never compared them.

The argument survives intact — no within-replicative comparison approaches the +0.653 within-OIS
ceiling, and the cleanest comparison (R2, no technical confound) is +0.054 — but the sentence as
written is false and a referee checking the two numbers against each other will see it.
Figure 2 and its caption state the accurate version. Suggested replacement wording: *"Two of the
three within-replicative comparisons (+0.054, +0.144) fall inside the cross-mechanism range
(−0.164 to +0.164); the third (+0.268) exceeds it only slightly and remains far below the
within-OIS ceiling of +0.653."*

Also recorded in `output/figures_manuscript/figure_flags.txt`.

### Carried forward from the prior audit — still visible in the source reports

These were reported as corrected in the draft; the **source reports still contain them** and are
not the basis of any figure here:

- `REPORT_enhancer_cross_mechanism.md:225` — 12,432 / 93.3%; stored value is 12,429 / 93.24%
  (Figure S1 uses the stored value and re-verifies it at run time).
- `REPORT_replicative_ceiling.md:234, 301` — IR vs OIS "+0.084"; stored value is +0.0835 → +0.083.
- `REPORT_replicative_ceiling.md:268` — CDKN2A distal BJ "+0.22"; stored value is +0.2149 →
  +0.215 (Figure 4 plots the stored value).
- `REPORT_ois_gates.md` §1.4 prose — "every condition peaks at K = 12"; ETS peaks at K = 6
  (Figure 6A plots and annotates the stored values).

# Data not available — panel reduced rather than reconstructed

**Figure 5B, AUC null.** The brief asked for the 200-permutation null distribution "for both
Spearman and AUC". `results_train.json` stores all 200 permutation **Spearman** values
(`control1_shuffled_labels.all_spearman`) but only summary statistics for AUC — mean, SD, min,
max and `n_beating`. The per-permutation AUCs were never written to disk. The AUC histogram is
therefore **not drawn**; the panel reports the stored AUC summary as text instead. Reconstructing
it would require re-running the 200 permutations on a GPU, which this task excludes.

# Files

| File | Purpose |
|---|---|
| `fig_style.py` | shared style, validated palette, size constants, save helper |
| `dataset_table.csv` | Figure 1 source table, assembled from the artifacts listed above |
| `figures_main.py` | Figures 1–7 |
| `figures_supp.py` | Figures S1–S5 |
| `figure_flags.txt` | mismatches detected at figure-build time |
| `figure*.png` / `figure*.pdf` | 300 dpi raster + vector, one pair per figure |
