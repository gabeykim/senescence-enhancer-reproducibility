# Audit: Five senescence chromatin/expression datasets for a Borzoi/gReLU activity-head fine-tune

**Date:** 2026-08-20
**Scope:** GEO metadata (SOFT/series-matrix) and processed supplementary files, plus the ENCODE
REST/JSON API, only. No SRA/FASTQ access anywhere in this audit.
**Primary target (deep audit):** GSE175533 (Calico Hayflick-limit multi-omics, WI-38).
**Secondary targets (lighter audit, full detail in linked sections):** GSE206402 + SuperSeries
family, GSE210285, GSE74238/GSE74324, ENCSR200OML/ENCSR978WIX.
**Reproducibility:** every number in this report comes from a script under `scripts/`, run in
numeric-prefix order; raw outputs are under `data/` and `output/`.

---

## Headline findings (read these first)

1. **Every accession named in this brief is real and live.** Unlike the prior audit in this
   project (which found GSE254358 does not exist), all 8 GEO accessions and both ENCODE
   accessions resolve, confirmed via NCBI `esummary` and the ENCODE REST API, not by title
   inspection alone. Full log: [`output/tables/accession_verification_log.txt`](output/tables/accession_verification_log.txt),
   script [`scripts/06_verify_accessions.sh`](scripts/06_verify_accessions.sh). One accession
   count discrepancy was caught along the way: GSE206402's SuperSeries (GSE206496) actually
   contains a **fourth, previously-unlisted subseries** (GSE205692, a 34-sample expression
   microarray) that the task brief didn't mention — see the GSE206402 section.

2. **GSE175533's ATAC-seq and RNA-seq arms are NOT split-sample pairs.** At the 6 PDL points
   where both assays exist, the RNA and ATAC samples carry different BioSample accessions and
   materially different harvest protocols (Trizol lysis on-plate vs. trypsinize-and-count
   100,000 cells). They come from the same continuous WI-38 culture lineage staged to the same
   nominal PDL, not from one dish split into two assays. Treat any ATAC↔RNA integration as
   "matched condition, unmatched aliquot."

3. **Radiation-induced senescence (RIS) — the closest available proxy for the project's eventual
   doxorubicin/therapy-induced validation model — has ZERO ATAC-seq data.** Confirmed from both
   the per-sample `Sample_library_strategy` field (100% RNA-Seq for all 18 RIS samples) and the
   series' own design text. The cell-density (CD) growth-arrest control is likewise RNA-seq
   only. Neither can be used to separate senescence-specific accessibility change from
   generic-arrest accessibility change — that separation is only possible at the transcript
   level in this dataset.

4. **A true per-replicate ATAC PCA/correlation heatmap could not be produced.** GEO deposits
   per-sample ATAC signal only as bigWigs bundled in a 14.8GB tar; individual bigWigs are not
   separately downloadable. Measured sustained throughput to NCBI's FTP server was well under
   1.7 MB/s during this audit, making even the smallest useful subset (the ~8.2GB WI-38
   replicative-senescence arm alone) an impractical multi-hour download for this session. This
   is substituted with (a) the real per-replicate RNA-seq PCA, which strongly separates by PDL,
   and (b) the authors' own limma-based per-peak differential statistics (already computed and
   deposited), used for the promoter-locus check, the differential-region count, and the BED
   export. Flagged prominently rather than silently worked around — see Task 4.

5. **Neither CDKN2A nor CDKN1A promoters show any significant accessibility change with PDL** in
   GSE175533's ATAC data (all `adjPval > 0.05` at every timepoint, for every promoter-annotated
   peak near either gene) — despite CDKN1A being very significantly induced at the RNA level
   (`padj = 7×10⁻³⁴`) and genome-wide ATAC accessibility changing massively elsewhere
   (119,823 significant peaks by PDL50). This is a real, reproducible disconnect between
   RNA induction and promoter-level chromatin opening at the two flagship senescence
   cell-cycle-arrest genes, flagged prominently per the task's own instruction.

6. **The training-region count is NOT the bottleneck.** 119,823 peaks are significantly
   differential (padj<0.05, |log2FC|≥1) in the deepest RS-PDL-vs-hTERT comparison alone; 128,866
   in the pooled late-PDL set exported to BED. This comfortably clears the "few thousand region"
   concern in the task brief — the opposite problem from what was worried about.

7. **The three ATAC-bearing series in this audit (GSE175533, GSE206402, GSE210285) use two
   different genome builds and three incompatible levels of peak-file availability** (a merged
   atlas + differential stats for GSE175533; per-replicate MACS2 calls with no merged atlas for
   GSE206402; zero peak files, bigWigs only, for GSE210285). Harmonizing them into one consensus
   peak set is not a quick join — see Task 5.

---

## PART 1 — PRIMARY TARGET: GSE175533 (WI-38 Hayflick-limit multi-omics)

Confirmed live: accession=GSE175533, 147 samples, BioProject PRJNA732700, public 2021/05/26.
Source: Chan M et al., *eLife* 2022;11:e70283.

### Task 1 — Existence and metadata

SOFT + series-matrix pulled from `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175533/`
(script [`10_fetch_gse175533.sh`](scripts/10_fetch_gse175533.sh)) and parsed directly, no
GEOparse, no HTML scraping (script [`11_parse_gse175533_metadata.py`](scripts/11_parse_gse175533_metadata.py)).
Output: [`output/tables/gse175533_sample_metadata.csv`](output/tables/gse175533_sample_metadata.csv),
147 rows, 0 unparsed.

**Full sample inventory (arm × assay):**

| Arm | Assay | n | What it is |
|---|---|---|---|
| RS_atac | ATAC-seq | 19 | WI-38 replicative senescence, 7 PDL points |
| RS_atac_hTERT | ATAC-seq | 16 | hTERT-immortalized counterpart, matched timepoints |
| RS_bulkRNA | RNA-Seq | 30 | WI-38 replicative senescence, 10 PDL points |
| hTERT_bulkRNA | RNA-Seq | 18 | hTERT-immortalized counterpart, 6 timepoints |
| RIS | RNA-Seq | 18 | Radiation-induced senescence (10 Gy X-ray) + mock-IR control |
| CD | RNA-Seq | 30 | Cell-density / contact-inhibition growth-arrest control |
| RS_scRNA | RNA-Seq | 11 | Single-cell RNA-seq (not used in this audit's bulk QC) |
| hTERT_scRNA | RNA-Seq | 5 | Single-cell RNA-seq, hTERT (not used) |
| **Total** | | **147** | |

**Extraction fields:** assay type (from `Sample_library_strategy`, confirmed independently of
title text — agrees 100%), cell line (WI-38 vs. WI-38 hTERT), PDL (from
`Sample_characteristics_ch1: population doublings` for RS/hTERT arms; genuinely absent for
RIS/CD, which are staged by day instead — not a parsing failure), timepoint code, treatment,
replicate letter.

**Parsing trap caught (same pattern as the sibling GSE220545 audit elsewhere in this project):**
`Sample_extract_protocol_ch1` and `Sample_growth_protocol_ch1` contain one shared paragraph
describing the RS/hTERT, RIS, *and* CD protocols concatenated identically on every single
RNA-seq sample, regardless of which arm it actually belongs to. Per-sample condition was never
inferred from this text — only from the structured `Sample_title` prefix, which is fully
regular and machine-parseable (0 of 147 titles failed to match a known pattern).

**Two internal naming inconsistencies flagged (not silently resolved):**
- The deposited `atac_peaks_sig.xlsx` labels one ATAC comparison "PDL45_v_htert2" while GEO's own
  sample titles for that exact TP6 timepoint say "PDL46_TP6" — treated as the same condition
  based on matching TP-index, not proven via a shared identifier.
- The deposited TPM table's RIS columns are headed "d0_A/B/C" while GEO's sample titles for the
  same 3 samples say "RIS_d2_no_xray" — almost certainly the same mock-IR control samples under
  two different day-numbering conventions (absolute culture day vs. day-relative-to-IR), but not
  independently confirmed via a shared sample ID.
- The RS_TC TPM table's deep-senescence-arm columns (PDL45/52/53) drop the "_TP7/9/10" suffix
  that GEO's sample titles carry for the same samples — cosmetic, but breaks naive exact-string
  joins between the two files.

### Task 2 — Design adequacy

Full log: [`output/tables/gse175533_design_adequacy_log.txt`](output/tables/gse175533_design_adequacy_log.txt)
(script [`17_design_adequacy_gse175533.py`](scripts/17_design_adequacy_gse175533.py)).

**ATAC inventory:**

| Arm | PDL | n reps | Flag |
|---|---|---|---|
| RS_atac | 20, 25, 30, 33, 37 | 3 each | OK |
| RS_atac | 46, 50 | 2 each | **<3 reps** |
| RS_atac_hTERT | 46, 51, 64, 73 | 3 each | OK |
| RS_atac_hTERT | 86, 109 | 2 each | **<3 reps** |

**RNA inventory:** every RS_bulkRNA (10 PDL points), hTERT_bulkRNA (6 TPs), RIS (6 day-points),
and CD (10 day-points) arm has exactly 3 replicates — **no replicate-count problem anywhere in
the bulk RNA data.** The replicate shortfall is confined to the two latest ATAC timepoints in
each of the RS and hTERT arms.

**(a) PDL points with BOTH ATAC and RNA, and are they the same biological sample?**
ATAC PDLs: {20, 25, 30, 33, 37, 46, 50}. RNA PDLs: {20, 25, 28, 33, 37, 45, 46, 50, 52, 53}.
**Shared: {20, 25, 33, 37, 46, 50} — 6 PDL points.** Important gotcha: the "TP" label does *not*
track the same PDL between assays past TP2 (ATAC TP3 = PDL30 but RNA TP3 = PDL28; ATAC TP7 =
PDL50 but RNA TP7 = PDL45) — any join must be done on the numeric PDL field, never on TP label.
**Same biological sample or parallel cultures? Parallel cultures** — see headline finding #2.

**(b) Recommended PROLIFERATING/SENESCENT cutoff, justified from the data:**
Marker trajectories (full table: [`output/tables/gse175533_rna_markers_RS.csv`](output/tables/gse175533_rna_markers_RS.csv))
and the ATAC differential-peak-count escalation (below) both show the same sharp inflection
between **PDL46 and PDL50**, not a smooth gradient across the whole PDL20–53 range:

| PDL | MKI67 (TPM) | LMNB1 (TPM) | IL6 (TPM) | ATAC peaks sig. vs hTERT |
|---|---|---|---|---|
| 20 | 27.0 | 43.9 | 0.62 | — |
| 25 | 17.8 | 28.2 | 0.67 | 6 |
| 33 | 11.5 | 17.4 | 0.70 | 3,867 |
| 37 | 8.8 | 15.8 | 1.06 | 20,269 |
| 46 | 8.1 | 13.6 | 1.60 | 36,013 |
| **50** | **0.21** | **0.73** | **8.52** | **119,823** |
| 52 | 1.83 | 2.76 | 8.98 | — |
| 53 | 0.59 | 1.14 | 29.87 | — |

MKI67 and LMNB1 both crash by >10-fold specifically between PDL46 and PDL50; IL6 jumps >5-fold
in the same interval; the ATAC differential-peak count more than triples in the same interval
(the single largest jump of any adjacent pair). Two independent assays agree on where the
transition happens. **Recommendation: PDL20–33 as the proliferating (OFF) reference, PDL50–53 as
the senescent (ON) arm, and PDL37–46 excluded as a transitional zone** if a clean binary
contrast is needed for training labels. Note the ATAC time course itself only reaches PDL50 (no
ATAC at PDL52/53), and PDL50 ATAC has only 2 replicates — the ON-state ATAC signal rests on a
thin replicate base even though the statistical contrast (via limma, deposited pre-computed) is
strong.

**(c) Is the CD (cell-density) arm usable to separate senescence-specific from generic-arrest
signal?** Only at the RNA level. CD is 30 RNA-seq samples (10 timepoints × 3 reps), **RNA-seq
only — no ATAC-seq exists for CD anywhere in the series.** A sequence-to-accessibility model
cannot use this dataset to ask "is this accessibility change senescence-specific or just
generic-growth-arrest" — that comparison is unavailable at the chromatin level here.

**(d) Do RIS samples include ATAC, or RNA only?** **RNA only — confirmed.** All 18 RIS samples
have `Sample_library_strategy = RNA-Seq`; zero have `ATAC-seq`. The series' own
`Series_overall_design` text states it explicitly: "Further controls for slowed growth (cell
density) and DNA damage (radiation-RIS) were included for bulk RNA-seq" — not ATAC-seq. This
directly limits the project's stated interest in a doxorubicin/therapy-induced validation
proxy: **RIS/DNA-damage-induced senescence has no accessibility signal in this dataset at all.**

### Task 3 — Processed data QC

All four expected supplementary files exist exactly as predicted in the task brief, plus a bonus
scRNA h5ad object (not used here) — script [`10_fetch_gse175533.sh`](scripts/10_fetch_gse175533.sh).
**`GSE175533_RAW.tar` (14.8GB) was deliberately not downloaded** — see headline finding #4.

**Two of the three xlsx files use "Strict OOXML" conformance**, which openpyxl/pandas silently
misread as zero-sheet empty workbooks (no error raised) rather than failing loudly. Fixed with a
namespace-URI rewrite + rezip (script [`12_fix_strict_ooxml.sh`](scripts/12_fix_strict_ooxml.sh));
the fix was verified correct, not just assumed — see the coordinate-system check below.

**ATAC peak atlas** (`GSE175533_atlas.bed.gz` → `atlas.bed`, BED6, no header):
- **363,470 total peaks**, genome build **Hg38** (confirmed two independent ways: explicit
  `Sample_data_processing = Genome_build: Hg38` on every sample, and the `atac_peaks_sig.xlsx`
  readme sheet stating "region: hg38 coordinates").
- Width distribution: min 74bp, median 351bp, mean 421bp, p90 756bp, p99 1,342bp, max 3,191bp.
- Chromosome distribution: roughly proportional to chromosome size (chr1 largest at 34,537
  peaks, chrY smallest at 33 peaks) — no obvious chromosome-level artifact.
- Promoter vs. distal (from the `annot` column in `atac_peaks_sig.xlsx`, cross-validated against
  atlas.bed coordinates — see coordinate check below): **9.9% promoter (36,133), 90.1% distal**
  (intron 47.1%, intergenic-proximal 24.7%, intergenic-distal 18.6%, exon 2.7%).
- **Coordinate-system check performed, not assumed:** the `region` field in `atac_peaks_sig.xlsx`
  ("chr1:629876-630048") is exactly the atlas.bed interval ("chr1\t629875\t630048") plus 1 on the
  start — confirmed 1-based-inclusive `region` strings vs. 0-based BED, spot-checked at two
  independent loci (chr1 and the CDKN2A region on chr9) before using this conversion for the BED
  export in Task 5.

**RNA (TPM table, `GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx`):**
- 33,514 genes, Salmon-quantified TPM (no raw counts deposited for bulk RNA — "library size" in
  the literal read-count sense is not computable from what GEO provides; every sample's TPM
  column sums to ~1,000,000 by definition, confirmed as a sanity check).
- Detected genes (TPM>0) per RS-arm sample: 18,831–23,004 (56–69% of the 33,514-gene table),
  stable across PDL with a mild upward drift at the latest timepoints (52/53) — consistent with
  transcriptional activation of new gene sets in senescence, not a QC problem.
- Full per-sample table: [`output/tables/gse175533_rna_qc_RS.csv`](output/tables/gse175533_rna_qc_RS.csv).

### Task 4 — Signal check

**RNA-seq (genuine per-replicate PCA + correlation):**
Figure: [`output/figures/gse175533_rna_pca_correlation.png`](output/figures/gse175533_rna_pca_correlation.png).
Script: [`13_rna_qc_signal_gse175533.py`](scripts/13_rna_qc_signal_gse175533.py). log2(TPM+1)
(substituted for the requested CPM — only TPM is deposited; TPM and CPM are both per-million
normalizations, so this is a like-for-like adaptation, flagged as such, not a silent swap).

- **PC1 explains 49.8% of variance and correlates with PDL at r=0.858.** PC2 (17.9%) mostly
  isolates the PDL45/52/53 "deep senescence secondary time course" samples as a distinct batch
  (visible as the green outlier cluster in the figure) — plausibly a real batch effect from that
  sub-experiment being generated separately, per the growth-protocol text.
- Replicates cluster tightly within PDL (mean within-PDL Spearman ρ=0.9515 vs. mean
  across-all-samples ρ=0.9407) — small but real and directionally consistent with genuine
  biological structure, not replicate noise dominating.

**ATAC-seq (substituted signal check — see headline finding #4 for why no literal per-replicate
PCA was run):** the differential-peak-count escalation with PDL (Task 2b table above; figure
[`output/figures/gse175533_atac_differential_escalation.png`](output/figures/gse175533_atac_differential_escalation.png))
is used as the macro-level ATAC signal-quality check instead. The monotonic, dramatic scaling —
6 significant peaks at PDL25 to 119,823 at PDL50 — is strong evidence of a real, PDL-tracking
accessibility signal; it would be very difficult to produce this exact escalation pattern from
noise or sample mislabeling.

**RNA senescence marker panel** (mean TPM, PDL20 → PDL53; full trajectories:
[`output/tables/gse175533_rna_markers_RS.csv`](output/tables/gse175533_rna_markers_RS.csv)):

| Gene | PDL20 → PDL53 | log2FC | Direction | Expected | Match? |
|---|---|---|---|---|---|
| CDKN2A (p16) | 45.1 → 63.5 | +0.49 | UP | up | OK (but noisy, non-monotonic; Wald linear-trend padj=0.41, not significant) |
| CDKN1A (p21) | 530 → 909 | +0.78 | UP | up | OK, very significant (Wald padj=7×10⁻³⁴) |
| LMNB1 | 43.9 → 1.14 | −5.26 | DOWN | down | OK, clean monotonic decline, sharpest drop PDL46→50 |
| MKI67 | 27.0 → 0.59 | −5.49 | DOWN | down | OK, clean monotonic decline, sharpest drop PDL46→50 |
| IL6 | 0.62 → 29.87 | +5.57 | UP | up | OK, sharpest rise PDL46→50 |
| CXCL8 (IL8) | 0.27 → 8.26 | +4.89 | UP | up | OK, but the rise is concentrated PDL50→52, later than IL6 |
| **SERPINE1 (PAI-1)** | 633 → 454 | **−0.48** | **DOWN** | up (SASP) | **MISMATCH — flagged prominently** |

**SERPINE1 does not behave as expected.** Raw TPM shows a mild net decline PDL20→PDL53, and the
DESeq2 linear-time model finds essentially no trend at all (log2FC=−0.091, padj=0.826) — not
"significantly down," but clearly not the expected SASP-like induction either. The full
trajectory fluctuates (peaks mid-course around PDL25/46, dips at the very end) rather than
showing a clean rise. This is a real deviation from the canonical SASP-factor expectation and is
called out here rather than smoothed over.

**Accessibility at the CDKN2A / CDKN1A promoters (the direct test of whether ATAC signal carries
the biology the model needs to learn):** full detail
[`output/tables/gse175533_atac_peak_qc_diff_markers_log.txt`](output/tables/gse175533_atac_peak_qc_diff_markers_log.txt).

**Neither promoter shows a significant accessibility change at any PDL timepoint.** CDKN2A has 3
promoter-annotated peaks nearby; all 5 PDL-vs-hTERT comparisons for all 3 peaks have
`adjPval > 0.18` (mostly >0.3–0.9), with log2FC magnitudes under 0.6 and inconsistent sign across
timepoints. CDKN1A's single promoter peak tops out at log2FC=+0.231, adjPval=0.156 at PDL50 — the
closest to significance of any of these loci, still not significant. **This holds despite CDKN1A
being one of the RNA data's most significantly-induced genes (padj=7×10⁻³⁴) and despite the
genome-wide accessibility landscape changing enormously (119,823 significant peaks elsewhere by
PDL50).** Two non-exclusive interpretations, neither resolvable from processed peak data alone:
(1) p16/p21 transcriptional induction in this system may not be driven primarily by
promoter-accessibility opening — regulation could act through other chromatin marks, TF binding
at distal elements, or post-transcriptional mechanisms; or (2) mosaic/heterogeneous induction
across the cell population (well documented for p16 in replicative senescence specifically) could
dilute a real per-cell promoter-opening signal below detection in bulk ATAC, the same way it can
dilute bulk RNA. **Practical implication: a model trained on this dataset's genome-wide
accessibility signal is not being handed a clean, learnable CDKN2A/CDKN1A promoter-opening
example** — whatever the model learns about senescence-associated accessibility, it is not
learning it at these two specific, most-canonical loci.

### Task 5 — Training-set construction feasibility

Full detail: [`output/tables/gse175533_atac_differential_summary.csv`](output/tables/gse175533_atac_differential_summary.csv),
script [`14_atac_peak_qc_diff_markers.py`](scripts/14_atac_peak_qc_diff_markers.py).

**Differential accessibility, using the authors' own deposited limma statistics** (RS-PDL vs.
matched-hTERT-timepoint, quantile-normalized peak counts — this is the only differential ATAC
signal GEO provides; per-replicate re-derivation was not feasible, see headline finding #4).
Threshold: `adjPval<0.05` and `|log2FC|≥1` (conventional, stated explicitly, not the only
defensible choice):

| Comparison | Peaks tested | Significant | Gained | Lost |
|---|---|---|---|---|
| PDL25 vs hTERT-TP2 | 363,470 | 6 | 4 | 2 |
| PDL33 vs hTERT-TP4 | 363,470 | 3,867 | 3,051 | 816 |
| PDL37 vs hTERT-TP5 | 363,470 | 20,269 | 10,536 | 9,733 |
| PDL46 vs hTERT-TP6 | 363,470 | 36,013 | 17,588 | 18,425 |
| PDL50 vs hTERT-TP7 | 363,470 | 119,823 | 58,131 | 61,692 |

**Union across all 5 comparisons: 134,774 distinct peaks.** Using the two "late" (senescent-proxy)
comparisons specifically, as recommended by the Task 2(b) cutoff: **128,866 peaks**, split
roughly evenly (62,168 gained / 66,698 lost accessibility in senescence) — a healthy,
non-degenerate direction balance for training a bidirectional accessibility-change predictor.
**This is two orders of magnitude above the "few thousand regions" underpowered threshold flagged
in the task brief — training-region count is not a limiting factor for this dataset.**
BED export (128,866 regions, hg38, direction and log2FC in the name field):
[`output/tables/gse175533_differential_accessibility_senescent_vs_hTERT.bed`](output/tables/gse175533_differential_accessibility_senescent_vs_hTERT.bed).

**Cross-dataset compatibility (GSE175533 vs. GSE206402 vs. GSE210285):**

| | GSE175533 (WI-38, RS) | GSE206402 (GM21, OIS) | GSE210285 (2BS, RS) |
|---|---|---|---|
| Genome build | **hg38** | **hg19/GRCh37** | **hg38** |
| Peak caller | MACS2, `-p 1e-1 --nomodel --shift -37 --extsize 73`, then IDR (0.05) across replicate pairs | MACS2 v2.2.7.1, `-p 1e-3 --nomodel`, then IDR → "master peaksets" (described, not deposited) | MACS2 (version unresolved; garbled processing text), params largely unrecoverable |
| Merged/atlas peak file deposited? | **Yes** — 363,470-peak atlas + differential stats | **No** — only 25 per-replicate narrowPeak files (75,969–178,175 peaks each) | **No peak file of any kind** — only 4 bigWig tracks |
| Differential stats deposited? | Yes (limma, 5 PDL-vs-hTERT comparisons) | No | No |

**Concretely, what harmonization would require:**
1. **Liftover GSE206402 hg19→hg38** to match the other two — not optional, a straight coordinate
   join across builds is wrong.
2. **Build a merged/consensus peak set for GSE206402** from its 25 per-replicate files (currently
   none exists despite the processing notes describing one) — e.g. re-running the IDR step the
   authors described but never deposited.
3. **Call peaks from scratch for GSE210285** — it contributes zero peak files, only signal tracks;
   this is the most work of the three.
4. **Reconcile differing peak-calling stringency even after that.** GSE175533 used a very
   permissive MACS2 p-value cutoff (1e-1) before IDR filtering; GSE206402 used a much stricter
   cutoff (1e-3) with no IDR-merged deposit. A naive overlap/Jaccard comparison of GEO-deposited
   peaks as-is would conflate genuine biological differences with pipeline-threshold artifacts.
5. The fully rigorous path is **re-processing all three from raw reads with one unified
   pipeline** (peak caller, parameters, genome build) — which requires SRA access, explicitly out
   of scope for this audit. A lighter interim path (liftover + re-call GSE210285's peaks with
   GSE175533-matched MACS2 parameters, then IDR-merge GSE206402's replicates) is feasible without
   raw reads but still nontrivial engineering, not a quick script.

---

## PART 2 — SECONDARY TARGETS (condensed; full detail in linked sections)

### GSE206402 + SuperSeries family (GM21, OIS) — [full section](output/gse206402_section.md)

25 ATAC-seq samples (empty-vector control at 3 timepoints, H-RAS-G12V overexpression at 7),
**7 of 10 condition×timepoint arms have only n=2 replicates** (only day-8 both arms and RAS-day-23
clear n≥3). Genome build **hg19**. Peak files are per-replicate MACS2 narrowPeak only — no
merged atlas exists despite processing notes describing one. **ATAC and RNA-seq (GSE206493) are
NOT paired** — the RNA arm has no empty-vector control at all (it's a POU2F2-knockdown
perturbation panel) and only 1 of its 2 timepoints overlaps the ATAC time course. The SuperSeries
actually has 4 assay subseries, not the 3 named in the task brief (a 34-sample microarray
subseries, GSE205692, was found via sample-count reconciliation). Every subseries has its own
separate BioProject — a naive "same BioProject ⇒ paired" heuristic fails completely here.

### GSE210285 (2BS, replicative senescence) — [full section](output/gse210285_section.md)

**"2 replicates per arm" confirmed exactly**: 4 samples total, `phenotype: Growing` vs.
`phenotype: Senescence`, 2 reps each — stated explicitly in `Series_overall_design`. Replicate
identity is confounded with sequencing platform (rep1=Illumina HiSeq 4000, rep2=BGISEQ-500 in
both arms), so platform and biological variance can't be statistically separated at this depth.
Genome build **hg38**, unambiguous. **No peak file exists in GEO for this series at all** —
despite processing notes describing a MACS2 call, only 4 bigWig signal tracks are deposited.
Passage/PDL is not reported anywhere. Verdict: usable only as a small, low-weight, single-context
addition once peaks are re-called from the bigWigs — not standalone.

### GSE74238 (ChIP) + GSE74324 (RNA) (IMR90, OIS) — [full section](output/gse74238_74324_section.md)

GSE74238 (22 samples): H3K27ac, BRD4, and input each have **exactly n=2 replicates** per
proliferating/quiescent/senescent condition — below the ≥3 bar, flagged explicitly. GSE74324 (72
samples) is **not** a simple 3-condition design — it's 8 perturbation arms (shRNA
knockdowns + JQ1) × P/Q/S × 3 replicates. Only the **shRen ("no treatment") arm — 9/72 samples —
is a clean, unperturbed P/Q/S baseline** directly comparable to the ChIP data; the other 75% of
samples are a mechanistic dependency panel, not a senescence-vs-proliferating resource. Genome
build **hg19**, confirmed identically in both series. Both series are genuine SubSeries of one
SuperSeries (GSE74328) with matched core biology, but GSE74238 has no ChIP counterpart for any of
GSE74324's knockdown arms — treat as paired only for the ~9–18 overlapping P/Q/S samples, not all
94 combined. Notable: GSE74238's 14 ChIP bigWigs (not downloaded, 3.79GB combined) are exactly
the kind of continuous per-base track a Borzoi-style model would want, if retrieval effort is
later invested.

### ENCSR200OML + ENCSR978WIX (ENCODE, IMR-90) — [full section](output/encode_imr90_section.md)

Clean, well-characterized proliferating-baseline ATAC-seq for IMR-90: 2 biological replicates
(asymmetric — rep1 nucleus-isolated, 95.2% mapped, 30.2M usable fragments; rep2 unspecified
fraction, only 77.2% mapped, 22.0M usable fragments, ENCODE's own audit flags it for low depth).
Genome build **GRCh38**. Complete standard ENCODE4 processed-output set (IDR peaks, conservative
IDR peaks, overlap-based default peaks — 160,824 distinct peaks after de-duplication — signal and
fold-change bigWigs). ENCSR978WIX is a ChromBPNet model trained directly on ENCSR200OML's own
alignments (not independent data) — useful as a methodological reference (a real
sequence→accessibility model already fit on this exact cell line) but not as a second empirical
accessibility measurement. **CDKN2A's promoter shows a real, reproducible ATAC peak in this
proliferating baseline** (chr9:21,993,965–21,996,073, present in every ENCODE peak set and
independently recovered by ChromBPNet), flagged prominently — see Task 6 below for what this
means for cross-line comparison. No passage number is reported for either replicate.

---

## PART 3 — TASK 6: VERDICT

### Is GSE175533 sufficient as the primary training set?

**Yes, for a differential-accessibility contrast — with real, specific caveats, not a clean yes.**

**What's strong:**
- 128,866 differential regions (senescent-proxy vs. hTERT), balanced direction split — an order
  of magnitude more than needed to avoid an underpowered sequence model.
- 5 biologically-replicated PDL timepoints (3 reps each) covering the meaningful early-to-mid RS
  range, with a clear, two-assay-corroborated inflection point (PDL46→50) to justify a
  proliferating/senescent split.
- Genome build is explicit and consistent (hg38) across every sample and every supplementary file
  in the series.
- RNA-seq signal is clean, strongly PDL-correlated (PC1 r=0.858), and 6/7 canonical markers move
  as expected.

**What's weak, specifically:**
- The two ATAC timepoints that matter most for the senescent (ON) end of the contrast — PDL46 and
  PDL50 — have only **2 replicates each**, not 3+.
- ATAC and RNA are not matched-aliquot pairs (parallel cultures only) — a model wanting to learn a
  joint accessibility→expression relationship at the same physical sample doesn't have that here.
- **RIS (DNA-damage-induced senescence) and CD (growth-arrest control) have zero ATAC data** — the
  dataset cannot support the project's stated interest in a therapy-induced-senescence proxy at
  the accessibility level, and cannot separate senescence-specific from generic-arrest chromatin
  change at all.
- **CDKN2A and CDKN1A promoters — the two most canonical senescence loci — show no significant
  accessibility change anywhere in the PDL course**, despite massive genome-wide accessibility
  change elsewhere and (for CDKN1A) highly significant RNA induction. A model evaluated
  specifically at these loci should not be expected to show a strong senescence-associated signal
  from this training data.
- SERPINE1 (a canonical SASP factor) does not move as expected — a second reminder that not every
  textbook marker is well-represented in this bulk, replicative-senescence-specific dataset.
- No true per-replicate ATAC PCA/correlation was possible in this audit (bandwidth-constrained);
  the differential-region set rests on the authors' own limma model, not an independently
  reproduced one.

**Recommended train/validation/test split — by chromosome, not by sample, to avoid leakage:**

| Split | Chromosomes | Rationale | Peaks (of 363,470 atlas) |
|---|---|---|---|
| **Test** | chr8, chr9 | chr9 contains CDKN2A — holding it out lets a trained model be evaluated, unseen, at exactly the locus this whole audit centers on | 32,606 (9.0%) |
| **Validation** | chr6, chr16 | chr6 contains CDKN1A, for the same reason | 32,055 (8.8%) |
| **Train** | all remaining autosomes (chr1–5, chr7, chr10–15, chr17–22) + chrX | — | ~298,809 (82.2%) |

Holding out chr9 and chr6 specifically means the eventual model's performance at CDKN2A and
CDKN1A is a genuine generalization test, not a memorized training example — directly useful given
finding #5 above. chrY (33 peaks total) is negligible and can be dropped or folded into train
without materially affecting any split.

### Cell-line transfer risk: WI-38 (training) vs. IMR-90 (wet-lab validation)

**Cannot be directly resolved from data alone — GSE175533 (WI-38) and ENCSR200OML (IMR-90) were
never profiled in the same experiment, on the same platform, or with the same peak-calling
pipeline** (GSE175533: hg38, MACS2 `-p 1e-1` + IDR at 0.05; ENCSR200OML: GRCh38, ENCODE4's
standard IDR/overlap pipeline — same assembly name, independently processed). The one concrete,
comparable data point available: **CDKN2A promoter accessibility is present and reproducible in
proliferating IMR-90 (ENCODE) but shows no significant PDL-dependent change in WI-38 (GSE175533)
at any timepoint** — these are different questions (baseline presence vs. differential change
across senescence progression), so this is not a direct contradiction, but it does mean the two
cell lines cannot be assumed to behave identically at this specific, biologically central locus.
Neither cell line's ATAC record reports a usable passage number (WI-38: available as numeric PDL;
IMR-90/ENCODE: `passage_number: null` on both replicates), so even a passage-matched comparison
isn't possible. **Flag this as a real, unquantified transfer risk**: a model trained on WI-38
replicative-senescence accessibility patterns is not guaranteed to transfer to IMR-90 without
validation on IMR-90 data itself, and no data source in this audit provides a matched
IMR-90-replicative-senescence time course to check against — ENCODE's IMR-90 record is a single
proliferating baseline, not a senescence time course.

---

## Assumptions log (consolidated — see each dataset's own section for full detail)

1. **GSE175533** ATAC↔RNA replicate-letter correspondence at matched PDL is by title parsing only
   (BioSample accessions differ, confirmed not the same aliquot).
2. **GSE175533** "PDL45_v_htert6" (atac_peaks_sig.xlsx) is treated as the same condition as GEO's
   "PDL46_TP6" ATAC samples, matched by TP-index, not a shared identifier.
3. **GSE175533** TPM table's RIS "d0" columns assumed equivalent to GEO's "RIS_d2_no_xray" samples
   (same sample count, same structural position as the pre-treatment control) — not confirmed via
   shared ID.
4. **GSE175533** Task 4's "CPM-normalize" instruction was executed as log2(TPM+1) since only TPM
   (not raw counts) is deposited for bulk RNA.
5. **GSE175533** ATAC signal check substitutes the authors' own differential statistics and an
   escalation-pattern analysis for a literal per-replicate PCA, which was infeasible given
   measured sub-1.7MB/s throughput to NCBI FTP for the necessary multi-GB bigWig downloads.
6. **GSE175533** BED coordinate conversion (1-based `region` string → 0-based BED) was verified
   empirically against atlas.bed at two loci, not assumed from convention alone.
7. Every secondary dataset's assumptions are logged in full in its own linked section (parsing
   heuristics, replicate-inference logic, garbled processing-text reconstructions, etc.) —
   not repeated here.

---

## Reproducing this audit

```bash
# Existence gate (run first, always)
bash scripts/06_verify_accessions.sh

# Primary target: GSE175533
bash scripts/10_fetch_gse175533.sh
python3 scripts/11_parse_gse175533_metadata.py
bash scripts/12_fix_strict_ooxml.sh data/GSE175533/GSE175533_atac_peaks_sig.xlsx data/GSE175533/GSE175533_atac_peaks_sig.transitional.xlsx
bash scripts/12_fix_strict_ooxml.sh data/GSE175533/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx data/GSE175533/GSE175533_TPM_table.transitional.xlsx
python3 scripts/13_rna_qc_signal_gse175533.py
python3 scripts/14_atac_peak_qc_diff_markers.py
python3 scripts/17_design_adequacy_gse175533.py

# Secondary targets (each independently reproducible; see their own sections for exact commands)
bash scripts/20_fetch_gse206402_family.sh && python3 scripts/21_fetch_gse206402_peaks.py && python3 scripts/22_parse_gse206402_metadata.py && python3 scripts/23_design_summary_gse206402.py
bash scripts/21_fetch_gse210285.sh && python3 scripts/22_parse_gse210285.py && python3 scripts/23_design_adequacy_gse210285.py
bash scripts/22_fetch_gse74238_74324.sh && python3 scripts/23_parse_gse74238_74324.py && python3 scripts/24_design_summary_gse74238_74324.py
bash scripts/23_fetch_encode_imr90.sh && python3 scripts/23_analyze_encode_imr90.py
```

Requires `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, `openpyxl`, `curl`, `unzip`,
`zip` on `PATH`. No GEOparse dependency anywhere (all SOFT/matrix files parsed directly).
`pyBigWig` was installed mid-audit but ultimately unused for GSE175533 (see headline finding #4);
it remains available for anyone extending this audit with better bandwidth.
