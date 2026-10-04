# OIS-scoped enhancer training set — H3K27ac primary target

**Date:** 2026-08-25 · local, no GPU · ~1.1 GB new downloads · ~25 min end to end
**Product:** `ois_enhancer_trainset.h5ad` — **3 tasks × 111,671 intervals**

---

## Headline

| | |
|---|---|
| Consensus H3K27ac regions | **112,099** (summit supported by ≥2 of 4 GM21 samples) |
| With a valid 524 kb hg38 window | **111,677** → **111,671** after dropping ambiguous lifts |
| **Distinct 524 kb windows needed** | **8,432 — a 13.2× compression** |
| Caching cost (fwd + RC, measured rate) | **0.95 h** instead of 12.5 h naive — **saves 11.6 h** |
| Cross-study OIS agreement reproduced | **+0.637** vs prior **+0.653** (Δ −0.016) |
| Split | train GM21 (skin) → validate IMR90 (lung), chr9+chr6 held out |
| Third dataset (PRJNA439280) | **EXCLUDE** — see §2 |

**The efficiency question has a large answer: 111,677 regions collapse to 8,432 distinct
trunk windows.** Borzoi predicts a 196,608 bp span from each 524,288 bp input, and enhancers
cluster, so one forward pass serves ~13 regions on average. This is measured by greedy
interval covering per chromosome, not estimated.

---

## 1. Verification and acquisition

All accessions resolve live: GSE206402 (uid 200206402, n=25), GSE205898 (200205898, n=30),
GSE206496 (200206496, n=101), GSE74238 (200074238, n=22), GSE74328 (200074328, n=94),
PRJNA439280 (BioProject 439280).

| Source | Assay | Deposited | Build | Used |
|---|---|---|---|---|
| GSE205898 | H3K27ac ChIP | **narrowPeak + bigWig** | hg19 | 4 K27ac bigWigs + 4 narrowPeaks |
| GSE74238 | H3K27ac ChIP | **bigWig only — no peaks** | hg19 | 6 K27ac bigWigs (PRO/SEN/QUI ×2) |
| GSE206402 | ATAC | narrowPeak + bigWig | hg19 | 8 ATAC bigWigs + 25 narrowPeaks |
| PRJNA439280 | ATAC + 4 histone marks | **raw reads only** | — | **excluded** |

**Every file used is hg19**, so all signal is quantified in hg19 and no coverage value is
ever read through a lifted coordinate. Liftover hg19→hg38 is used only to emit
Borzoi-compatible intervals: **112,078 / 112,099 summits lifted (99.98%)**, 21 failed.

**GSE74238 deposits no peak calls at all** — only bigWigs. Region definition therefore comes
entirely from GSE205898, which is stated as a bias in §3 rather than glossed.

---

## 2. 🔴 Third dataset (PRJNA439280) — recommend EXCLUDE

**Decision: exclude. Do not process.** Four independent reasons, in order of severity:

1. **Disk is a hard blocker.** 58 SRA runs of ChIP + ATAC is roughly **150–200 GB** of
   FASTQ. **7.5 GB free** on this machine. It cannot be staged, let alone processed.
2. **The Zenodo route does not rescue it.** Record 3731264 *does* exist (DOI
   10.5281/zenodo.3731264, "AP-1 Imprints a Reversible Transcriptional Programme of
   Senescent Cells") — the survey's report was accurate. But its 14 files are 13
   supplementary tables plus **`WORKFLOW.tar.gz` at 58.4 GB**. There are no standalone
   bigWigs or peak files. 58 GB is also infeasible here.
3. **The cheap tables don't contain what's needed.** I downloaded and inspected
   `TABLE_S2.xlsx` (34 MB) — it holds H3K27ac/H3K4me1/H3K4me3/H3K27me3/ATAC **peak calls**
   at OIS T0 / 72h / 144h with columns `Chr, Start, End, Width, Position, Distance, Gene`.
   **Coordinates and annotation only — no signal values, no fold changes.** A
   presence/absence pseudo-response from peak calls would not be comparable to the
   bigWig-coverage log2FCs used for GM21 and IMR90, so it would degrade rather than extend
   the region set. `TABLE_S6` is siRNA expression data, not chromatin.
4. **Processing cost if pursued anyway:** ~175 GB download, an hg38 aligner index, then 58
   alignments at ~30–60 min each on 8 cores → **30–60 h of local compute**, plus dedup,
   coverage and peak calling. On a many-core cloud box this is a few hours of wall time but
   real money and a day of setup.

**What excluding it costs:** a third independent OIS dataset, and specifically a WI-38
(third tissue) replicate that would let the ceiling be checked three ways rather than two.
**What it does not cost:** the split. The train/validate design needs two independent
studies and has them (GM21 skin, IMR90 lung).

**If you want it later, the sensible order is:** rent a cloud box with disk, process the
H3K27ac libraries only (~12–16 of the 58 runs), and use it as a *third held-out test set*
rather than as training data. I did not start any of this and am not proceeding without
your go-ahead.

---

## 3. Region definition

**Source:** GSE205898 H3K27ac narrowPeaks. Summits (narrowPeak column 10 offset) clustered
within **250 bp**; a cluster is kept if supported by **≥2 of the 4 samples**; the region is
**1,000 bp centred on the cluster's median summit**.

| Sample (n=2 per arm) | Peak summits |
|---|---|
| EV Rep1 | 80,382 |
| EV Rep2 | 99,615 |
| RAS_D18 Rep1 | 132,065 |
| RAS_D18 Rep2 | 119,821 |
| **Consensus (≥2 of 4)** | **112,099** |

IMR90 signal coverage over those regions: 109,321–111,240 of 112,099 per sample (97.5–99.2%),
so the GM21-defined region set is well populated in the validation study.

### 🔴 Width, and how it relates to the 200 bp design target

**These are not the same coordinate, and conflating them would put designed sequence on a
nucleosome.** H3K27ac marks the nucleosomes *flanking* a regulatory element and is depleted
over the nucleosome-free region where transcription factors actually bind. So:

- the **H3K27ac summit** is where the acetylation signal is maximal — the right place to
  **measure**;
- the **nucleosome-free region**, located by the **ATAC summit**, is the right place to
  **design** 200 bp.

Measured directly: the **median distance from a K27ac summit to the nearest GM21 ATAC summit
is 245 bp**, and 92,616 / 112,099 regions (82.6%) have an ATAC summit inside the 1 kb window.
The offset is real and roughly one nucleosome-plus-linker, exactly as the biology predicts.

The build therefore stores **both**: `region_center_hg38` (K27ac summit — the measurement
anchor, what `.X` is computed over) and `design_anchor_hg38` (nearest ATAC summit — where a
200 bp element should be placed), with `design_anchor_dist`. The 1 kb measurement width is
chosen to span the flanking acetylated nucleosomes; it is deliberately *not* 200 bp, because
a 200 bp window centred on a K27ac summit would sit in the acetylation trough.

---

## 4. Targets

Per-region log2 fold change, computed **within each study independently** from mean-1
normalized bigWig coverage over the 1 kb windows, pseudocount 0.10.

| Task | Replicates | n regions | mean | sd | range |
|---|---|---|---|---|---|
| **GM21_RAS_vs_EV** *(train)* | **2 RAS vs 2 EV** | 112,099 | +0.057 | 0.761 | −3.26 … +3.29 |
| **IMR90_SEN_vs_PRO** *(validate)* | **2 SEN vs 2 PRO** | 111,063 | +0.128 | 0.965 | −3.18 … +4.84 |
| **IMR90_SEN_vs_QUI** *(validate)* | **2 SEN vs 2 QUI** | 111,066 | +0.117 | 1.067 | −3.88 … +5.19 |
| GM21_ATAC_RAS_vs_EV *(filter only)* | 4 RAS vs 4 EV | 112,099 | +0.176 | 0.557 | −2.77 … +3.32 |

**Every arm is n=2 except the ATAC filter.** That is the binding limitation on this build and
it is repeated with each number rather than footnoted.

### The quiescent arm is highly informative — keep it

| Quantity | Value |
|---|---|
| Spearman(SEN_vs_PRO, SEN_vs_QUI) | **+0.880** |
| Spearman(SEN_vs_PRO, QUI_vs_PRO) | **+0.036** |
| QUI vs PRO | mean +0.011, **sd 0.460** (vs 0.965 for SEN vs PRO) |

Read together: **the senescence response is nearly the same whichever arrested control is
used (+0.880), while quiescence itself barely moves H3K27ac (sd 0.46, mean ≈ 0) and does not
resemble a weaker senescence (+0.036).** So the H3K27ac programme captured here is
senescence-specific, not generic growth arrest. This matches the earlier finding that a
trained model placed quiescent cells near proliferating rather than senescent.

Both IMR90 contrasts are retained as separate validation tasks.

---

## 5. Harmonization check — the +0.653 reproduces

| Threshold | n | Spearman | Pearson |
|---|---|---|---|
| \|lfc\| ≥ 0 | 111,063 | +0.620 | +0.645 |
| ≥ 0.25 | 103,276 | +0.624 | +0.646 |
| **≥ 0.5 (primary)** | 83,716 | **+0.637** | +0.656 |
| ≥ 1.0 | 39,576 | +0.657 | +0.701 |

**Prior test +0.653 → this build +0.637, difference −0.016.** The region definitions are not
identical (this build centres on H3K27ac summits and drops ATAC-only regions; the prior used
midpoint-centred anchors from a combined ATAC+K27ac consensus), so exact equality was not
expected. The cross-study OIS agreement reproduces at the same magnitude and is flat across
thresholds, as before. **Harmonization is consistent with the earlier test — no unresolved
discrepancy.**

---

## 6. Batch check, and why no pooled correction is applied

PCA on 108,174 complete-case regions × 10 H3K27ac samples:

| Space | PC | Var | R²(study) | R²(status) | Dominated by |
|---|---|---|---|---|---|
| raw pooled | PC1 | 43.8% | 0.448 | **0.777** | **STATUS** |
| raw pooled | PC2 | 25.9% | 0.483 | 0.338 | study |
| per-study centered | PC1 | 52.8% | 0.000 | **0.908** | STATUS |

**Senescence status dominates PC1 even on raw pooled signal** (R² 0.777 vs 0.448) — markedly
unlike the expression work, where raw PC1 was 99.7% study. Batch structure is present on PC2
but does not dominate the leading component.

### The tradeoff, stated

The split is cross-study **by design**, so study identity and validation fold are perfectly
confounded. Any correction *fitted across both studies* would (a) use validation-set
information at training time, and (b) remove precisely the between-study variation the
holdout exists to test.

**Recommendation: do not apply a pooled batch correction.** The targets are already
per-study log2 fold changes — within-study contrasts in which the study-level offset cancels
by construction. That is the correct correction here: computed inside each study
independently, never mixing them. The pooled PCA above is reported only to show what pooling
raw coverage *would* have done.

---

## 7. Accessibility filter — optional, default OFF

| Threshold | Differential H3K27ac | Also differential ATAC | Of those, same direction |
|---|---|---|---|
| \|log2FC\| ≥ 0.5 | 56,949 | **25,002 (43.9%)** | **93.0%** |
| \|log2FC\| ≥ 1.0 | 20,695 | 5,890 (28.5%) | 100.0% |

Direction agreement of **93.0%** matches the prior test's 93.3% — consistent. Note the
asymmetry: only 43.9% of differential-H3K27ac regions also move in accessibility, whereas the
prior test found 65% of differential-*accessibility* regions also move in H3K27ac. Acetylation
change without accessibility change is the commoner event, which is a further argument for
H3K27ac as the primary target.

**Region counts both ways:** filter OFF → **111,671** intervals (default); filter ON
(`var['passes_atac_filter']`, |ATAC lfc| ≥ 0.5) → **39,202** intervals. Stored as a column so
a training run can subset without a rebuild.

---

## 8. The object

```
ois_enhancer_trainset.h5ad
  .X    (3, 111671) float32   signed log2 fold change, range [-3.879, +5.191]
        165,810 of 335,031 values are negative -> act_func=None is mandatory
  .obs  task, study, cell_line, tissue, mechanism, control_arm, n_replicates, role
  .var  chrom, start, end (524,288 bp hg38), region_center_hg38,
        region_start/end_hg38 (the 1 kb measurement window),
        summit_hg19, design_anchor_hg38, design_anchor_dist,
        n_samples_supporting, passes_atac_filter, split
  .uns  contract, target_transform, split_strategy, scope,
        accessibility_filter, window_covering, provenance
```

| Task | Study | Cell line / tissue | Control arm | n | Role |
|---|---|---|---|---|---|
| GM21_RAS_vs_EV | GSE205898 | GM21 / skin fibroblast | empty vector | 2v2 | **train** |
| IMR90_SEN_vs_PRO | GSE74238 | IMR90 / lung fibroblast | proliferating | 2v2 | validation |
| IMR90_SEN_vs_QUI | GSE74238 | IMR90 / lung fibroblast | quiescent | 2v2 | validation |

All interval widths are exactly 524,288. Splits: **train 100,070 / test 11,601** (chr9+chr6).

### Regions lost

| Reason | Count |
|---|---|
| Liftover failed | 21 |
| Window runs before chromosome start | 154 |
| **Window runs past chromosome end** | **247** |
| Ambiguous lift (distinct hg19 summits → one hg38 coordinate, segmental duplications) | 6 |
| **Total lost** | **428 of 112,099 (0.38%)** |

The 6 ambiguous lifts were dropped rather than arbitrarily keeping one of each pair.

### Recommended split direction and why

**Train on GM21, validate on IMR90.** Three reasons:

1. **Region definition already comes from GM21** (the only study with deposited peaks).
   Making GM21 the training source keeps the region-definition bias inside the training
   fold, where it is least harmful; using IMR90 to train would mean validating on regions
   defined by the validation study.
2. **GM21 has matched ATAC in the same cells**, which is what supplies the accessibility
   filter and the 200 bp design anchor — both training-time concerns.
3. **IMR90 carries the quiescent arm**, which is most valuable at validation, where the
   specificity check (does the model separate senescence from growth arrest?) actually
   bites.

---

## 9. Efficiency — the number that sets the GPU budget

| | |
|---|---|
| Regions with valid windows | 111,677 |
| **Distinct 524 kb windows required** | **8,432** |
| Compression | **13.2×** |
| Forward-only caching (0.2018 s/pass, measured on RTX 4090) | **0.47 h** |
| Forward + reverse-complement | **0.95 h** |
| Naive one-window-per-region (fwd+RC) | 12.5 h |
| **Saved** | **11.6 h ≈ $8.70 at $0.75/h** |

Computed by greedy interval covering per chromosome over the predicted 196,608 bp span, not
estimated. Per-chromosome breakdown in `window_covering_per_chrom.csv`.

*(The covering was computed over the 111,677 valid regions; the 6 ambiguous-lift regions were
dropped afterwards, which cannot increase the window count.)*

---

## 10. Sanity checks — all pass

### CDKN2A (chr9, 59 regions: 8 promoter, 51 distal) — the enhancer premise holds

| | GM21_RAS_vs_EV | IMR90_SEN_vs_PRO | IMR90_SEN_vs_QUI |
|---|---|---|---|
| promoter (≤2 kb) | **−0.919** | **−1.341** | **−1.419** |
| **distal** | **+0.172** | **+0.421** | **+0.150** |

Promoter down, distal up, **in both independent studies and against both control arms** —
reproducing the −1.03 / −1.07 promoter and +0.24 / +0.49 distal pattern that motivated the
enhancer pivot.

### CDKN1A (chr6, 37 regions: 9 promoter, 28 distal)

| | GM21 | IMR90 vs PRO | IMR90 vs QUI |
|---|---|---|---|
| promoter | +0.219 | −0.042 | −0.264 |
| distal | +0.100 | +0.198 | +0.247 |

Weaker and less consistent than CDKN2A — the promoter does not fall here — but distal is
positive in all three tasks.

### SASP loci — all up, in all three tasks

| Gene | n regions | GM21 mean (max) | IMR90 vs PRO mean (max) | IMR90 vs QUI mean (max) |
|---|---|---|---|---|
| IL6 | 12 | +0.280 (+1.115) | +0.238 (+1.331) | +0.200 (+1.330) |
| **CXCL8** | 17 | **+1.093** (+2.144) | **+2.528** (+3.823) | **+2.552** (+3.430) |
| **IL1A** | 34 | +0.718 (+2.682) | **+1.581** (+4.071) | **+1.587** (+4.241) |
| **MMP3** | 24 | **+1.081** (+1.974) | **+2.278** (+3.816) | **+2.272** (+4.079) |

Every SASP locus gains H3K27ac in every arm, and the gains survive the quiescent control —
so they are senescence-associated, not arrest-associated.

### Split integrity

- chr9 + chr6 → **11,601 regions marked `test`**, excluded from training.
- Training task: `GM21_RAS_vs_EV`. Validation tasks: `IMR90_SEN_vs_PRO`, `IMR90_SEN_vs_QUI`.
- **Sample-level overlap between studies: none** (disjoint GSMs, different labs).
- **On "no region in both sets":** the cross-study holdout is on **data source**, not on
  regions. The same regions are deliberately scored in both studies — that is the entire
  question being asked (does a model fit to GM21's response at an element predict IMR90's
  response at that same element?). Region-level overlap is therefore expected and is *not*
  leakage. The leakage hazard that does exist is that region **definition** came from GM21
  peaks; it is quantified in §3 (IMR90 coverage 97.5–99.2% of those regions).

---

## 11. Caveats, stated plainly

1. **Every arm is n=2.** Two replicates per condition in both studies. This bounds every
   number here and cannot be fixed with deposited data.
2. **Region definition is GM21-only**, because GSE74238 deposited no peaks. An IMR90-specific
   enhancer with no GM21 H3K27ac peak is invisible to this build.
3. **The GM21 senescent arm is RAS D18 only.** Later timepoints were excluded because this
   SuperSeries' subject is *escape* from OIS and D32+ showed an attenuated signature in
   earlier work.
4. **Two OIS studies, not three** — see §2.
5. **`.X` is signed.** A softplus head will silently floor the negative half; `act_func=None`
   is recorded in `uns['contract']['head']` and is mandatory.

---

## 12. Outputs

All at **`/Users/gabeykim/Downloads/Senescence/output/ois_enhancer_trainset/`**

| File | Contents |
|---|---|
| `ois_enhancer_trainset.h5ad` | the training set, 3 × 111,671 |
| `region_responses_ois.csv` | per-region responses + coordinates + filter flag |
| `regions_hg19_hg38.csv` | region definitions, both builds, ATAC design anchors |
| `signal_matrix_ois.csv` | raw per-region coverage, all 18 bigWig samples |
| `window_covering_per_chrom.csv` | the 8,432-window covering, per chromosome |
| `batch_pca_stats.csv` | PCA variance / R² decomposition |
| `locus_sanity_checks.csv`, `sasp_checks.csv` | CDKN2A/CDKN1A and SASP |
| `region_stats.json`, `step2_stats.json` | counts, liftover rates |
| `STEP1_LOG.txt`, `STEP2_LOG.txt` | full run logs |
| `figures/ois_batch_pca.png` | batch structure before/after centering |
| `scripts/70_build_ois_regions.py`, `scripts/71_ois_trainset.py` | build scripts |
