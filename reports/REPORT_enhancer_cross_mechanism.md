# Cross-mechanism kill test on ENHANCER responses

**Date:** 2026-08-24 · run locally, no GPU · ~5.2 GB downloaded, ~35 min end to end

---

## Verdict

**Enhancer responses diverge across induction mechanisms at least as badly as gene
expression did — and the divergence is now measured against a high, well-established
within-mechanism ceiling, so it cannot be blamed on measurement noise.**

| Comparison | What it is | Spearman | n |
|---|---|---|---|
| **(b) GM21 OIS H3K27ac vs IMR90 OIS H3K27ac** | **within-mechanism CEILING** | **+0.653** | 88,530 |
| **(d) GM21 OIS ATAC vs GM21 OIS H3K27ac** | within-study, cross-assay | **+0.607** | 36,560 |
| **(a) WI-38 replicative ATAC vs GM21 OIS ATAC** | **cross-mechanism** | **−0.123** | 93,989 |
| **(c) WI-38 replicative ATAC vs GM21 OIS H3K27ac** | cross-mechanism, cross-assay | **−0.118** | 32,049 |

Swapping the **assay** (ATAC → H3K27ac) while holding mechanism fixed costs ~0.05 of
correlation. Swapping the **mechanism** while holding assay fixed costs ~0.78 and flips the
sign. **Mechanism is the axis that matters; assay is not.**

**A unified enhancer target is not defensible. The enhancer pivot must be scoped to a
single induction mechanism from the outset.**

---

## 1. What was usable, and what was not

| Series | Cell line / mechanism | Assay | Deposited | Used? |
|---|---|---|---|---|
| GSE175533 | WI-38, replicative | ATAC | peak atlas + **deposited per-region log2FC** (363,470 regions, hg38) | ✅ |
| GSE206402 | GM21, OIS | ATAC | 25 narrowPeak + 25 bigWig (hg19) | ✅ 8 bigWigs |
| GSE205898 | GM21, OIS | H3K27ac / H3K4me1 ChIP | narrowPeak + bigWig (hg19) | ✅ 4 K27ac bigWigs |
| GSE74238 | IMR90, OIS | H3K27ac / BRD4 ChIP | **bigWig only — no peaks ever deposited** | ✅ 4 K27ac bigWigs |
| GSE210285 | 2BS, replicative | ATAC | **4 bigWigs, ZERO peak files**, n=2/arm | ❌ **EXCLUDED** |

**GSE210285 is confirmed unusable and excluded**, as the brief anticipated. Its
`filelist.txt` contains one TAR and four bigWigs — no peak file of any kind — and n=2 per
arm. There is nothing to build a consensus region set from, and it is the *only* second
replicative dataset available. That has a consequence, stated in §7.

**Only per-sample files were downloaded** (5.2 GB), not the RAW archives — GEO serves
individual supplementary files under `/geo/samples/`, so GSE205898's 9.0 GB, GSE206402's
8.2 GB and GSE74238's 3.8 GB tars were all avoided.

---

## 2. Harmonization — and why liftover turned out not to be the problem

**All three OIS studies are hg19; only GSE175533 is hg38.** Unlike the expression work,
peak coordinates cannot sidestep this by being symbol-keyed. Two choices followed:

**Anchor space is hg19.** Every bigWig is hg19, so anchoring there keeps liftover out of
the measurement path entirely — it is used *only* to carry anchors into hg38 to look up
GSE175533's deposited response. Anchoring in hg38 would have put an hg38→hg19 lift in
front of every one of the 3.8 M coverage queries.

**Fixed-width 1 kb recentering** on peak midpoints. Peak widths differ by caller and
pipeline between these studies, so raw intervals are not comparable units; a fixed width
makes per-base coverage comparable and removes peak width as a confounder. It also means
only the *midpoint* needs lifting.

**Consensus criterion:** a segment is retained if called in **≥ 2 samples** of that assay
(sweep-line depth over per-sample merged intervals, so no single sample can self-support).
Anchors within 500 bp of each other are merged.

| | |
|---|---|
| GM21 ATAC consensus segments (≥2 of 8 samples) | 180,099 |
| GM21 H3K27ac consensus segments (≥2 of 4 samples) | 122,971 |
| **Anchors after 1 kb recentering + dedup** | **237,824** |
| — in GM21 ATAC consensus | 180,078 |
| — in GM21 H3K27ac consensus | 112,244 |
| — in **both** (supports comparison d) | 54,498 |
| — matched into a WI-38 peak after lifting | **136,846 (57.5%)** |

### Liftover was not lossy — quantified, not assumed

| Check | Result |
|---|---|
| Midpoints lifting hg19→hg38 | **237,628 / 237,824 = 99.92%** (196 failed, 0.08%) |
| Span integrity (both ends lifted, 5,000-anchor sample) | 4,996 / 5,000 = 99.92% |
| — of those, lifted span within 50 bp of 1 kb | **4,994 / 4,996 = 99.96%** |

The brief's caution is correct in general — the ~99.9% figure is for point positions, not
wide intervals — but at 1 kb width the span survives in 99.96% of cases. **Liftover does
not degrade this region set**, so the correlations below are not being reported on a
damaged set.

---

## 3. Response definitions and normalization

| Study | Response (log2, senescent ÷ proliferating) |
|---|---|
| **WI-38 replicative ATAC** | `PDL50_v_hTERT7 − mean(PDL25_v_hTERT2, PDL33_v_hTERT4, PDL37_v_hTERT5)`. Each arm is already a log2FC against its own matched hTERT control, so the subtraction cancels the hTERT baseline. **PDL45 excluded as transitional** (prior audit put the inflection at PDL46–50). |
| **GM21 OIS ATAC** | RAS D18 + D23 (n=4) vs pBABE D8 + D32 (n=4) |
| **GM21 OIS H3K27ac** | RAS D18 (n=2) vs EV (n=2) |
| **IMR90 OIS H3K27ac** | senescent (n=2) vs proliferating (n=2) |

**RAS D32+ excluded** from GM21 senescent arms: this SuperSeries' subject is *escape* from
OIS, and the expression build measured a markedly attenuated signature at D32–56.

**Normalization:** each bigWig sample's per-anchor mean coverage is divided by that
sample's mean across all anchors (every sample → mean 1). Pseudocount 0.10 added before the
log ratio so near-zero anchors cannot produce unbounded fold changes. This is a within-
sample depth normalization; it cannot fix genuine pipeline differences, which is why §4 runs
first.

**Stated asymmetry:** comparison (a) mixes a deposited DESeq-style log2FC (WI-38) with a
coverage ratio (GM21). There is no way around this without discarding the best available
WI-38 estimate. It is precisely why (b) and (d) are load-bearing rather than optional —
and why §6 tests the WI-38 construction directly.

---

## 4. Batch check — study dominates raw signal, and is fully removed by centering

PCA on the 206,604-anchor × 16-sample complete-case matrix:

| Space | PC | Var | R²(study) | R²(status) | Dominated by |
|---|---|---|---|---|---|
| raw pooled | PC1 | 47.5% | **0.969** | 0.022 | **STUDY** |
| raw pooled | PC2 | 16.8% | 0.093 | **0.784** | STATUS |
| raw pooled | PC3 | 12.5% | 0.879 | 0.052 | STUDY |
| **per-study centered** | **PC1** | **41.5%** | **0.000** | **0.878** | **STATUS** |
| per-study centered | PC2 | 19.1% | 0.000 | 0.042 | STATUS |

Raw PC1 is 96.9% explained by study — the same hazard as the expression work (which was
99.7%). **But it is cleanly removable here**: after per-study centering, R²(study) drops to
**0.000** and PC1 becomes 87.8% senescence status.

Importantly, **the correlations in §5 are computed on per-study log2 fold changes, which are
already within-study contrasts** — the study offset cancels in the ratio by construction.
The batch structure is therefore not what drives the results below. Note also that status
already dominates raw PC2 (78.4%), which it never did in the expression data.

---

## 5. The headline analysis

Primary threshold: |log2FC| ≥ 0.5 in at least one axis. Sensitivity across all thresholds:

| Comparison | \|lfc\|≥0 | ≥0.25 | **≥0.5 (primary)** | ≥1.0 |
|---|---|---|---|---|
| **(a)** WI-38 rep ATAC vs GM21 OIS ATAC | −0.124 (n=131,093) | −0.124 (n=119,811) | **−0.123 (n=93,989)** | −0.117 (n=46,977) |
| **(b)** GM21 OIS K27ac vs IMR90 OIS K27ac | +0.640 (n=111,204) | +0.643 (n=104,873) | **+0.653 (n=88,530)** | +0.678 (n=47,662) |
| **(c)** WI-38 rep ATAC vs GM21 OIS K27ac | −0.088 (n=43,499) | −0.098 (n=40,127) | **−0.118 (n=32,049)** | −0.115 (n=15,118) |
| **(d)** GM21 OIS ATAC vs GM21 OIS K27ac | +0.553 (n=54,498) | +0.566 (n=49,333) | **+0.607 (n=36,560)** | +0.702 (n=13,428) |

**Every comparison is flat across thresholds** — no threshold choice changes any conclusion.
Pearson tracks Spearman throughout (a: −0.108, b: +0.672, c: −0.112, d: +0.650 at primary).

### The ceiling is high, which is what makes the rest interpretable

Comparison **(b)** reaches **ρ = +0.653** between **different labs, different cell lines**
(GM21 skin vs IMR90 lung fibroblasts), **different ChIP protocols, and different peak
pipelines** — with one side (IMR90) having no deposited peaks at all and being quantified
purely from bigWig coverage. That is a demanding positive control, and it passes decisively.

**So the cross-mechanism values are not measurement noise.** The apparatus can detect
agreement at ρ ≈ 0.65 when the mechanism is shared. It reports ρ ≈ −0.12 when it is not.

### The cross-mechanism correlations are ~zero and mildly NEGATIVE

ρ = −0.123 and −0.118. Both are highly significant only because n ≈ 30,000–94,000; the
effect size is negligible. **This is the same result as the expression test (+0.093 /
+0.142), if anything slightly worse, because the sign is now marginally negative rather
than marginally positive.**

The scatter plots make it visually unambiguous: (b) and (d) show clear diagonal structure;
(a) and (c) are structureless clouds.

---

## 6. Is the negative sign an artifact of the WI-38 response construction?

Both WI-38 comparisons land at ≈ −0.12 against two *different* GM21 assays, which points at
the WI-38 response vector rather than at either GM21 side. The construction was therefore
tested directly — four alternative definitions, primary threshold, same anchors:

| WI-38 response definition | ρ vs GM21 ATAC | ρ vs GM21 K27ac | n |
|---|---|---|---|
| PDL50 − mean(25,33,37) *(used)* | **−0.123** | **−0.118** | 93,989 |
| PDL50_v_hTERT7 alone (no difference-of-differences) | −0.096 | −0.109 | 97,863 |
| PDL50 − PDL37 only | −0.155 | −0.121 | 94,136 |
| PDL45+PDL50 − mean(25,33,37) | −0.150 | −0.170 | 81,806 |

**The negative sign is robust** (−0.096 to −0.170) and does not depend on the
difference-of-differences, on which hTERT controls are used, or on whether PDL45 is
included. It is not a construction artifact.

I am **not** claiming replicative senescence actively closes what OIS opens. ρ = −0.12 is
near zero; the honest statement is *no shared structure, with a small negative tilt that
survives every reformulation I tried*.

---

## 7. ⚠️ The limitation that bounds this conclusion

**There is no within-replicative ceiling.** GSE210285 was the only second replicative
dataset and it has no peaks, so comparison (b) establishes a ceiling for **OIS only**. A
WI-38-specific measurement problem therefore cannot be fully excluded by internal evidence.

What argues against it: the WI-38 response correlates at ≈ −0.12 with *both* GM21 ATAC and
GM21 K27ac — two different assays, two different pipelines — and survives four
reformulations; and GSE175533's deposited log2FC is a replicate-aware differential table
from the original authors, not something reconstructed here.

What would settle it: a second replicative-senescence ATAC or H3K27ac dataset with
deposited peaks, to build a within-replicative ceiling. **Until that exists, the strict
claim is "replicative and OIS enhancer responses do not correlate," not "replicative
enhancer responses are measurable and opposite."**

---

## 8. Accessible ≠ active — H3K27ac should be the primary target

Within GM21, where both assays exist on 54,498 shared anchors:

| Threshold | Differential in ATAC | Also differential in H3K27ac | Same direction |
|---|---|---|---|
| \|log2FC\| ≥ 0.5 | 20,515 | **13,330 (65.0%)** | 12,432 (60.6% of all; **93.3% of the co-differential**) |
| \|log2FC\| ≥ 1.0 | 5,339 | 3,376 (63.2%) | 3,376 (63.2%; **~100% of co-differential**) |

**~35% of differential-accessibility regions show no matching H3K27ac change.** Those are
the poised/open-but-silent regions the brief anticipated. When both assays *do* move, they
agree in direction ~93% of the time.

**Conclusion: accessibility alone would mislabel roughly one third of candidate regions as
active when they are merely open. H3K27ac should be the primary target, with accessibility
as a supporting filter.** Note this also means the assay choice is not free — but §5 shows
it costs far less than the mechanism choice (ρ 0.607 vs 0.653 for swapping assay, versus
0.653 → −0.12 for swapping mechanism).

---

## 9. Sanity check — CDKN2A and CDKN1A, promoter vs distal

This is a direct test of the enhancer premise: prior work established both **promoters** are
flat for accessibility, so the pivot rests on **distal** regions moving.

**CDKN2A** (chr9:21,967,751–21,995,300 hg19 ± 200 kb), 77 anchors:

| Region class | n | WI-38 ATAC | GM21 ATAC | GM21 K27ac | IMR90 K27ac |
|---|---|---|---|---|---|
| promoter-proximal (≤2 kb) | 10 | −0.200 | −0.020 | **−1.029** | **−1.068** |
| **distal** | 67 | −0.100 | **+0.125** | **+0.240** | **+0.490** |

**CDKN1A** (chr6:36,644,237–36,655,116 hg19 ± 200 kb), 68 anchors:

| Region class | n | WI-38 ATAC | GM21 ATAC | GM21 K27ac | IMR90 K27ac |
|---|---|---|---|---|---|
| promoter-proximal (≤2 kb) | 9 | +0.576 | +0.286 | +0.080 | −0.134 |
| **distal** | 59 | −0.082 | **+0.399** | **+0.173** | **+0.186** |

Two things, and they point in opposite directions for the pivot:

1. **The enhancer premise is confirmed — for OIS.** At CDKN2A, promoter H3K27ac goes *down*
   (−1.03, −1.07) while distal H3K27ac goes *up* (+0.24, +0.49) in both independent OIS
   studies. Distal regulation is where the OIS signal lives, exactly as the pivot assumes.
2. **The premise is not confirmed for replicative.** WI-38 distal is flat-to-slightly-
   negative at both loci (−0.100, −0.082) while every OIS measure is positive. The same
   mechanism split seen genome-wide reappears at the two canonical loci.

---

## 10. Answers to the four questions asked

**Do enhancer responses diverge across mechanisms the way gene responses did?**
**Yes — equally badly, arguably marginally worse.** ρ = −0.123 (ATAC vs ATAC) and −0.118
(ATAC vs H3K27ac), against +0.093/+0.142 for gene expression. Chromatin does not rescue the
divergence; it reproduces it.

**How does cross-mechanism compare to the within-mechanism ceiling?**
Ceiling **+0.653** across different labs and cell lines. Cross-mechanism **−0.12**. The gap
is **~0.78 of correlation, with a sign flip.** Because the ceiling is high and measured on
the harder configuration (different labs, different cell lines, one side with no deposited
peaks), the low cross-mechanism value is a statement about biology, not about noise.

**Is accessibility alone adequate, or is H3K27ac necessary?**
**H3K27ac should be primary.** 35% of differential-accessibility regions show no H3K27ac
change — open but not activated. Where both move they agree 93% of the time, so
accessibility is a reasonable *filter* but a lossy *target*.

**Should the enhancer pivot target a single mechanism, or is a unified target defensible?**
**Single mechanism, decided from the outset. A unified enhancer target is not defensible.**
Pooling replicative and OIS enhancer responses would average two essentially uncorrelated
(slightly anti-correlated) programmes — the same failure that cost the expression model
0.095 AUC and its shuffled-label control when the contrast was unified.

---

## 11. What this means for the project

The enhancer pivot is **not** killed — but its scope is now fixed rather than open:

- **Within OIS, the signal is strong and reproducible across labs and cell lines
  (ρ = +0.65).** That is a materially better foundation than anything the expression route
  produced, and it is measured on two genuinely independent studies.
- **The mechanism-divergence claim now extends from expression to chromatin**, on
  ~94,000 regions with a high positive control. That strengthens the planned preprint's
  central claim considerably: it is not an artifact of gene-level quantification or of
  microarray/RNA-seq harmonization, because it reproduces in an entirely different data
  modality with a different failure profile.
- **Design against OIS specifically**, using H3K27ac as the target with accessibility as a
  filter. GM21 and IMR90 give two independent OIS studies that agree at ρ = 0.65 — a real
  cross-study training and validation split, which the expression route never had.
- **Do not spend further effort on a unified senescence target** in any modality. Three
  independent tests now say the same thing: expression ρ ≈ 0.1, model transfer at chance,
  enhancer ρ ≈ −0.12.

Before building, the one gap worth closing is a **second replicative dataset with deposited
peaks**, to establish a within-replicative ceiling. Without it, the claim is bounded as in §7.

---

## 12. Outputs

`output/enhancer_test/`

| File | Contents |
|---|---|
| `region_responses.csv` | 237,824 anchors × 15 cols — coordinates, liftover status, WI-38 annotation, all four responses |
| `anchors.csv` | anchor definitions, hg19 + lifted hg38 midpoints, assay membership |
| `signal_matrix.csv` | raw per-anchor mean coverage, 16 bigWig samples |
| `cross_mechanism_correlations.csv` | all comparisons × all thresholds |
| `batch_pca_stats.csv` | PCA variance and R² decomposition |
| `locus_sanity_checks.csv` | CDKN2A / CDKN1A promoter vs distal |
| `harmonization_stats.json`, `step2_stats.json` | counts, liftover rates, active-vs-open |
| `STEP1_LOG.txt`, `STEP2_LOG.txt` | full run logs |

`output/figures/enhancer_cross_mechanism_scatter.png` — the four scatter panels
`output/figures/enhancer_batch_pca.png` — batch structure before/after centering

Scripts: `scripts/50_fetch_enhancer_signal.sh`, `scripts/51_build_enhancer_regions.py`,
`scripts/52_enhancer_cross_mechanism.py`

### Assumptions recorded

1. **Peak-caller differences** between studies are handled by fixed-width recentering, not
   by trusting deposited peak boundaries. Peak *scores* are never used — only bigWig
   coverage — except for GSE175533, which supplies its own differential statistics.
2. **GSE175533's response is on a different footing** from the other three (deposited
   replicate-aware log2FC vs coverage ratio). Tested four ways in §6.
3. **GM21 senescent = RAS D18+D23**, excluding D32+ as post-escape, carried over from the
   expression build's measured finding.
4. **IMR90 arms** are the 2 proliferating and 2 senescent H3K27ac bigWigs; the quiescent
   arm exists in that series but is not used here.
5. **Pseudocount 0.10** after mean-1 normalization; results are stable across the four
   thresholds tested, and the pseudocount only affects near-zero anchors.
6. **≥2-sample consensus support** — a union of single-sample calls would be dominated by
   one-off peaks.
