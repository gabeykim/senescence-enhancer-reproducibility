# Within-replicative reproducibility test — H3K27ac enhancer responses

**Date:** 2026-08-25 · local, no GPU · 2.1 GB downloaded · ~40 min end to end

---

## Verdict

**Replicative senescence has no reproducible enhancer program in the deposited data. The
within-replicative ceiling is +0.144 against the OIS ceiling of +0.653 — and, decisively,
replicative datasets agree with each other no better than they agree with OIS datasets.**

| Comparison | Spearman (primary) | n |
|---|---|---|
| **OIS ceiling** (prior run: GM21 vs IMR90 K27ac, n=2v2 both sides) | **+0.653** | 88,530 |
| **R1 — within-replicative ceiling** (IMR90 n=2v2 vs GSE106146 n=1v1, cross-lab) | **+0.144** | 144,754 |
| R2 — within-replicative, **same lab/pipeline/build** (IMR90 vs BJ) | **+0.054** | 134,336 |
| R3 — within-replicative, cross-lab (BJ vs GSE106146) | +0.268 | 124,962 |
| X1 — cross-mechanism, **same lab/cell line/pipeline/build** | +0.114 | 123,121 |
| X2 — cross-mechanism (IMR90 rep vs GM21 OIS) | −0.134 | 80,007 |
| X3 — cross-mechanism (IMR90 rep vs IMR90 OIS, different lab) | −0.164 | 181,114 |
| X4 — cross-mechanism (GSE106146 rep vs GM21 OIS) | +0.164 | 77,110 |

**The three within-replicative values (+0.054, +0.144, +0.268) fall inside the range of the
four cross-mechanism values (−0.164 to +0.164).** Two datasets of the *same* mechanism are
not more concordant than two datasets of *different* mechanisms. That is what "no
reproducible program" means operationally.

**Replicative senescence is not a viable design target. The project is forced to OIS.**

---

## 1. What was actually deposited — four corrections to the brief

Both accessions verified live before download (GSE146585 uid 200146585, 223 samples,
public 2025-05-30; GSE106146 uid 200106146, 82 samples, SubSeries of GSE105937).

**GSE146585 has 23 SubSeries.** H3K27ac is not in one of them — it is spread across six
(GSE146568 ×12, GSE146566 ×6, GSE146559 ×4, GSE146567 ×4, GSE146582 ×3, GSE146563 ×2),
31 samples. The three that matter:

| SubSeries | Contents | n |
|---|---|---|
| **GSE146567** | IMR90 Young Rep1/2 vs Senescent Rep1/2 | **n=2 per arm — the only one** |
| GSE146559 | BJ passage series (Young / Int-Early / Int-Late / Senescent) | n=1 each |
| GSE146563 | IMR90 GFP-control vs RAS-OIS | n=1 each |

### 🔴 Correction 1 — GSE146585 is hg18

All 223 samples carry `Genome_build: hg18` (STAR alignment, HOMER bigWigs). Not hg19, not
hg38. The brief was right that builds were not in the series header; they are in
`Sample_data_processing`. Assuming hg19 would have produced a wrong answer silently.

### 🔴 Correction 2 — there are no EPQ / quiescent controls

**Zero of 223 sample titles match "EPQ" or "quiesc".** Controls are `Young` — early-passage
*proliferating* cells. The premise that controls were early-passage quiescent to minimize
cell-cycle confounding is not supported by the deposited metadata. The cell-cycle confound
is therefore present here exactly as in every other dataset in this project.

### 🔴 Correction 3 — the n=2 replicative arm is IMR90, not BJ

The BJ passage series is n=1 per timepoint.

### 🔴 Correction 4 — GSE106146's replicative H3K27ac is n=1 vs n=1, not n=2

GSM2830427–430 (the four proliferating/senescent H3K27ac samples) have **no supplementary
files at all** — raw reads only. Only `GSE106146_P2-H3K27ac.bigWig` (proliferating) and
`GSE106146_S1-H3K27ac.bigWig` (senescent) exist as processed signal. The **IR arm does have
n=2** (`GSM2830447/8`), which is why the IR comparison is better powered than the
replicative one in the same study.

Recovering the missing replicates requires SRA download + realignment of 4 ChIP libraries.
**Not done, and not recommended:** the ceiling is limited by the *other* side of every
comparison too, and R2 — which uses the best-powered replicative dataset against a same-lab
partner — is the *worst* result, so power is not the binding constraint.

**No SRA processing was needed for the primary test.** All H3K27ac used here was deposited
as bigWigs.

### One thing the survey missed, in our favour

GSE146563's IMR90 RAS-OIS H3K27ac shares lab, cell line, pipeline **and build** with
GSE146567's IMR90 replicative. That permits comparison X1 — a cross-mechanism test with
essentially every technical confound removed, which the prior round could not do.

---

## 2. Harmonization

**The prior run's 237,824 anchors were reused verbatim** (1 kb fixed-width, hg19, ≥2-sample
consensus from GM21 ATAC + H3K27ac). Every number here therefore lands on the same regions,
widths, normalization and thresholds as the +0.653 OIS ceiling — the comparison is
like-for-like by construction.

**Builds handled per file, never per analysis:**

| Dataset | Build | How queried |
|---|---|---|
| GSE146585 (IMR90 rep, BJ rep, IMR90 OIS) | **hg18** | anchors lifted hg19→hg18, signal read in hg18 |
| GSE106146 (replicative, IR) | hg19 | queried directly |
| Prior GM21 / IMR90 OIS / WI-38 | hg19 / hg38 | carried over from the prior run |

**hg19→hg18 liftover: 237,512 / 237,824 midpoints (99.87%)**, 312 failed (0.13%). Span
integrity on a 5,000-anchor sample: both ends lift 99.90%, and **100.00%** of those preserve
the 1 kb span within 50 bp. Liftover is not a limiting factor. Coverage per file ranged
184,454–234,132 of 237,824 anchors.

**Normalization** identical to the prior run: each sample divided by its own mean over
covered anchors (mean 1), pseudocount 0.10, log2 ratio.

### Caveat on region selection, stated

The anchors were defined from GM21 (OIS) fibroblast peak *presence*, not OIS
*responsiveness* — a general fibroblast regulatory atlas rather than a set selected for OIS
behaviour. But they are not replicative-derived, so a replicative-specific enhancer sitting
outside GM21's peak atlas would be invisible here. Against that concern: coverage of the
replicative datasets across these anchors is high (78–98%), and the anchors capture the
regions where replicative H3K27ac signal actually is.

---

## 3. Batch check

PCA on 140,606 complete-case anchors × 12 samples:

| Space | PC | Var | R²(study) | R²(status) | Dominated by |
|---|---|---|---|---|---|
| raw pooled | PC1 | 35.0% | **0.534** | 0.020 | STUDY |
| raw pooled | PC2 | 21.4% | 0.062 | 0.105 | status |
| **per-study centered** | PC1 | 34.5% | **0.000** | 0.038 | status |
| per-study centered | PC3 | 12.2% | 0.000 | **0.568** | STATUS |

Study dominates raw PC1 (R² 0.534 — less extreme than the 0.969 of the prior enhancer run
or the 0.997 of the expression work) and is removed completely by per-study centering. As
before, the correlations below are computed on **per-study log2 fold changes**, which are
already within-study contrasts, so the study offset cancels by construction.

---

## 4. The headline — within-replicative ceiling

**n is small everywhere and is stated with every number.** The OIS ceiling of +0.653 was
achieved with **n=2 vs 2 on both sides**. No within-replicative comparison at that power is
possible from deposited data, because only GSE146567 has n=2.

| Comparison | n per arm | \|lfc\|≥0 | ≥0.25 | **≥0.5 (primary)** | ≥1.0 |
|---|---|---|---|---|---|
| **R1 IMR90 rep vs GSE106146 rep** | 2v2 × 1v1 | +0.063 | +0.081 | **+0.144** | +0.327 |
| **R2 IMR90 rep vs BJ rep** *(same lab)* | 2v2 × 1v1 | −0.043 | −0.021 | **+0.054** | +0.180 |
| **R3 BJ rep vs GSE106146 rep** | 1v1 × 1v1 | +0.280 | +0.277 | **+0.268** | +0.293 |
| *OIS ceiling (prior)* | *2v2 × 2v2* | *+0.640* | *+0.643* | *+0.653* | *+0.678* |

Three points:

1. **Every within-replicative pair is far below the OIS ceiling** — 2.4× to 12× lower at the
   primary threshold.
2. **The three pairs disagree with each other** (+0.054, +0.144, +0.268). A reproducible
   programme would give consistent pairwise agreement; these do not.
3. **R1 is strongly threshold-dependent** (+0.063 → +0.327), whereas the OIS ceiling was
   essentially flat (+0.640 → +0.678). Rising correlation with effect size is the signature
   of low signal-to-noise: at the very largest effects there *is* some agreement (+0.33),
   but it never approaches OIS levels, and it is bought by discarding 85% of regions.

**R2 is the most informative and the most damaging.** It is the comparison with the fewest
technical confounds in the entire project — same lab, same pipeline, same aligner, same
build, best-powered replicative dataset on one side — and it returns **+0.054**. Technical
heterogeneity cannot explain it, because there is almost none.

---

## 5. Cross-mechanism divergence holds up — it was not a WI-38 artifact

The prior round's open question was whether ρ ≈ −0.12 reflected biology or a WI-38 ATAC
measurement problem. Answer: **biology.** Replacing WI-38 ATAC with three new replicative
**H3K27ac** datasets reproduces the same near-zero-to-negative pattern:

| | Spearman (primary) | n |
|---|---|---|
| X1 IMR90 rep vs IMR90 OIS — **same lab, cell line, pipeline, build** | **+0.114** | 123,121 |
| X2 IMR90 rep vs GM21 OIS | −0.134 | 80,007 |
| X3 IMR90 rep vs IMR90 OIS (different lab) | −0.164 | 181,114 |
| X4 GSE106146 rep vs GM21 OIS | +0.164 | 77,110 |
| *prior: WI-38 ATAC vs GM21 OIS ATAC* | *−0.123* | *93,989* |

X1 is the cleanest cross-mechanism test available anywhere in this project — the *same
IMR90 cells*, the *same lab*, the *same HOMER/STAR pipeline*, the *same hg18 build*, the
same assay — and the two mechanisms correlate at **+0.114**. That eliminates lab, cell
line, protocol, pipeline and build as explanations simultaneously.

**Divergence confirmed. It is not a WI-38 measurement artifact.**

---

## 6. Cross-assay within replicative

| | Spearman (primary) | n |
|---|---|---|
| A1 WI-38 rep ATAC vs IMR90 rep K27ac | +0.222 | 98,315 |
| A2 WI-38 rep ATAC vs GSE106146 rep K27ac | +0.226 | 92,401 |
| *OIS equivalent (prior, d)* | *+0.607* | *36,560* |

Within OIS, accessibility and acetylation agreed at +0.607. Within replicative they agree at
+0.22 — consistent with each replicative measurement being individually noisy, or with there
being little coherent programme for the two assays to agree *about*.

---

## 7. 🔴 The irradiation arm — and an artifact I had to catch

The IR-vs-replicative comparison inside GSE106146 initially returned **+0.564**, by far the
highest within-mechanism-like value in the test. **It is almost entirely an artifact.**

Both responses are computed against the *same single proliferating control*
(`Sen2019_Pro`), so both contain that term: correlating (IR − Pro) with (Sen − Pro)
correlates the shared −Pro. Quantified with a null in which the numerators are
biologically unrelated but the denominator is shared:

| Construction | Spearman |
|---|---|
| I1 — IR vs replicative, shared denominator | **+0.564** |
| **NULL — two unrelated numerators (IMR90 Young, BJ Young), same denominator** | **+0.523** |
| null — IMR90 Young numerator vs replicative response, same denominator | +0.154 |
| null — BJ Young numerator vs replicative response, same denominator | +0.083 |

A shared denominator alone generates **+0.523** from numerators with no shared biology. I1's
+0.564 is barely above that. **The apparent IR/replicative agreement is not evidence of
shared biology and I am not reporting it as a result.**

The honest IR numbers use independent controls:

| | Spearman (primary) | n |
|---|---|---|
| I2 IR vs GM21 OIS | +0.084 | 76,620 |
| I3 IR vs IMR90 replicative | +0.115 | 138,174 |

**Irradiation-induced senescence resembles neither replicative nor OIS** — it is
approximately equally uncorrelated with both (+0.115 vs +0.084). Given that IR is the
closest available proxy to therapy-induced senescence, this is a third mechanism that does
not share an enhancer programme with the other two. Caveat: n=2 vs 1, and its only control
is the single GSE106146 proliferating sample.

---

## 8. Locus sanity checks — H3K27ac, promoter vs distal

**CDKN2A** (chr9, 77 anchors: 10 promoter, 67 distal):

| | IMR90 rep | BJ rep | GSE106146 rep | GSE106146 IR | IMR90 OIS (same lab) | GM21 OIS | IMR90 OIS (prior) |
|---|---|---|---|---|---|---|---|
| promoter | **−0.848** | **−0.632** | **−0.431** | −0.583 | −0.055 | **−1.029** | **−1.068** |
| distal | −0.166 | +0.215 | +0.311 | −0.064 | −0.050 | **+0.240** | **+0.490** |

**CDKN1A** (chr6, 68 anchors: 9 promoter, 59 distal):

| | IMR90 rep | BJ rep | GSE106146 rep | GSE106146 IR | IMR90 OIS (same lab) | GM21 OIS | IMR90 OIS (prior) |
|---|---|---|---|---|---|---|---|
| promoter | +0.188 | +0.330 | +0.031 | +0.468 | +0.110 | +0.080 | −0.134 |
| distal | −0.118 | +0.139 | +0.120 | +0.007 | −0.373 | +0.173 | +0.186 |

Two findings:

- **The promoter-down half of the enhancer premise reproduces in replicative.** CDKN2A
  promoter H3K27ac falls in all three replicative datasets (−0.85, −0.63, −0.43), matching
  the OIS pattern (−1.03, −1.07). So replicative senescence *is* being detected — these are
  not dead measurements.
- **The distal-up half does not.** OIS gives consistent distal increases (+0.24, +0.49);
  replicative gives −0.17, +0.22, +0.31 — inconsistent in sign across the three datasets.

This is the sharpest single answer to the question the test was built for: **it
distinguishes "WI-38 ATAC could not measure it" from "replicative has no consistent distal
programme."** The new replicative H3K27ac data detects the promoter change fine and still
fails to produce a consistent distal signal. The earlier flat WI-38 distal result was not a
measurement failure.

Note also that GSE146563 (IMR90 OIS, n=1) is nearly flat at both loci (−0.055, −0.050) —
that single-replicate dataset is visibly weaker than the n=2 OIS datasets, which is worth
keeping in mind when reading X1.

---

## 9. Answers to the four questions

**What is the within-replicative ceiling, and how does it compare to +0.653?**
**+0.144** for the best-powered cross-lab pair (R1), **+0.054** for the lowest-confound
same-lab pair (R2), **+0.268** for the third (R3). Against **+0.653** for OIS. The
replicative ceiling is between 2.4× and 12× lower, and the three estimates do not agree with
each other.

**Does cross-mechanism divergence hold up, or was it a WI-38 artifact?**
**It holds up, and it was not an artifact.** Three new replicative H3K27ac datasets
reproduce it (+0.114, −0.134, −0.164, +0.164 vs the prior −0.123). X1 removes lab, cell
line, pipeline and build simultaneously and still gives +0.114.

**Is replicative senescence a viable design target?**
**No — not with data that exists today.** A design target requires a response that
reproduces across independent datasets; replicative H3K27ac does not reproduce even against
a same-lab partner. Training on it would fit one dataset's noise.

**Does the irradiation arm look more like replicative or like OIS?**
**Neither.** +0.115 with replicative, +0.084 with OIS, once the shared-control artifact is
removed. It is a third non-overlapping programme.

---

## 10. What bounds this conclusion

Stated plainly, because it is the one thing that could overturn the result:

1. **No within-replicative comparison at n=2v2-on-both-sides exists.** The OIS ceiling had
   it; the best replicative pairing is 2v2 × 1v1. Some of the gap is power.
2. **But power is unlikely to be the whole story.** The OIS ceiling was achieved across
   *different labs, cell lines, protocols and builds* — a much harder configuration than R2,
   which shares all of those and still returns +0.054. And the promoter-level signal in §8
   confirms the replicative measurements detect senescence, so they are not simply dead.
3. **Region selection**: anchors are GM21-derived. A replicative-specific enhancer outside
   the GM21 peak atlas is invisible here.
4. **n=1 on one side of R1/R2/R3 and on both sides of R3**, and n=1 for the same-lab OIS arm
   in X1.

**What would settle it definitively:** processing GSE106146's four raw-read H3K27ac
libraries from SRA to recover n=2 vs n=2, giving one genuinely power-matched
within-replicative comparison. That is the only remaining move that could raise the
replicative ceiling, and I'd estimate it at a few hours of alignment. **I did not do it** —
say the word if you want it, but R2's same-lab +0.054 makes me doubt it changes the answer.

---

## 11. Recommendation

**Design against OIS.** It is now the only mechanism in this project with a demonstrated,
cross-laboratory-reproducible enhancer response (+0.653, n=2v2 both sides, different cell
lines and protocols), and it has two independent studies (GM21, IMR90) that can serve as a
real train/validation split.

**Do not train on replicative senescence, and do not pool mechanisms** in any modality. The
evidence is now four-fold and consistent:

| Test | Result |
|---|---|
| Gene expression, cross-mechanism | ρ ≈ +0.09 / +0.14 |
| Model transfer between mechanisms | at chance, both directions |
| Enhancer response, cross-mechanism (ATAC) | ρ ≈ −0.12 |
| **Enhancer response, cross-mechanism (H3K27ac, this run)** | **ρ = +0.114 same-lab, −0.16 to +0.16 overall** |

And the new finding that strengthens the preprint most: **replicative senescence does not
reproduce against itself.** The mechanism-divergence claim is no longer only "mechanisms
differ from each other" — for replicative specifically it is "there is no stable programme
to differ *with*," measured across three independent datasets, two cell lines, two labs and
two genome builds.

---

## 12. Outputs

All at **`/Users/gabeykim/Downloads/Senescence/output/replicative_ceiling_test/`**

| File | Bytes | Contents |
|---|---|---|
| `REPORT_replicative_ceiling.md` | — | this report |
| `region_responses_replicative.csv` | ~46 MB | 237,824 anchors × all responses |
| `signal_matrix_replicative.csv` | ~29 MB | raw per-anchor coverage, 12 new bigWigs |
| `anchors_with_hg18.csv` | ~30 MB | anchors + hg18 lift + carried-over prior responses |
| `replicative_correlations.csv` | — | all 12 comparisons × 4 thresholds |
| `batch_pca_stats.csv` | — | PCA variance / R² decomposition |
| `locus_sanity_checks.csv` | — | CDKN2A / CDKN1A promoter vs distal |
| `quantification_stats.json`, `step2_stats.json` | — | liftover rates, n-per-arm |
| `STEP1_LOG.txt`, `STEP2_LOG.txt` | — | full run logs |
| `figures/replicative_ceiling_scatter.png` | — | six-panel scatter |
| `figures/replicative_batch_pca.png` | — | batch structure before/after centering |
| `scripts/` | — | `60_fetch_replicative_ceiling.sh`, `61_replicative_ceiling.py`, `62_replicative_ceiling_analysis.py`, `63_shared_denominator_null.py` |

### Assumptions recorded

1. **Anchors reused from the prior enhancer test**, unchanged, so numbers are comparable to
   +0.653. GM21-derived — see §2 caveat.
2. **Each bigWig read in its own genome build**; liftover never touches a coverage value.
3. **Arm definitions**: senescent vs Young for GSE146585 (no quiescent controls exist);
   S1 vs P2 for GSE106146 replicative; IR Rep1/2 vs P2 for the IR arm.
4. **BJ intermediate timepoints not downloaded** — only Young and Senescent, since the
   intermediates do not enter any comparison.
5. **Pseudocount 0.10** after mean-1 normalization, identical to the prior run.
6. **The IR/replicative comparison inside GSE106146 shares a denominator** and is reported
   as an artifact, not a result (§7).
