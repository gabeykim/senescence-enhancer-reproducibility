# Senescence expression training set from GSE175533

**Date:** 2026-08-22
**Scope:** Local data curation only. No GPU instance, no model loading, no SRA access.
**Build script:** [`scripts/30_build_senescence_trainset.py`](../../scripts/30_build_senescence_trainset.py)
(single script, fully reproducible, all judgment calls as module-level constants).

## Verdict

**The training set is ready to fine-tune on.** All six sanity checks pass, all
structural contract checks pass, and the biology is correct in the direction the
literature demands. Two blocking-adjacent caveats are recorded below (the
quiescent class is not a textbook quiescence contrast; the ON arm has n=9
samples), neither of which prevents training but both of which constrain what
can be claimed from the result.

Primary artifact: **`senescence_trainset.h5ad`** — 33 tasks × 19,652 genes,
float32, hg38, chromosome-split with chr9+chr6 fully held out.

---

## Step 1 — Acquire

Verified against the live GEO listing rather than assumed. All three expected
files are present, and all three cached copies **byte-match** the remote:

| File | Bytes | Status |
|---|---|---|
| `GSE175533_atlas.bed.gz` | 4,191,769 | matches remote |
| `GSE175533_atac_peaks_sig.xlsx` | 68,111,647 | matches remote |
| `GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx` | 46,113,681 | matches remote |

(GEO also hosts `GSE175533_RAW.tar` and `GSE175533_sceasy_hay.h5ad.gz`; neither
is needed here and neither was downloaded. No SRA access.)

### Genome build — confirmed hg38 on every file used, three independent ways

| Evidence source | Value |
|---|---|
| SOFT `Sample_data_processing` | `Genome_build: Hg38` |
| SOFT RNA alignment text | "Salmon generated index based on Hg38" |
| `atac_peaks_sig.xlsx` readme sheet | "region \| hg38 coordinates" |
| Annotation used (GENCODE v44) | GRCh38 |

**No hg19 file is used anywhere. No liftover is required** and
`hg19ToHg38.over.chain.gz` is not needed for this build. No builds are mixed.

### Data defect found and handled: Excel date corruption

**27 gene symbols in the GEO-deposited TPM workbook have been destroyed by
Excel's date auto-conversion** — the classic MARCH1–11 → `2021-03-01..11`,
SEPT1–14 → `2021-09-01..14`, DEC1 → `2021-12-01` corruption. They arrive as
`datetime.date` objects, not strings.

These are **dropped and counted, not guessed at**: `2021-03-01` could be MARCH1
or MARC1, and the canonical MARCH*/SEPT* symbols are otherwise absent from the
table, so there is no in-file evidence to disambiguate. **No marker gene is
affected.** This is a defect in the deposited file, not in this pipeline, and
anyone else using this workbook will hit it.

---

## Step 2 — Arms (every judgment call stated, with reasoning)

All 96 RNA samples across four arms were assigned. Full per-sample table:
[`sample_assignments.csv`](sample_assignments.csv).

| Arm | Class | Role | n |
|---|---|---|---|
| RS | **proliferating (OFF)** | train_eligible | **15** |
| RS | **senescent (ON)** | train_eligible | **9** |
| RS | transition_ambiguous | excluded | 6 |
| CD | **quiescent** | train_eligible | **9** |
| CD | cd_not_yet_arrested | excluded | 21 |
| RIS | ris_senescent | holdout_ris | 15 |
| RIS | ris_baseline | holdout_ris | 3 |
| hTERT | htert_immortalized | excluded_control | 18 |

**ON (n=9):** PDL50 ×3, PDL52 ×3, PDL53 ×3
**OFF (n=15):** PDL20, PDL25, PDL28, PDL33, PDL37 (×3 each)
**QUIESCENT (n=9):** CD day 4, 7, 10 (×3 each)

### PDL boundary — `OFF_MAX_PDL=37`, `ON_MIN_PDL=50`

The prior audit located the inflection **between PDL46 and PDL50**. OFF is set
at ≤37 rather than ≤46 to leave deliberate margin below a boundary whose exact
position is uncertain to within a few PDL, so the OFF arm is unambiguously
pre-transition.

### Transition zone (PDL45, PDL46) — **EXCLUDED** (`TRANSITION_POLICY="exclude"`)

Ambiguous by construction and the single most likely source of label noise. A
binary head trained with genuinely intermediate samples forced to one side
learns a blurred boundary — and the blur lands exactly where the biology of
interest lives. Costs 6 samples; OFF still retains 15 and ON 9.
Configurable: `assign_off` / `assign_on` available.

### hTERT — **EXCLUDED from training, retained as evaluation control** (`HTERT_POLICY="exclude"`)

Tempting as 18 extra OFF samples, but they are *immortalized*, not merely
young. Constitutive telomerase is itself a major transcriptional perturbation,
and these cells were carried to PDL 46–109, far outside the WI-38 range. A head
trained with hTERT as OFF risks learning "telomerase-positive immortalized
cell" as its OFF signature rather than "young proliferating primary
fibroblast" — a direct confound with the intended axis. They remain valuable as
a **negative control at evaluation** (a senescence head should score them low)
and are emitted in the sample table for that use.

### Cell density — **THIRD CLASS, late timepoints only** (`CD_POLICY="third_class"`, `CD_QUIESCENT_DAYS=[4,7,10]`)

**This decision was corrected after measuring the arm rather than assuming it.**
The CD arm is a **progressive arrest time course**, not a uniform quiescent
state:

| CD day | D1 | D1.5 | D2 | D2.5 | D2.75 | D3 | D3.5 | **D4** | **D7** | **D10** |
|---|---|---|---|---|---|---|---|---|---|---|
| mean MKI67 TPM | 24.2 | **53.1** | 46.3 | 32.1 | 23.2 | 17.6 | 11.5 | **6.0** | **0.9** | **4.5** |

Early CD timepoints are *actively proliferating* — D1.5 MKI67 (53.1) is more
than 3× the proliferating PDL20–37 mean (15.5). Calling them "quiescent" would
have been simply wrong. Restricting to D4/D7/D10 gives a defensible
arrested-but-not-senescent class.

**Caveat recorded, not smoothed over:** late-CD cells show **higher** IL6
(42.3 vs 15.8) and CXCL8 (29.8 vs 5.5) than the senescent arm — contact-inhibited
dense cultures mount their own inflammatory response. So this third class
teaches *"senescence vs dense arrested culture"*, **not** the textbook
senescence-vs-quiescence contrast. A head trained on it should not be claimed to
separate senescence from quiescence in general. Set `CD_POLICY="exclude"` if
that ambiguity is unacceptable.

### RIS — **HELD OUT as a transfer test** (`RIS_POLICY="holdout"`)

RIS is the closest available proxy to the project's eventual
doxorubicin/therapy-induced validation model. Its value as a held-out
generalization probe far exceeds its value as ~15 extra training samples.
Training on it would destroy the one chance this dataset offers to ask *"does a
head trained on replicative senescence transfer to DNA-damage-induced
senescence?"* — the exact question the wet-lab plan depends on. Emitted with
`role=holdout_ris`. (Note `d0_A/B/C` are the pre-irradiation mock baseline, not
senescent.)

---

## Step 3 — Intervals

TSS-centered 524,288 bp windows, per the measured contract. TSS taken from
**GENCODE v44 (GRCh38)**: gene start for `+` strand, gene end for `−` strand.
Where a symbol maps to multiple gene records, the longest-span record on a
standard chromosome is kept.

| Stage | Count |
|---|---|
| Genes in TPM table (after date-corruption drop) | 33,487 |
| Genes with TPM>0 in ≥1 ON/OFF sample | 27,499 |
| **Intervals built** | **19,652** |
| Dropped: no GENCODE annotation | 7,722 |
| Dropped: window past chromosome end | 77 |
| Dropped: window before chromosome start | 48 |

### Chromosome-boundary handling — **125 genes dropped, explicitly**

Genes whose 524,288 bp window would run past a chromosome end (77) or before
position 0 (48) are **dropped and counted, never clipped**. Clipping would shift
the TSS off-center and silently break the contract — a window that isn't
TSS-centered is worse than a missing window. Full list:
[`dropped_genes.csv`](dropped_genes.csv).

### The 7,722 "no annotation" drops — characterized, not just accepted

**90.3% (6,974) are clone-based identifiers** (`AC000032.1`, `AL######.#`,
`AB015752.1` …) — unnamed loci from the older annotation the submitters
quantified against. They have no named TSS to center on and are legitimately
unusable. The remaining 748 are named symbols absent from GENCODE v44, mostly
**HGNC renames** (AARS→AARS1, ADSS→ADSS2, CARS→CARS1, EPRS→EPRS1…). Recovering
those would require an HGNC previous-symbol mapping. At **<3% of expressed
genes with zero marker genes affected**, this is a quantified, acceptable loss.

### Accessibility as a region prior (optional, off by default)

Accessibility is used **only to decide which regions are worth including —
never as a regression target.** 19,944 genes are nearest-neighbours of
significant late-PDL differential peaks (of 363,470 atlas peaks).

| Variant | Genes | train | val | test | CDKN2A/CDKN1A retained |
|---|---|---|---|---|---|
| **Default (no restriction)** | **19,652** | 15,910 | 2,002 | 1,740 | both |
| `--restrict-to-diff-accessibility` | 15,120 | 12,267 | 1,531 | 1,322 | both |

**Does it change the set meaningfully?** It removes 4,532 genes (23%) — a real
but not dramatic reduction that still leaves ample volume, and both benchmark
loci survive. **Default is OFF**, because the prior audit established that
CDKN2A/CDKN1A promoters show *no* significant accessibility change despite
strong CDKN1A RNA induction — i.e. accessibility demonstrably fails to flag
the very loci that matter most here. Filtering by it risks discarding
expression-informative genes for exactly the wrong reason. Both variants are
built and shipped; use the restricted one only if training cost dominates.

---

## Step 4 — Chromosome split

| Split | Chromosomes | Genes |
|---|---|---|
| **test** | **chr9, chr6** | **1,740** |
| **val** | chr8, chr16, chr20 | 2,002 |
| **train** | chr1–5, chr7, chr10–15, chr17–19, chr21, chr22, chrX, chrY | 15,910 |

Lists are **disjoint by construction** (asserted in code — gReLU's `split()`
returns an empty train set silently if they aren't). chr9 and chr6 are held out
**entirely**, so CDKN2A (chr9) and CDKN1A (chr6) performance is a genuine
generalization test rather than memorization. Verified: chr9/chr6 genes appear
in the `test` split and nowhere else.

---

## Step 5 — Format

`senescence_trainset.h5ad`, matching the measured Task 4 contract:

| Contract item | Value in artifact |
|---|---|
| `.X` shape (n_tasks, n_intervals) | **(33, 19652)** |
| dtype | **float32** |
| `.obs` | sample, arm, cls, pdl, day (33 rows) |
| `.var` | chrom, start, end, strand, tss, gene_type, **split** |
| Interval width | **exactly 524,288 bp**, all 19,652 |
| Genome | hg38 |
| `.uns` | contract, target_transform, arm_policy, arm_policy_rationale, splits, provenance |

The full reasoning for every step-2 judgment call is embedded in
`.uns["arm_policy_rationale"]`, so the artifact carries its own provenance.

### Target normalization: **log1p(TPM)** — stated explicitly

**What was applied:** `log1p(TPM)`. Measured output range **0.000 to 9.970**.

**Why, and why not Borzoi's own transform:** Borzoi's native targets live in a
variance-stabilized space (per-bin `**0.75` power + soft-clip + scale — the
chain Task 2 had to invert). That chain is **not** reproduced here, and
deliberately so: it applies to *per-bin coverage tracks*, whereas this head
predicts a *per-gene scalar*. Applying a bin-level transform to a gene-level
scalar would be cargo-culting.

What *does* carry over is the **principle**. Raw TPM spans ~6 orders of
magnitude and is dominated by a handful of very-high-expression genes; Task 2
measured that prediction error already concentrates in low-expression genes, and
an untransformed target would worsen that by letting a few housekeeping genes
dominate MSE. `log1p` is monotonic, defined at zero, and — critically — keeps
targets **non-negative**, matching the measured `ConvHead` **softplus** output
activation. A z-scored target would place mass below zero that softplus
structurally cannot emit; that is the specific failure mode where a model trains
happily and predicts nonsense.

---

## Step 6 — Sanity checks

### Marker checks (ON vs OFF) — **ALL PASS**, run as a gate before building

The build **exits non-zero and writes nothing** if any gated marker fails.

| Gene | Expected | ON mean TPM | OFF mean TPM | log2FC | p | Result |
|---|---|---|---|---|---|---|
| **CDKN1A** | UP | 1159.69 | 814.98 | **+0.51** | 1.78e-02 | **PASS** |
| **IL6** | UP | 15.79 | 0.77 | **+4.34** | 2.92e-03 | **PASS** |
| **CXCL8** | UP | 5.45 | 0.18 | **+4.86** | 3.01e-03 | **PASS** |
| **LMNB1** | DOWN | 1.54 | 24.93 | **−4.01** | 7.49e-07 | **PASS** |
| **MKI67** | DOWN | 0.88 | 15.47 | **−4.13** | 1.09e-06 | **PASS** |
| CDKN2A | report only | 69.11 | 57.92 | +0.25 | 2.12e-01 | see below |
| SERPINE1 | report only | 466.74 | 780.39 | −0.74 | 4.20e-07 | see below |

**CDKN1A is significantly higher in ON — the arm assignment is correct.**

### CDKN2A — reported as requested

**CDKN2A expression rises only +0.25 log2FC and is NOT significant (p=0.21).**
This is consistent with the prior audit's finding that its *promoter
accessibility* is flat across the whole PDL course: at this locus, neither
accessibility nor expression moves convincingly in this dataset. Two readings
are compatible with the data and cannot be separated here: (a) p16 induction in
replicative senescence is well documented to be mosaic across a cell
population, so a bulk measurement dilutes it; (b) WI-38 replicative senescence
may be more p21-driven than p16-driven at these PDLs. **Practical consequence:
do not expect a head trained on this data to predict CDKN2A induction, and do
not treat CDKN2A as a validation target for it.** CDKN1A is the usable
p-pathway benchmark here — which is precisely why holding out both chr9 and
chr6 matters.

**SERPINE1 goes DOWN (−0.74, p=4.2e-07)** despite being a canonical SASP
factor — the same non-canonical behaviour flagged in the earlier audit,
reproduced here. Reported, not explained away.

### Senescent vs quiescent — the third class does *not* separate cleanly

| Gene | ON mean | Quiescent mean | log2FC | p |
|---|---|---|---|---|
| CDKN1A | 1159.70 | 885.52 | +0.39 | 0.065 |
| IL6 | 15.79 | 32.29 | −1.03 | 0.153 |
| CXCL8 | 5.45 | 21.75 | **−1.99** | **0.040** |
| LMNB1 | 1.54 | 10.76 | −2.79 | 0.002 |
| MKI67 | 0.88 | 3.80 | −2.10 | 0.009 |

**3/5 separate at p<0.05, but the SASP markers do not separate in the expected
direction** — CXCL8 is significantly *higher* in quiescent. As noted in step 2,
dense contact-inhibited cultures mount their own inflammatory program. The
quiescent class is retained because distinguishing senescence from arrest is
genuinely valuable, but **what it actually teaches is "senescence vs dense
arrested culture," and the report says so rather than claiming more.**

### Structural checks — all pass

| Check | Result |
|---|---|
| Interval count == target count | 19,652 == 19,652 ✓ |
| Duplicate genes in `.var` | 0 ✓ |
| Genes appearing in >1 split | **0** ✓ |
| Chromosomes spanning >1 split | 0 ✓ |
| CDKN2A on chr9, split=test | ✓ held out |
| CDKN1A on chr6, split=test | ✓ held out |
| NaN fraction in `.X` | 0.000000 ✓ |
| Distinct interval widths | `[524288]` only ✓ |
| Targets non-negative (softplus-compatible) | ✓ |

Independently re-verified by a fresh read-back of the written `.h5ad`, not
just from build-time state.

---

## Ready to fine-tune? **Yes — with two constraints on what can be claimed**

**Ready:** the artifact matches the measured contract exactly, the biology is
correct and gated, splits are leak-free, and both benchmark loci are held out.

**Constraint 1 — ON arm is n=9 (3 PDL points × 3 replicates).** That is thin,
and the three PDL points are serially related (same culture lineage), so the
effective independent-sample count is closer to 3 than 9. Expect the head to be
easy to overfit; use the val split aggressively and treat any single-run result
with suspicion.

**Constraint 2 — the quiescent class is not a clean quiescence contrast** (see
above). If the project's claim depends on "distinguishes senescence from
quiescence," this data cannot support it as-is; rebuild with
`CD_POLICY="exclude"` and make the claim narrower.

**Not blocking, but plan for it:** RIS is held out deliberately as the
transfer test — the first real check of whether replicative-senescence training
generalizes to DNA-damage-induced senescence, which is the bridge to the
doxorubicin wet-lab plan. Run it before trusting any transfer claim.

## Artifacts

| File | Contents |
|---|---|
| `senescence_trainset.h5ad` | **Primary artifact** — 33 × 19,652, default policies |
| `senescence_trainset_daPrior.h5ad` | Accessibility-prior-restricted variant, 33 × 15,120 |
| `sample_assignments.csv` | All 96 RNA samples with arm/class/role |
| `marker_sanity_checks.csv` | ON-vs-OFF marker statistics |
| `senescent_vs_quiescent_markers.csv` | Third-class discrimination statistics |
| `dropped_genes.csv` | Every dropped gene with its reason |
| `REPORT_trainset.md` | Raw build log (measured numbers, generated by the script) |

**Reproduce:**
```bash
python3 scripts/30_build_senescence_trainset.py                              # default
python3 scripts/30_build_senescence_trainset.py --restrict-to-diff-accessibility \
        --out-prefix senescence_trainset_daPrior                             # variant
```
