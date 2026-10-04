# Audit: GSE254358 (ATAC-seq) + GSE220545 (RNA-seq) for a senescence sequence-to-activity fine-tune

**Date:** 2026-08-20
**Scope:** GEO metadata (SOFT/series-matrix) and processed supplementary files only. No SRA/FASTQ access.
**Reproducibility:** all commands live in `scripts/`, run in numeric order; outputs land in `data/` and `output/`.

---

## Headline finding (read this first)

**GSE254358 does not exist.** It is not a private/embargoed record and not a path-construction
mistake — it is absent from the GEO database entirely. This was confirmed three independent ways:

| Check | Result |
|---|---|
| `esearch db=gds term=GSE254358[ACCN]` | `count: 0` |
| `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254358/` | HTTP 404 |
| `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/` (parent range folder) | HTTP 200 — the range folder is real, so this isn't a path bug |

I also checked GSE220545's BioProject (`PRJNA910178`) for a linked companion series and found
none — GSE220545 is the only GEO record under that BioProject, with no SuperSeries/relations
pointing anywhere else. And I checked six plausible single-digit-transposition typos of
`GSE254358` (`GSE254538`, `GSE245358`, `GSE253458`, `GSE254385`, `GSE245538`, `GSE253548`); all
resolve to real but unrelated series (mouse myocarditis, MLL-AF9 leukemia, U87MG TMZ resistance,
osteosarcoma FFPE — none is ATAC-seq, senescence, or IMR90). That's not proof no typo exists, but
it rules out the simple cases. **Confidence: high.**

**Practical consequence:** every ATAC-seq-side deliverable in this audit (Tasks 1, 3, 4 for the
ATAC arm) is unexecutable — there is no metadata and no supplementary file to pull. Task 2's
pairing question is answered trivially: the two series cannot be "genuinely paired" or even
"two independent experiments to stitch together," because only one of them is a real experiment.
**If you have a corrected accession for the ATAC-seq series, send it and I'll re-run this audit
against it** — the RNA-seq-side pipeline below is fully reusable, and the ATAC pipeline is
scripted but was never exercised against real data (see `scripts/04_atac_status.py`).

Everything below this point covers what *could* be audited: GSE220545 (RNA-seq), which is real
and reachable, but whose actual design is narrower than the "senescent vs proliferating, matched
to the ATAC series" framing implies (see Task 2).

---

## Task 1 — Metadata extraction

Parsed directly from `GSE220545_family.soft` (hand-rolled SOFT parser, `scripts/02_parse_metadata.py`
— no GEOparse dependency, no HTML scraping). Output: `output/tables/sample_metadata.csv`.

| GSM | Replicate | Timepoint | Cell line | Genotype | Treatment | Status (per protocol cross-ref) |
|---|---|---|---|---|---|---|
| GSM6806685 | A | Day 1.5 (30h) | IMR90 | ER-RAS | 4OHT | uncommitted — RAS-activated, NOT senescent, NOT an untreated control |
| GSM6806686 | B | Day 1.5 (30h) | IMR90 | ER-RAS | 4OHT | uncommitted |
| GSM6806687 | C | Day 1.5 (30h) | IMR90 | ER-RAS | 4OHT | uncommitted |
| GSM6806688 | A | Day 3 (72h) | IMR90 | ER-RAS | 4OHT | committed/senescent |
| GSM6806689 | B | Day 3 (72h) | IMR90 | ER-RAS | 4OHT | committed/senescent |
| GSM6806690 | C | Day 3 (72h) | IMR90 | ER-RAS | 4OHT | committed/senescent |

- **Cell line/type:** IMR90 (normal human diploid fibroblast), single line, single genotype
  (ER-RAS fusion construct) throughout.
- **Induction method:** oncogene-induced senescence (OIS) via 4-hydroxytamoxifen (4OHT)
  activation of a conditional ER-RAS<sup>G12V</sup> transgene — not replicative, not
  doxorubicin/therapy-induced, not irradiation.
- **Timepoint:** two levels only, 30h and 72h post-4OHT.
- **Passage number:** **not reported** anywhere in the SOFT record.
- **Donor:** **not reported** — IMR90 is a named cell line/single lot, not a donor cohort; GEO
  gives no donor/lot/passage identifier to check batch effects against.
- **Replicate identifier:** explicit in `Sample_title` ("Replicate A/B/C").

**Important parsing note (flagged, not silently resolved):** GEO gives this series no
per-sample senescence-status field. The only place "committed" vs "uncommitted" is defined is one
series-level sentence in `Sample_treatment_protocol_ch1`, copied identically onto all six samples
("Uncommitted control cells were treated ... 30h and committed cells were treated ... 72h").
A naive keyword search against that field would call *every* sample "ambiguous," because the
word "committed" is literally a substring of "uncommitted" and the sentence mentions both
regardless of which arm the sample belongs to. The script excludes that boilerplate field from
per-row status inference and instead encodes the 30h→uncommitted / 72h→committed mapping as an
explicit, hand-verified cross-reference (see `MANUAL_STATUS_CROSSREF` in
`scripts/02_parse_metadata.py`). This is the kind of thing that would silently mislabel a
training set if GEOparse or a naive characteristics-parser were trusted blindly.

**Assumption flagged:** GSE220545's only supplementary file
(`GSE220545_72h_vs_30h_DESEQ.txt`) has no per-sample GSM column headers — it has replicate-letter
columns (`A30h`, `B30h`, ... `C72h`). I mapped these to GSMs by matching replicate letter + timepoint
against `Sample_title` (e.g. "Replicate A" + "Day 1.5" → `A30h` → GSM6806685). This is the only
internally-consistent mapping (3 letters × 2 timepoints = 6 columns = 6 GSMs) but it is not stated
explicitly in the GEO record, so treat it as an assumption, not a GEO-confirmed fact.

---

## Task 2 — Design adequacy

Full log: `output/tables/design_adequacy_log.txt` (script: `scripts/05_design_adequacy.py`).

### Replicate counts

| Series | Condition | n | Flag |
|---|---|---|---|
| GSE220545 | Day 1.5 / uncommitted (30h) | 3 | OK (≥3) |
| GSE220545 | Day 3 / committed (72h) | 3 | OK (≥3) |
| GSE254358 | — | 0 | **Series does not exist** |

GSE220545 clears the ≥3-biological-replicates bar for both arms it has. But it only has
**two** conditions, and — this is the important part — **neither is a naive proliferating
baseline**. Both arms are 4OHT-treated, RAS-activated IMR90 cells; the difference between them is
30h vs 72h of RAS signaling, not "senescent" vs "never-induced, cycling fibroblast." The series'
own summary is explicit about this: it studies cells "not yet senescent but committed towards
that cell fate," i.e. a *commitment time-course*, not a senescent-vs-proliferating case-control
design. If the fine-tuning task genuinely needs a clean senescent-vs-proliferating contrast, this
series' 30h arm is a weak proxy for "proliferating" — it's RAS-activated at 30h, not vehicle-treated
or untreated.

### Are the two series paired?

No — because there's only one series. GSE254358 doesn't exist (see headline finding), and
GSE220545 is not part of any SuperSeries and has no companion series under its BioProject
(`PRJNA910178`, single GEO record: itself). There is no shared sample ID, donor, or BioProject to
check because there is nothing on the other side to check against. **Confidence: high** — this
rests on direct existence and relational queries against NCBI's own database, not inference.

---

## Task 3 — Processed data QC

### RNA-seq (GSE220545)

Only one supplementary file exists for this series: `GSE220545_72h_vs_30h_DESEQ.txt`
(DESeq2 output for the 72h-vs-30h contrast). GEO holds **no standalone per-sample counts/TPM
matrix** for this series — `Sample_supplementary_file_1` is `NONE` on every GSM. The DE table's
last six columns, however, carry DESeq2 size-factor-normalized per-sample counts
(`A30h, B30h, C30h, A72h, B72h, C72h`), so a usable expression matrix was reconstructed from
those columns (58,219 genes × 6 samples). **This is a reconstruction, not a deposited matrix** —
flagged as an assumption above.

| GSM | Column | Status | Normalized library size | Detected genes (>0 counts) | % of 58,219 genes |
|---|---|---|---|---|---|
| GSM6806685 | A30h | uncommitted | 13,734,790 | 20,563 | 35.3% |
| GSM6806686 | B30h | uncommitted | 13,863,010 | 20,243 | 34.8% |
| GSM6806687 | C30h | uncommitted | 13,942,680 | 20,364 | 35.0% |
| GSM6806688 | A72h | committed | 13,810,530 | 20,472 | 35.2% |
| GSM6806689 | B72h | committed | 13,692,510 | 20,926 | 35.9% |
| GSM6806690 | C72h | committed | 13,896,860 | 20,523 | 35.3% |

Library sizes are tight across all six samples (13.69M–13.94M, <2% spread) and detected-gene
counts are stable (20,243–20,926). No sample stands out as an outlier by depth or gene-detection
rate. **Caveat:** "library size" here is the sum of already-normalized counts, not raw sequencing
depth — GEO provides no raw read counts or FASTQ-derived depth for this series, and per the task
constraints SRA was not queried to get real depth.

### ATAC-seq (GSE254358)

Not executable — no series, no peak calls, no count matrix, no per-sample anything. See headline
finding.

---

## Task 4 — Signal check

Full numeric log: `output/tables/rna_qc_signal_log.txt`. Figure: `output/figures/rna_pca_and_correlation.png`.
Script: `scripts/03_rna_qc_signal.py`.

### RNA-seq PCA / correlation

- CPM-normalized (on top of the already DESeq2-normalized counts) → log2(CPM+1) → PCA on
  20,710 genes detected in ≥3/6 samples.
- **PC1 explains 80.6% of variance and cleanly separates committed (72h) from uncommitted (30h)**
  samples — mean PC1 separation between groups is ~74 units against a mean within-group spread of
  ~5 (ratio ≈ 14). PC2 (7.2%) and PC3 (5.5%) show no clear condition structure.
- Replicates cluster tightly within condition (see scatter plot) — no swapped or outlier sample.
- Spearman correlation on the same log2-CPM matrix: mean within-30h ρ = 0.976, mean within-72h
  ρ = 0.978, mean between-group ρ = 0.959. Real but modest separation in correlation space — that's
  expected for two timepoints of the same continuous induction process rather than two
  qualitatively distinct states.

**Interpretation:** the transcriptomic signal is real, reproducible, and dominates the variance —
this is a clean, well-executed RNA-seq experiment for the contrast it actually measures (early vs
late RAS-activation commitment). It is not, on its own, evidence of a "senescent vs proliferating"
signal, because — again — there is no proliferating arm in this series.

### ATAC-seq

Not executable — see headline finding.

### Senescence marker check

| Gene | Expected | log2FC (72h vs 30h) | Fold | padj | Observed | Match? |
|---|---|---|---|---|---|---|
| CDKN2A (p16) | up | **−0.59** | 0.66× | 0.004 | **DOWN** | **✗ MISMATCH** |
| CDKN1A (p21) | up | +0.75 | 1.68× | 2.3e-08 | UP | ✓ |
| LMNB1 | down | −1.18 | 0.44× | 1.8e-20 | DOWN | ✓ |
| MKI67 | down | −1.40 | 0.38× | 1.6e-06 | DOWN | ✓ |
| IL6 | up | +3.93 | 15.2× | 6.8e-20 | UP | ✓ |
| CXCL8 (IL8) | up | +4.28 | 19.4× | 2.6e-21 | UP | ✓ |

**5/6 markers move as expected. CDKN2A (p16) moves the wrong way — flagged prominently, as
instructed.** Two explanations are consistent with what's in the metadata, and I can't
distinguish them from processed data alone:

1. **Biologically plausible, not a misannotation.** p21 (CDKN1A) is a fast, p53-driven response to
   RAS/DNA-damage signaling and typically rises within hours to a few days; p16 (CDKN2A)/RB-pathway
   engagement is classically slower and is often what makes senescence *irreversible* rather than
   what initiates arrest. At 72h the series' own summary says cells are "committed towards" — not
   yet at — full senescence. A p21-up/p16-not-yet-up pattern at this early timepoint would be
   textbook-consistent with that framing, not a red flag.
2. **Sample-annotation risk.** Given how fragile the 30h/72h → uncommitted/committed mapping was
   to reconstruct (Task 1 note), a swapped timepoint label would also produce exactly this kind of
   single-marker mismatch against an otherwise-clean 5/6 signal. I have no independent evidence for
   a swap (PCA/correlation both cleanly separate on the labeled groups, and 5 other markers agree
   with the stated direction), but it can't be ruled out from processed data alone.

Given that MKI67, LMNB1, IL6, CXCL8, and CDKN1A all move correctly and the PCA separates cleanly
along the labeled axis, explanation (1) is more likely than (2) — but flagging it prominently per
the task's own instruction, since a single-marker mismatch on the canonical senescence gatekeeper
gene is exactly the kind of thing that should give a model-training pipeline pause, not be
smoothed over.

---

## Task 5 — Verdict

**Too thin, on two independent grounds — one fatal, one design-level.**

1. **Fatal: no ATAC-seq data exists.** GSE254358 is not a real accession. A "sequence-to-activity"
   fine-tune on top of Borzoi/gReLU needs chromatin accessibility (or another quantitative
   regulatory readout) as a training signal — RNA-seq alone doesn't supply that. Until a real,
   resolvable ATAC-seq accession is provided, there is no accessibility arm to fine-tune against
   at all. This isn't a depth or replicate problem; it's a missing dataset.

2. **Design-level, even setting (1) aside:** GSE220545, the one series that does exist, is not a
   senescent-vs-proliferating experiment. It's a 2-timepoint OIS commitment time-course (30h vs
   72h post-RAS-activation) with no untreated/vehicle proliferating control. Six samples, 3
   replicates per arm — replicate depth is fine — but the *contrast itself* doesn't match the
   "senescent vs proliferating" framing the task specifies. A model trained on this contrast would
   learn "early vs late RAS-pathway commitment in one fibroblast line," which is related to but
   narrower than general senescence-vs-proliferation biology, and it comes from a single cell line
   / single induction method / single genotype, so there's no route to check whether a learned
   effect generalizes beyond ER-RAS-driven OIS in IMR90.

**What's specifically missing, if this is to be salvaged:**
- **A real ATAC-seq accession** — re-run `scripts/01–05` against it once you have one; the
  RNA pipeline in this repo is a template for what to check (peak counts/widths, replicate
  overlap/Jaccard, CPM+log PCA, marker-gene direction check).
- **A true proliferating/untreated control arm** for the RNA-seq side (or a different RNA-seq
  series that has one) — without it, "senescent vs proliferating" is being approximated by
  "72h RAS-active vs 30h RAS-active," which is a materially different biological contrast.
- **Induction-method and cell-type diversity** — this data covers exactly one context (OIS,
  ER-RAS, IMR90). If the model is meant to generalize across senescence triggers (replicative,
  therapy-induced, irradiation) or cell types, none of that is represented here.
- **Matched pairing** — even if a real GSE254358-equivalent is found, it would need to be checked
  for the same cell line, comparable passage/timepoint, and ideally shared donor material before
  treating the two omics layers as a matched multi-modal pair rather than two datasets stitched
  together after the fact.

---

## Assumptions log (consolidated)

1. `GSE220545_72h_vs_30h_DESEQ.txt` column letters (A/B/C) map 1:1 to `Sample_title` "Replicate
   A/B/C" by timepoint — internally consistent, not GEO-confirmed.
2. That file's `A30h…C72h` columns are DESeq2 size-factor-normalized counts, not raw counts;
   "library size" and detected-gene counts in Task 3 are computed on that normalized scale.
3. `ko.sd` / `wt.sd` column headers in the same file look like unedited leftovers from a generic
   DESeq2 template (knockout/wild-type) rather than series-specific labels — not used in any
   computation, noted for transparency only.
4. Senescence status (committed/uncommitted) for GSE220545 was cross-referenced by hand from the
   series-wide `Sample_treatment_protocol_ch1` sentence, since GEO has no per-sample status field
   for this series — see Task 1 note.
5. GSE254358 was treated as genuinely nonexistent after three independent negative checks
   (eutils, FTP, BioProject cross-reference) plus a typo sweep — not merely "temporarily
   unreachable."

## Reproducing this audit

```
bash scripts/01_fetch_geo_metadata.sh   # pulls SOFT/matrix/suppl from GEO FTP for both accessions
python3 scripts/02_parse_metadata.py    # -> output/tables/sample_metadata.csv
python3 scripts/03_rna_qc_signal.py     # -> output/tables/rna_*.csv, output/figures/rna_pca_and_correlation.png
python3 scripts/04_atac_status.py       # documents why the ATAC arm is not executed
python3 scripts/05_design_adequacy.py   # -> output/tables/design_adequacy_log.txt
```

Requires `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, and `curl` on `PATH`. No
GEOparse dependency (SOFT parsed directly — see `scripts/02_parse_metadata.py` docstring for why).
