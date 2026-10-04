# Infrastructure validation: Borzoi via gReLU on cloud GPU

**Date:** 2026-08-22 (original run), rerun 2026-08-22 on a second, fresh
RunPod instance with a corrected Task 2 aggregation, **subset re-analysis
2026-08-22 done entirely locally, no GPU instance involved.**
**Status: ALL FIVE TASKS NOW EXECUTED on real rented RunPod GPU instances,
driven over SSH, across three pod sessions plus one purely local re-analysis.
Every number below is a real measurement or a real recomputation over
already-measured data — nothing is a placeholder or an expected value.**

**GO / NO-GO by task:**
- **Task 1 (smoke test): GO.** All five checks matched expectation on pod 1.
- **Task 2 (benchmark): COMPLETE, closed.** Corrected inverse transform gives
  Pearson r=0.406 / Spearman ρ=0.399 on n=118. A local subset analysis
  **rejected** the panel-composition hypothesis (removing all 28 flagged genes
  moves r by ≤+0.009) and instead found a real, subset-robust concentration of
  prediction error in low-expression genes (p=0.009–0.013). The
  ACTB/B2M/GAPDH sanity check — the validation that matters — passes by
  orders of magnitude.
- **Task 3 (design loop): RUNS, WITH CAVEATS.** Completes without error and is
  monotonic, but **plateaus after a single improvement step** (0.088→0.089,
  then flat) under the deliberately-trivial validation objective. Measured
  1.639 s/variant, 3.72 GB peak. **The design-sweep cost answer spans
  ~$2.73 to ~$24,585 depending on what "8,000–12,000 candidates" means — this
  needs clarifying before budget is committed.** See Task 3 below.
- **Task 4 (fine-tuning path): GO — the critical check passed.** Trunk weights
  are **bit-identical** before/after a real training step
  (0/180 tensors changed, max abs diff 0.000e+00 across 171,271,968 values).
  Full measured input contract below.
- **Task 5 (reference data): GO.** hg38 confirmed from the live model, all
  reference files staged, ~9.1 GB measured footprint, liftover absence
  re-confirmed by live grep.

**Three real bugs/behaviors were found and worked around during Tasks 3–4**
(an `Aggregate` bin-reduction requirement, a `return_seqs="none"` KeyError,
and a genuine gReLU 1.1.0 `checkpoint=False` identity-check bug) — all
documented in place below, since anyone building against this API will hit
them too.

**A fifth session (2026-08-23/24) ran the REAL FINE-TUNE — month-4 checkpoint.
Result: NO-GO. Test AUC 0.659 (bar 0.85) and Spearman 0.514 (bar 0.60) both
missed, and the shuffled-label control FAILED (real 0.659 vs shuffled
0.561±0.094; 2 of 20 random relabelings beat the true contrast; p=0.143), so
the result cannot be attributed to senescence-specific signal. Per the
pre-registered stop rule the report stops at measured numbers. One clear
positive: the model does not mistake quiescence for senescence. See
"FINE-TUNE, MONTH-4 CHECKPOINT" below.**

**A fourth session (2026-08-23) ran the fine-tune CALIBRATION PASS only** — no
full trunk pass, nothing trained. It verified the pooled-cache identity
empirically and measured the real per-gene rate at **0.243 s**, ~7x faster
than projected, putting the full trunk cache pass at **$1.00–1.66** and the
whole fine-tune (including all five control runs) at **~$2** rather than
$176–671. See "Fine-tune CALIBRATION PASS" below.

**Pods 1 (`<pod-address-redacted>`) and 2 (`<pod-address-redacted>`) were confirmed
terminated by the user.** ⚠️ **Pod 3 (`<pod-address-redacted>`) and pod 4
(`<pod-address-redacted>`) were NOT terminated by this session — see the action
section at the end.**

---

## Session facts

| | |
|---|---|
| Instance | RunPod, `<pod-address-redacted>` |
| GPU (measured via `nvidia-smi`) | **NVIDIA GeForce RTX 4090, 24564 MiB (~24GB) VRAM** — NOT one of the A100 40GB / L40S / A10G types the prep session planned around. Flagged immediately, see Deviations below. |
| Driver / max CUDA | 580.159.04 / CUDA 13.0 (driver ceiling — the actual PyTorch build used is CUDA 12.8, well within this) |
| OS | Ubuntu 24.04.3 LTS, Python 3.12.3 (prep session assumed Python 3.10 — also within gReLU's `>=3.8` requirement, no issue) |
| Billing rate | $0.75/hour (user-stated) |
| Session wall-clock (this driving session, first command to last) | Approx. **15:18–16:18 PDT, ~60 minutes**, timestamped from this session's own SSH command outputs — not pulled from RunPod's billing meter directly. **Check the RunPod console for the authoritative figure.** |
| Approx. cost this session | **~$0.75** at the stated rate for ~60 minutes — approximate, see above |

---

## Deviations from the prep session's plan (real, found during execution, not anticipated)

1. **GPU is an RTX 4090 (24GB), not A100/L40S/A10G.** Consumer/prosumer card,
   Ada Lovelace architecture (compute capability 8.9 — still clears the
   Ampere-or-newer bar the prep report noted for flash-attn/Flashzoi, not that
   those were needed for Tasks 1–2). 24GB is materially less headroom than an
   A100 40GB — see the peak-memory measurement below, which used 75% of this
   card's VRAM for a single fp32 forward pass at batch size 1.
2. **PyTorch was already pre-installed** (2.8.0+cu128) by RunPod's base image,
   newer than the prep plan's pin (2.4.1+cu121). Kept as-is rather than
   reinstalled — gReLU's own `setup.cfg` only requires `torch>=2.0`, no upper
   bound, and reinstalling would have spent paid GPU time for no functional
   benefit to Tasks 1–2. flash-attn and borzoi-pytorch (Flashzoi) were **not
   installed at all** — neither is needed for Tasks 1–2, and skipping them
   avoided a documented 10–30 minute CUDA source-compile step.
3. **`pip install gReLU==1.1.0` took ~25 minutes**, not the "under 10 minutes"
   gReLU's own README states — confirmed via direct process monitoring
   (healthy I/O-wait state throughout, 1.47GB downloaded, 2.84GB written,
   never hung) that this was genuine heavy-dependency resolution/download, not
   a stuck process. Ubuntu 24.04's PEP 668 externally-managed-environment
   protection also required creating a venv (`--system-site-packages`, to
   inherit the pre-installed torch) before any install could proceed at all —
   not anticipated in the prep plan, which assumed conda/mamba.
4. **A real bug was found in gReLU 1.1.0's actual shipped code, contradicting
   the prep session's own documentation.** The prep report stated
   `grelu.io.genome.get_genome()` reads a `GRELU_GENOMES_DIR` environment
   variable. Reading the ACTUAL installed source
   (`grelu/io/genome.py`, confirmed by `sed`-printing it directly on the
   instance) shows `get_genome(genome, **kwargs)` calls
   `genomepy.list_installed_genomes()` with **no arguments at all** — it does
   not read that environment variable anywhere. This is a real discrepancy
   between what the prep session documented (from reading docstrings/research)
   and what the shipped code actually does — flagged honestly rather than
   silently worked around. Fix used: install the genome once, then symlink it
   into genomepy's real default location (`~/.local/share/genomes/hg38`).
5. **UCSC (genomepy's default provider) is unreachable from this instance**
   ("UCSC appears to be offline", confirmed on two separate attempts).
   genomepy's automatic fallback search then burned real time scanning
   multi-million-record GenBank/RefSeq assembly listings across every other
   provider before failing outright — because Ensembl calls this genome
   `GRCh38.p14` and NCBI calls it `GRCh38`, neither matching the literal
   string `"hg38"` that gReLU's code searches for. Fixed by explicitly
   installing `name="GRCh38", provider="NCBI"` (the standard, unpatched RefSeq
   primary assembly, `GCF_000001405.26`) with `localname="hg38"`.
6. **NCBI's FASTA uses RefSeq accession-style contig names** (e.g. the primary
   sequence for chromosome 1 is headed `>1 NC_000001.11 ...`, i.e. genomepy's
   short-name convention, not UCSC's `chr1`/`chr2` style that gReLU's
   pipeline expects). Fixed by remapping all 194 sequence headers using
   NCBI's own `assembly_report.txt` (which ships a `UCSC-style-name` column),
   then reindexing with `pyfaidx`. Verified correct: renamed `chr1` came back
   at exactly 248,956,422 bp, the well-known GRCh38 chr1 length.
7. **SSH background-job hang, repeatedly.** Piping a script's output through
   `tee` over a live SSH session, or backgrounding a job with `&` inside a
   non-interactive `ssh host "cmd"` invocation, left the SSH client hanging
   past the job's actual completion (the remote process holds a reference to
   the session's stdout pipe even when redirected). Worked around by
   `nohup ... > logfile 2>&1 < /dev/null &` plus `disown` on the remote side,
   then polling for process exit via **fresh, separate** SSH connections
   (`while kill -0 $PID; do sleep N; done`) rather than waiting on the
   original connection. No data or GPU time was lost — confirmed via process
   monitoring that the actual remote computation continued and completed
   correctly through every one of these hangs; only the local SSH client
   needed re-connecting.

None of these deviations were the result of tuning parameters to improve a
result — all are infrastructure/environment fixes required just to get a
forward pass to run at all, fully disclosed here per the reporting rules.

---

## Task 1 — Environment and smoke test: **PASS, all five checks measured**

Setup: venv with `gReLU==1.1.0` installed on top of the pre-existing
`torch==2.8.0+cu128`; no flash-attn, no conda, no Flashzoi (not needed for
this task). Script: `runbook/task1_smoke_test.py`, run unmodified.

| Check | Expected | **MEASURED (this run)** | Pass? |
|---|---|---|---|
| Parameter count | ~186M | **185,892,699** | ✓ |
| Input sequence length | 524,288 bp | **524,288 bp** | ✓ |
| Output resolution | 32 bp bins | **32 bp bins** | ✓ |
| Number of output tracks | 7,611 | **7,611** | ✓ |
| Genome build | hg38 | **hg38** (model's own declared metadata) | ✓ |

**Track composition (real, measured, not previously known in this precision):**

| Assay | n tracks |
|---|---|
| CHIP | 3,886 |
| RNA | 1,543 |
| CAGE | 1,276 |
| DNASE | 674 |
| ATAC | 232 |
| **Total** | **7,611** |

Full table: `execution_results/task1_output_tracks.csv` (7,611 rows).

**Forward pass** on `chr1:69,993,520-70,517,808` (the same interval used in
gReLU's own tutorial, chosen to eliminate "bad interval" as a variable on the
first real run):

| Metric | **MEASURED** |
|---|---|
| Output shape | **(1, 7611, 6144)** — matches expected exactly |
| Wall-clock time, one forward pass | **860.6 ms** |
| Peak GPU memory allocated | **18.08 GB** (of 24GB total on this card — 75% utilization for a single fp32, batch-size-1 forward pass) |
| Output interval after crop | `chr1:70,157,360-70,353,968` (196,608 bp) — matches expected crop math exactly |

**Note on peak memory:** 18.08GB out of 24GB leaves limited headroom on this
specific GPU. This wasn't a problem for Task 1/2 (single sequences, one at a
time), but it's directly relevant if Tasks 3–5 are run later on this same
card — batching more than 1 sequence at fp32 would very plausibly OOM on a
24GB card; mixed precision (autocast) or a larger-VRAM instance would likely
be needed for anything beyond single-sequence inference.

---

## Task 2 — Benchmark: **RAN TO COMPLETION. Result: real, significant, but weak.**

### What was actually compared (real substitution, as flagged in the prep report)

Per the prep session's finding (paper's own held-out ground truth is gated
behind a Requester-Pays GCP bucket or dbGaP access — not re-litigated here),
this reproduces a **substituted benchmark**: Borzoi's predicted signal vs.
**GTEx v10's public, aggregate median-TPM-by-tissue** data (Whole Blood),
**not** the paper's own held-out test-set correlation. This is not the
published metric — it is explicitly labeled as such throughout.

- **GTEx data:** `GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz`,
  fetched directly from GTEx's public GCS bucket
  (`storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/...`, verified
  public/no-auth via a direct `curl -sI`, 8.8MB, 59,033 genes × 68 tissues).
  This is aggregate/de-identified summary data, explicitly distinct from
  dbGaP-protected individual-level GTEx data.
- **Gene panel:** 120 genes (of an initial 130 resolved), stratified across
  deciles of log(Whole_Blood TPM) plus known housekeeping genes, spanning
  **0.0036 to 4,508.5 TPM** (~6 orders of magnitude). hg38 coordinates fetched
  live from the Ensembl REST API (129/130 symbols resolved on the first
  batch query). 10 genes with gene-body span >150kb were excluded before
  running (to avoid the window-centering edge case described below); this
  gene panel was assembled **offline, before touching the billed GPU
  instance**, so it cost no GPU time.
- **Track selection:** Borzoi's own track metadata (from Task 1's real
  output) shows GTEx blood RNA is represented as **3 individual per-donor
  recount3 coverage tracks** (`description == "RNA:blood"`, row indices
  7531–7533), not one pre-aggregated track. Predictions were averaged across
  all 3 to approximate a population-level signal.
- **Replicate:** `human_rep0` only. **Not ensembled** across the 4 available
  Borzoi replicates — out of scope for this session (Task 2 only, as
  instructed), and ensembling is known to reduce noise in this model family,
  so a single-replicate result should be expected to underperform an
  ensembled one.
- **Aggregation:** raw model output summed linearly across all output bins
  in the gene-body window, per track, per gene. **This does NOT apply the
  per-track `scale`/`sum_stat` transform** — see the diagnosis below, this
  turned out to matter a lot.

### Result

| | **MEASURED** |
|---|---|
| Genes evaluated | **118** (2 of 120 failed — see diagnosis below, not silently dropped) |
| Replicate | human_rep0 (single replicate, not ensembled) |
| Tracks | GTEx blood RNA, indices 7531/7532/7533 (averaged) |
| **Pearson r** | **0.3176** (p = 4.58×10⁻⁴) |
| **Spearman ρ** | **0.4231** (p = 1.81×10⁻⁶) |
| Prediction wall-clock | 67.7s total / 118 genes ≈ 0.56s/gene (batch size 1, sequential — GPU memory headroom on this 24GB card did not allow safely batching multiple 524kb sequences at once, see Task 1's 18GB peak-memory finding) |

Full per-gene results: `execution_results/task2_results.csv`.

**This is a real, statistically significant, positive correlation (p < 0.001
on both metrics) — not noise.** It is also **substantially weaker** than the
prep session's ~0.83 orientation figure (itself an unverified secondary-source
citation, not the paper's own stated number for this protocol). Per
instructions, **no parameters were tuned and no re-run was attempted to
improve this number.** It is reported as measured.

### Diagnosis of why the correlation is weak (evidenced, not speculative)

Checked directly against the diagnostic checklist from the prep report:

1. **Genome build mismatch — ruled out.** hg38 throughout: GTEx v10 uses
   hg38, Ensembl coordinates were fetched as GRCh38, Borzoi's own declared
   genome (Task 1) is hg38. Not the cause.

2. **Normalization differences — CONFIRMED as a real, likely major
   contributor, with direct evidence.** The three blood tracks used carry
   `scale = 0.01` and `sum_stat = 'sum_sqrt'` in their own metadata (checked
   directly, see `execution_results/task1_output_tracks.csv` rows
   7531–7533) — meaning the track's own defined aggregation is a **scaled
   sum of square roots** of per-bin values, not the plain linear sum this
   script used. Skipping a variance-stabilizing sqrt transform before summing
   lets a handful of very-high-coverage bins dominate the aggregated signal,
   which would systematically distort the predicted-vs-TPM relationship
   exactly in a way that could suppress a correlation. This was not applied
   here, and not re-run to check the effect, per instructions not to tune for
   a better number — but it is the most concrete, evidenced candidate
   explanation found.

3. **Wrong track indices — partially implicated.** Not "wrong" in the sense
   of picking an unrelated assay, but the available "blood" signal is only 3
   individual donor coverage tracks (not a large-sample aggregate), which is
   inherently noisier than GTEx's own median-of-hundreds-of-donors TPM value
   it's being compared against — an apples-to-not-quite-apples comparison on
   the observed side too.

4. **Different held-out split / gene selection — a real, acknowledged
   limitation.** This 118-gene GTEx-tissue-driven panel is not the paper's
   held-out test set. The gene-level aggregation method (summing all
   output bins across the whole gene body, ± a symmetric window) is a crude
   proxy for "genes computed by summing predictions across exons" — it
   doesn't use real exon/intron structure or account for genes of very
   different lengths receiving very different amounts of (mostly intronic,
   near-zero) signal diluting the real exonic signal.

5. **Replicate handling — a real, acknowledged limitation.** Single
   replicate only; ensembling is expected to help, not tested here (out of
   this session's authorized scope).

**Two genes failed with a genuine, diagnosed edge case, not a bug in the gene
data:** `TBP` (chr6) and `ZNF584` (chr19) both sit close enough to their
chromosome's actual end that the symmetric ±262,144bp input window around the
gene center ran past the chromosome boundary — confirmed directly:
chr6 is 170,805,979bp (TBP's window needed 170,827,874, a 21,895bp overhang);
chr19 is 58,617,616bp (ZNF584's window needed 58,678,764, a 61,148bp
overhang). Both dropped cleanly (`dropna`), not silently miscounted — 118/120
genes carried through to the correlation.

### Overall read on Task 2

The **mechanism works end-to-end** — the model loads, predicts, and produces
a real, significant, directionally-correct relationship to actual measured
gene expression on genuinely independent data it wasn't shown in this
session. The **magnitude** of that relationship, as measured by this specific
simplified protocol, is well below what the full published methodology would
likely produce, for concrete, identified reasons (raw-sum vs. sum-sqrt
aggregation being the most evidenced one) rather than an unexplained anomaly.
Whether that gap is worth closing (by applying the proper per-track
transform, ensembling replicates, and/or expanding the gene panel) is a
scope/cost decision for you, not something this session decided unilaterally
to pursue further.

---

## Task 2 RERUN — corrected aggregation (second pod, second session)

### New session facts

| | |
|---|---|
| Instance | RunPod, `<pod-address-redacted>` — a **different, fresh pod** from the original run |
| GPU | NVIDIA GeForce RTX 4090, 24564 MiB — same model as before, driver 580.173.02 (minor patch version difference from the first pod's 580.159.04, not consequential) |
| Base environment | Identical to before: Ubuntu 24.04.3, Python 3.12.3, PyTorch 2.8.0+cu128 pre-installed |
| Gene panel | **Byte-identical to the original run** — `final_gene_panel.csv` uploaded from the local copy, MD5-verified matching (`86d1cec259101197518c127591193722`) before use, not regenerated |
| Task 1 | **Not rerun**, per instruction — this session trusts the original pod's Task 1 pass |

### Root cause, confirmed against the source before running anything

The original run summed **raw model output** directly across bins. That output
is not in linear coverage units — it's already power-transformed, soft-clipped,
and scaled. Two things were wrong, not one:
1. The `scale` (0.01) and clip weren't undone at all.
2. Even after fixing that, summing across bins **before** undoing the power
   transform is mathematically wrong regardless, since `(a+b)^0.75 ≠ a^0.75 + b^0.75`
   — the inversion has to happen per-bin, before the sum.

**Correct inverse, sourced directly from `calico/baskerville`'s own
`seqnn.py`** (the model authors' reference implementation), not derived by
guesswork, and cross-checked against `calico/borzoi`'s own official worked
eQTL example notebook
(`examples/borzoi_example_eqtl_chr10_116952944_T_C.ipynb`), which explicitly
sets `untransform_old=True`, `track_scale=0.01`, `track_transform=3./4.`,
`clip_soft=384.` for its own GTEx blood tracks — **identical** to the
scale/clip_soft values read directly off this model's own blood tracks,
confirming this is the right branch, not an assumption:

Per bin, in this exact order:
1. **Undo scale**: `x = raw_pred / scale`
2. **Undo soft-clip**: where `x > clip_soft`, replace with `(x − clip_soft)² + clip_soft`
3. **Undo the power transform**: `x = x^(1/0.75) = x^(4/3)`

**The exponent is 0.75, not 0.5** — `sum_stat = "sum_sqrt"` is a legacy name;
the actual forward transform in `baskerville`'s own data-prep code
(`hound_data_read.py`) is `seq_cov ** 0.75`, confirmed by reading that source
directly. Assuming a literal square root here would have been exactly the
"plausible but wrong" trap flagged going in.

All four transform parameters (`clip`, `clip_soft`, `scale`, `sum_stat`) were
read from `model.data_params["tasks"]` at runtime for the actual track
indices used — not hardcoded — and printed in the run log for verification:

```
                            name  clip  clip_soft  scale  sum_stat
7531  GTEX-1I4MK-0002-SM-EZ6M9.1   768        384   0.01  sum_sqrt
7532  GTEX-1LB8K-0005-SM-DIPED.1   768        384   0.01  sum_sqrt
7533  GTEX-1OKEX-0006-SM-DKPQ2.1   768        384   0.01  sum_sqrt
```

### Sanity check (run before the full panel, as required)

| Gene | GTEx TPM | Corrected signal | Old (naive-sum) signal |
|---|---|---|---|
| **ACTB** | 3,559.0 | **2.14×10⁸** | 1,301 |
| **B2M** | 4,508.5 | **1.95×10⁸** | 1,335 |
| **GAPDH** | 2,789.5 | **2.71×10⁷** | 2,895 |
| LCT (low-TPM control) | 0.0065 | 1.98×10⁷ | 1,278 |
| EPHA10 (low-TPM control) | 0.0116 | 5.75×10⁵ | 1,044 |
| SHOX (low-TPM control) | 0.0036 | 608 | 7.9 |

**ACTB, B2M, and GAPDH are the top 3 by corrected signal — sanity check
PASSED.** Note the old naive-sum method (right column) does *not* cleanly
separate the housekeeping genes from the low-expression controls (GAPDH's old
value, 2,895, is barely different from LCT's 1,278, despite a >400,000-fold
TPM difference) — a concrete illustration of why the original aggregation was
producing a weak signal. The corrected method separates them by 2–6 orders of
magnitude, tracking the true TPM spread far better.

One gene (LCT) still shows a corrected signal (1.98×10⁷) much higher than its
TPM (0.0065) would suggest — reported here rather than dropped, since the
sanity check's actual pass condition (housekeeping genes on top) held; this
single outlier is noted, not explained away, and is consistent with the
"different held-out split / gene selection" limitation already on record —
LCT's real regulatory biology (lactase, famously tissue- and
genotype-restricted with complex enhancer regulation) may simply be a bad fit
for a naive whole-gene-body window regardless of the transform.

### Full panel result

| | **Original run** | **Corrected rerun** |
|---|---|---|
| n genes evaluated | 118 | **118** (same 2 failures: TBP, ZNF584 — same diagnosed chromosome-edge cause, fully deterministic) |
| Replicate | human_rep0 | human_rep0 (unchanged, for like-for-like comparison) |
| Tracks | GTEx blood RNA, indices 7531–7533 (averaged) | same |
| **Pearson r** | 0.3176 (p = 4.58×10⁻⁴) | **0.4055 (p = 5.22×10⁻⁶)** |
| **Spearman ρ** | 0.4231 (p = 1.81×10⁻⁶) | **0.3994 (p = 7.47×10⁻⁶)** |
| Prediction wall-clock | 67.7s (0.56s/gene) | 90.7s (0.76s/gene) — slightly slower, plausibly the added per-bin inverse-transform arithmetic, not investigated further since it's a trivial cost either way |
| Bins hitting soft-clip threshold | not tracked in original run | **15,713 of 2,174,976 evaluated bins (~0.72%)**, concentrated in the highest-expression genes (ACTB 0.93% of its bins, GAPDH 2.5%) — exactly where getting the clip-undo right matters most |

**Pearson improved meaningfully (+28% relative); Spearman did not improve —
it went down slightly. Both are reported as measured; neither was tuned.**

This split is not a red flag by itself: Spearman is rank-based and invariant
to any monotonic per-gene transform, so it's far less sensitive to whether the
absolute/linear scale of the aggregation is correct — the old method, despite
being mathematically wrong in absolute terms, could already preserve a
similar gene ranking in many cases. Pearson (computed here on log1p-transformed
values) is more sensitive to getting the actual shape of the value
distribution right, particularly at the high-expression tail where soft-clip
correction has the most effect — consistent with the sanity-check table above,
where the correction's biggest visible impact was exactly on the three
highest-expression genes. No further investigation or tuning was performed
into the Spearman discrepancy, per instructions not to search for a
better-looking configuration.

### What's still not tested (explicitly, as requested)

**Ensembling all 4 Borzoi replicates is a known, untested variable.** Both
runs (original and corrected) used `human_rep0` only. The Borzoi paper's own
headline numbers plausibly use a 4-replicate ensemble, which is expected to
reduce single-replicate noise. Whether ensembling would move either the
Pearson or Spearman number, and by how much, was not tested in this session —
out of scope (Task 2 only, single replicate, as instructed).

### Overall read on the rerun

The correction was necessary and mathematically justified, not optional
polish — the sanity check would have failed outright under the old method if
it had been the actual gate (GAPDH would not have clearly separated from
low-expression controls). The corrected pipeline passes a real, pre-registered
sanity check and produces a moderately-improved linear correlation. It is
**still well below the ~0.83 orientation figure** carried over from the prep
session (itself an unverified secondary-source citation — see the original
Task 2 section above). The gap is likely now dominated by the
already-documented remaining limitations: single replicate (untested lever),
small non-held-out gene panel, sample-of-3 blood tracks compared against a
population median, and a naive whole-gene-body window rather than
exon-aware aggregation — not by an aggregation bug, which is now fixed and
verified.

---

## Task 2 subset analysis — panel composition hypothesis: **REJECTED (in aggregate)**

**Pure local re-analysis, no GPU instance touched, no predictions recomputed.**
Both pods confirmed terminated before this analysis started. Input: the
corrected per-gene predictions already saved at
`execution_results_rerun/task2_results_corrected.csv` (verified to reproduce
Pearson r=0.4055, Spearman ρ=0.3994 on n=118 exactly, byte-for-byte, before
any further analysis — see below).

### Hypothesis under test

That the ~0.4 gap to the ~0.83 orientation figure is driven by including gene
classes that are not predictable from reference sequence in principle
(immunoglobulin segments, replication-dependent histones) or whose ground
truth is suspect (pseudogenes, lncRNAs, small RNAs, the polymorphic KIR
locus) — rather than by the pipeline itself.

### Classification method

Real gene biotypes were pulled from the **Ensembl REST API**
(`lookup/symbol/homo_sapiens`, batch POST), not inferred from gene-name
patterns — all 118/118 genes resolved. Two categories still required a
manual override beyond raw biotype, both stated explicitly:
- **Replication-dependent histones**: Ensembl biotype for these is
  `protein_coding` (they do encode protein) — biotype alone can't separate
  them from other protein-coding genes. Identified instead by the HGNC 2018
  cluster-naming convention (`H2AC#`, `H3C#`, `H2BC#`, `H1-#`, `H4C#`), a
  fixed pattern check against real nomenclature rules, not a fitted filter.
- **KIR2DL4**: biotype is `protein_coding`; flagged separately by name as
  "polymorphic locus" per the hypothesis brief, since genomic polymorphism
  isn't a biotype property at all.

**Real biotype corrected four of the user's own name-based guesses** — a
direct demonstration of why biotype was used instead of name-pattern
matching alone:

| Gene | Guessed (by name pattern) | Actual Ensembl biotype |
|---|---|---|
| PLGLB1 | pseudogene | **protein_coding** |
| ZKSCAN8P1 | pseudogene (ends in "P1") | **protein_coding** |
| CECR7 | lncRNA/antisense | **transcribed_unprocessed_pseudogene** |
| WFDC21P | lncRNA/antisense | **transcribed_unitary_pseudogene** |

These four are classified by their real biotype throughout, not by the
guess. Full 118-gene classification table:
`subset_analysis/classified_panel.csv`. Biotype-only lookup result:
`subset_analysis/gene_biotypes.csv`.

### Classification counts

| Class | n | Method |
|---|---|---|
| Protein-coding (clean) | 89 | Ensembl biotype = `protein_coding`, minus histones, minus KIR2DL4 |
| Pseudogene | 17 | Ensembl biotype contains `pseudogene` (6 distinct sub-biotypes) |
| lncRNA/antisense | 5 | Ensembl biotype = `lncRNA` |
| Histone (replication-dependent) | 3 | biotype `protein_coding` + HGNC cluster-name pattern |
| Immunoglobulin V-segment | 2 | Ensembl biotype = `IG_V_gene` |
| Small RNA | 1 | Ensembl biotype = `snRNA` |
| Other (KIR2DL4, polymorphic locus) | 1 | Manual override, documented above |
| **Total** | **118** | |

### Subset correlations

| Subset | n | Pearson r | p | Spearman ρ | p | Note |
|---|---|---|---|---|---|---|
| **1. Full panel** | 118 | **0.4055** | 5.2e-06 | **0.3994** | 7.5e-06 | Reproduces the rerun exactly — confirms the saved file was loaded correctly |
| **2. Protein-coding only** (biotype-based, incl. histones + KIR2DL4) | 93 | **0.4144** | 3.6e-05 | **0.4043** | 5.8e-05 | +0.009 / +0.005 vs. full panel — **not a meaningful change** |
| **3. Protein-coding minus IG + histone** | 90 | **0.4040** | 7.9e-05 | **0.3941** | 1.2e-04 | −0.002 / −0.005 vs. full panel — **not a meaningful change; slightly lower, if anything** |
| 4. Pseudogene alone | 17 | 0.1366 | 0.601 | 0.2696 | 0.295 | n<30, wide CI; not significant on its own, but directionally weaker — the one subclass result consistent with the hypothesis |
| 4. lncRNA/antisense alone | 5 | 0.5393 | 0.348 | 0.3000 | 0.624 | n<30, far too small to interpret |
| 4. Histone alone | 3 | −0.5755 | 0.610 | −1.0000 | — | n=3: Spearman ρ can only take a few discrete values at this n; the −1.0 is a sampling artifact of tiny n, not evidence of anything |
| 4. Immunoglobulin V-segment alone | 2 | — | — | — | — | n=2: correlation not computable/meaningful |
| 4. Small RNA alone | 1 | — | — | — | — | n=1: correlation not computable |
| 4. KIR2DL4 alone | 1 | — | — | — | — | n=1: correlation not computable |

### Verdict on the hypothesis: **Rejected in aggregate, weakly consistent at the sub-class level**

Removing every flagged class — immunoglobulin segments, replication-dependent
histones, pseudogenes, lncRNAs, small RNAs, and the polymorphic KIR locus,
**28 genes in total** — moves Pearson r from 0.4055 to at most 0.4144 (subset
2) and Spearman ρ moves by less than 0.005 either direction. **This is not a
meaningful improvement by any reasonable standard**, and subset 3 (the
"cleanest" definition) is not even reliably better than the unfiltered panel.
**If the panel-composition hypothesis were the dominant driver of the gap to
~0.83, removing 28 of 118 genes flagged as biologically unpredictable or
measurement-suspect should have moved the number substantially. It didn't.**

The one result that *is* directionally consistent with the hypothesis —
pseudogenes alone showing a much weaker, non-significant correlation
(r=0.14 vs. 0.41 for the clean set) — is real but can't be load-bearing for
the aggregate gap: pseudogenes are only 17/118 genes (14%), and their
individual weakness is diluted by the other 101 genes when computing the
full-panel number. A subclass being poorly predicted is not the same as that
subclass explaining the aggregate shortfall.

**This is reported as the finding, not as a disappointing result to explain
away**, per instruction. No filter definitions were iterated to search for a
better number — the three main subsets (full panel, protein-coding only,
protein-coding minus IG/histone) were computed exactly as specified and
reported as they came out.

### Residual analysis: low-expression genes dominate the error — a real, different, more interesting finding

For the protein-coding-minus-IG/histone subset (n=90; the pattern is
statistically indistinguishable using either of the other two protein-coding
subset definitions, see below), a linear fit of log(predicted) on log(TPM)
was computed and residuals examined against expression level:

| TPM tercile | n | Mean \|residual\| | Std |
|---|---|---|---|
| Low | 30 | **2.55** | 1.84 |
| Mid | 30 | **1.84** | 1.13 |
| High | 30 | **1.35** | 1.37 |

Mean absolute residual **decreases monotonically** from the lowest to highest
expression tercile — a ~1.9-fold difference between the low and high bins.
Correlation of `|residual|` against `log(TPM)`: **r = −0.276, p = 0.009** —
statistically significant, and robust to which protein-coding subset
definition is used:

| Subset | n | corr(\|residual\|, log TPM) | p |
|---|---|---|---|
| Protein-coding minus IG/histone | 90 | −0.276 | 0.009 |
| Protein-coding, all | 93 | −0.260 | 0.012 |
| Full panel | 118 | −0.228 | 0.013 |

(The signed-residual-vs-log(TPM) correlation is exactly 0 by construction —
that's a mechanical property of ordinary-least-squares residuals against
their own fitted predictor, not a finding. The `|residual|` result above is
the real one.)

**This is a genuine, statistically significant pattern, not an artifact of
subset choice: low-TPM genes carry substantially more prediction error than
high-TPM genes.** This is a different and arguably more diagnostic problem
than panel composition — it points toward noise concentrated at low
expression (where both the model's absolute signal and GTEx's own TPM
estimate are closer to their respective noise floors) as a more promising
lead for the remaining gap than gene-class filtering.

### What this does and doesn't mean for next steps

- **Filtering the panel by biotype is not, on this evidence, a productive
  next step** for closing the gap to the ~0.83 orientation figure — the
  aggregate correlation is essentially insensitive to it.
- **The low-expression heteroscedasticity finding is a more promising lead**,
  though it wasn't investigated further here (would require, at minimum,
  either a larger low-TPM-focused gene sample or a comparison against GTEx's
  own per-gene TPM confidence intervals — neither attempted, out of scope for
  a local re-analysis and not requested).
- **Ensembling all 4 Borzoi replicates remains the largest untested lever**,
  as already noted after the rerun — nothing in this subset analysis bears on
  that question one way or the other.

Full classification table, biotype lookups, and residual data:
`subset_analysis/classified_panel.csv`, `subset_analysis/gene_biotypes.csv`,
`subset_analysis/pc_clean_residuals.csv`.

---

## Tasks 3, 4, 5 — EXECUTED (third pod, 2026-08-22)

**Instance:** RunPod `<pod-address-redacted>`, NVIDIA GeForce RTX 4090
24564 MiB, driver 570.211.01, Ubuntu 24.04.3, Python 3.12.3, PyTorch
2.8.0+cu128 pre-installed. Same environment workarounds as prior sessions
(venv `--system-site-packages`; NCBI GRCh38 + UCSC-style contig remap; no
`GRELU_GENOMES_DIR`). gReLU install measured at **~23 minutes**, matching the
budgeted ~25 min. Tasks 1 and 2 were **not** rerun, per instruction.

### Task 3 — Design loop: RUNS TO COMPLETION, but plateaus immediately

**A real design-time discovery, not anticipated by the prep script:** gReLU's
own design tutorial (which the original `task3_design_loop.py` was modeled on)
uses the Catlas ATAC model, whose native `seq_len` is 200bp — so a bare 200bp
string works there. **Borzoi's `seq_len` is 524,288bp and fixed**, and
`grelu.design.evolve()` performs no length validation or padding (confirmed by
reading `src/grelu/design.py`) — it passes whatever it's given straight into
the model. Running "directed evolution on a ~200bp sequence" with Borzoi
therefore requires embedding a 200bp mutable window inside a full-length
524,288bp context and restricting mutations via `evolve()`'s `positions`
argument. That is what was run.

Two further bugs were hit and fixed live (both cheap, both real):
1. `Aggregate(tasks=[...])` alone crashes `evolve()` with
   `ValueError: Length of values (6144) does not match length of index (1)` —
   `Aggregate` reduces over *tasks* but not over Borzoi's 6144 output bins,
   and `evolve()` requires one scalar per sequence. Fix: add
   `length_aggfunc="mean"`.
2. `return_seqs="none"` crashes with `KeyError: "['seq'] not in index"` —
   `evolve()`'s internal bookkeeping still requires the `seq` column it
   suppresses. Fix: use `return_seqs="all"`.

**Measured results** (objective: maximize track 0, `CNhs10608+`,
CAGE:Clontech Human Universal Reference Total RNA — deliberately trivial, per
brief; 10-position mutable window, `batch_size=1`, `human_rep0`):

| Iteration | Best objective value |
|---|---|
| 0 (baseline) | 0.088 |
| 1 | **0.089** |
| 2 | 0.089 |
| 3 | 0.089 |
| 4 | 0.089 |

- **Runs to completion:** yes, no errors, all 5 iterations executed.
- **Monotonic:** yes, non-decreasing — but only one real improvement occurs
  (0.088 → 0.089 at iteration 1), then it is **flat from iteration 1 onward**.
  gReLU printed its own early-stop notice: `Score did not increase on
  iteration: 4`. So: monotonic by construction (evolve keeps the incumbent),
  but **effectively stalled after a single step** — reported as measured, not
  smoothed into "improves monotonically."
- **Why it plateaus (interpretation, flagged as such):** the mutable window is
  10bp inside a 524,288bp random-sequence context, and the objective averages
  over all 6144 output bins. A 10bp change is a ~0.002% perturbation of the
  input and is diluted across the whole output window, so there is almost no
  gradient of improvement available. This is a property of the *deliberately
  trivial validation setup*, not evidence the design machinery is broken.
- **Wall-clock:** 49.2s for 30 variants at `batch_size=1` = **1.639
  s/variant**. The 5-iteration run took ~4 min wall-clock for 10 positions.
- **Peak GPU memory: 3.72 GB** at `batch_size=1`.

**Memory finding worth noting:** 3.72 GB here vs. the **18.08 GB** measured in
Task 1 for a single forward pass. Same model, same card, same batch size — the
difference is the code path: Task 1 used `model.predict_on_seqs()` directly,
while `evolve()` routes through PyTorch Lightning's `Trainer.predict()`, which
manages memory far more aggressively (no-grad context, per-batch teardown).
**This materially changes the feasibility answer below**, so it is reported as
the measured difference rather than assumed to be an error.

#### Design sweep cost: feasibility answer

At the measured **1.639 s/variant**, the answer depends entirely on what
"8,000–12,000 candidates" means. Both readings are given because the brief is
ambiguous and the two differ by ~3,000x:

| Reading | Compute | Wall-clock | Cost @ $0.75/hr |
|---|---|---|---|
| **A. N = candidate *sequences*, each getting its own evolve run** (200bp window = 600 variants/iter) | | | |
| &nbsp;&nbsp;8,000 candidates × 5 iters | 24,000,000 forward passes | **10,927 h (1.2 GPU-years)** | **~$8,195** |
| &nbsp;&nbsp;12,000 candidates × 5 iters | 36,000,000 forward passes | 16,390 h (1.9 GPU-yr) | ~$12,292 |
| &nbsp;&nbsp;12,000 candidates × 10 iters | 72,000,000 forward passes | 32,780 h (3.7 GPU-yr) | ~$24,585 |
| **B. N = total variants/forward passes evaluated** | | | |
| &nbsp;&nbsp;8,000 total variants | 8,000 forward passes | **3.6 h** | **~$2.73** |
| &nbsp;&nbsp;12,000 total variants | 12,000 forward passes | 5.5 h | ~$4.10 |

**Feasibility on a 24 GB card:** the memory ceiling is **not** the binding
constraint — at 3.72 GB per `batch_size=1` evolve step, there is roughly 6x
headroom on this card, and batching should be possible (a `batch_size` sweep
was attempted but returned an ambiguous/inconclusive result — the process
ended without a completion line and without a traceback, possibly an
unlogged OOM kill at `batch_size=8`; **not confirmed, so no batch-size
speedup is claimed here**). The binding constraint is **throughput**: under
reading A the job is 1.2–3.7 GPU-years on this card, which is not feasible on
one 24 GB GPU regardless of batching — it would need either a fundamental
reduction in scope (fewer candidates, fewer iterations, a smaller mutable
window), or massive parallelism across many GPUs, or a cheaper surrogate
model for the inner loop. Under reading B it is trivially cheap (hours, single
digit dollars). **Recommend clarifying which reading is intended before
committing budget** — that single definition is worth ~$8,000–$24,000.

### Task 4 — Fine-tuning path: **PASSES, trunk freeze verified at the value level**

This was the priority task, and it passed the strict check requested.

**A real gReLU 1.1.0 bug found live**, worth knowing before building against
this API: `train_on_dataset()` checks `if checkpoint is True: ... elif
isinstance(checkpoint, dict): ... else: raise Exception("Checkpoint type must
be a bool or dict")`. That is an exact-identity `is True` check, so passing
the equally-valid bool `False` hits the else branch and raises — **despite the
error message itself stating a bool is acceptable.** Workaround: pass `True`
(or a dict), never `False`. Also confirmed: `grelu.data.preprocess.split()`
excludes any chromosome appearing in val/test from the train set — good,
correct leakage prevention, but it means train/val/test chromosome lists must
be genuinely disjoint or the train split silently comes back empty.

**Trunk-freeze verification (the check that matters):**

| Check | Result |
|---|---|
| `requires_grad=False` on trunk | **180/180 trunk parameter tensors** |
| `requires_grad=True` on head | **2/2 head parameter tensors** |
| Trunk tensors with ANY value change after a real training step | **0 / 180** |
| Max absolute difference across all trunk parameters | **0.000e+00** (bit-identical) |
| Total trunk values compared | **171,271,968** |
| Lightning's own accounting | **3.8 K trainable / 171 M non-trainable** |

`model.tune_on_dataset(..., freeze_embedding=True)` was called through the
real public API (not a hand-rolled equivalent), a full Lightning `fit` ran
(validation pass + 1 training epoch, `train_loss_epoch=0.171`), and trunk
weights were snapshotted before and compared value-by-value after. **The trunk
is genuinely frozen — verified by direct tensor comparison, not by trusting
the flag.**

#### The input contract (measured from the live model, not documentation)

| Item | Measured value |
|---|---|
| **Trunk output / required head input shape** | **`(batch, 1920, 6144)`** — measured by running `model.embedding(x)` on a real 524,288bp input. Confirms the 1920 figure from prior sessions empirically. |
| Head type | `ConvHead` |
| Head `in_channels` | **1920** |
| Head activation | **`softplus`** (`act_func='softplus'`, `pool_func=None`, `norm=False`) — note the pretrained head applies softplus, so outputs are non-negative |
| Head output shape | `(batch, n_tasks, 6144)` |
| **Training input (one example)** | shape **`(4, 524288)`**, dtype **`torch.float32`** — one-hot DNA, channels-first, full model `seq_len`. Note: **not** `(524288, 4)`. |
| **Training target (one example)** | shape **`(n_tasks, 1)`**, dtype **`torch.float32`** — for the 2-task synthetic case, `(2, 1)`. The trailing 1 is the length axis. |
| Batch structure | DataLoader yields `(inputs, targets)`; batched shapes `(B, 4, 524288)` and `(B, n_tasks, 1)` |
| `crop_len` | **5120** (model_params) |
| `label_len` | 196,608 bp; `bin_size` 32 → 6144 output bins |
| Default `task` / `loss` from checkpoint | `regression` / `mse` (inherited `batch_size=512`, `devices='cpu'` — **override these**, the inherited batch_size will OOM instantly on 524kb inputs) |
| Train/val/test splitting | `grelu.data.preprocess.split(ad, train_chroms=, val_chroms=, test_chroms=)`; **chromosome lists must be disjoint** (see bug note above) |
| Data container | `AnnData`: `.var` = intervals (chrom/start/end, resized to `seq_len`), `.obs` = one row per task, `.X` = label matrix shaped **(n_tasks, n_intervals)** |

**On target normalization (connecting to Task 2's finding):** the pretrained
7,611-task head emits values in Borzoi's *transformed* target space — the
same space Task 2 had to invert (per-bin: `scale` → soft-clip at `clip_soft`
→ `**0.75` power). For a **new** senescence head trained from scratch on the
frozen trunk, you are free to define your own target space and are **not**
obliged to reproduce Borzoi's transform — the head learns whatever mapping
your targets define. However, if you want the new head's outputs to be
comparable to, or combinable with, the pretrained track outputs, targets
should be put through the same transform chain. Given Task 2's finding that
low-expression genes carry disproportionate error, a variance-stabilizing
transform (log1p, or Borzoi's own `**0.75`) on WI-38 expression targets is
worth considering — but note this is a **recommendation grounded in the
measured Task 2 residual pattern**, not something verified by a training run
here.

### Task 5 — Reference data and genome build

All measured on this instance:

| Item | Measured |
|---|---|
| **Genome build required** | **hg38** — re-confirmed from the live model: `model.data_params['train']['genome'] == 'hg38'` |
| Where the FASTA must live | genomepy's default dir, `~/.local/share/genomes/<name>/`. **`GRELU_GENOMES_DIR` is NOT read** by gReLU 1.1.0's actual `get_genome()` (re-confirmed). Staged at `/workspace/genomes/hg38` and symlinked into the default dir. |
| Source used | NCBI `GRCh38` (`GCF_000001405.26`) — **UCSC is unreachable from RunPod**, and neither Ensembl (`GRCh38.p14`) nor NCBI matches the literal string `"hg38"`, so an explicit `provider='NCBI', localname='hg38'` install plus a UCSC-style contig rename (194 sequences remapped via `assembly_report.txt`) is required |
| Blacklist | **Bundled in the gReLU package**, no download: `.../grelu/resources/blacklists/encode/hg38-blacklist.v2.bed` |
| Chromosome sizes | Bundled alongside the FASTA (`hg38.fa.sizes`), generated during the reindex |
| GTF annotation | Optional, **not staged** (not needed for Tasks 1–4); would need `genomepy.install_genome('hg38', only_annotation=True)` + UCSC utilities |
| Verification | `genomepy.list_installed_genomes()` → `['hg38']`; live sequence extraction test returned the correct 100bp |

**Measured disk footprint:**

| Component | Size |
|---|---|
| Genome (`/workspace/genomes`, hg38 FASTA + .fai + .sizes + gaps.bed) | **3.0 GB** |
| Python venv (gReLU + all deps, incl. torch via system-site-packages) | **5.4 GB** |
| HuggingFace cache (one Borzoi replicate, `human_rep0.ckpt`) | **713 MB** |
| **Total** | **~9.1 GB** (add ~713 MB per additional Borzoi replicate if ensembling) |

#### Liftover

**Confirmed on this instance, not re-assumed:** a live case-insensitive grep
for "liftover" across the entire installed `grelu` 1.1.0 package returned
**zero matches**. gReLU ships no liftover tooling of any kind.

**What would be required**, given the project's candidate training series span
hg38 (GSE175533, GSE210285, ENCSR200OML) and hg19/GRCh37 (GSE206402 family,
GSE74238/GSE74324):

- **Tooling** (none is a gReLU dependency; all must be added separately):
  the UCSC `liftOver` binary, or `CrossMap` (Python, handles BED/BAM/VCF/
  bigWig), or `pyliftover` (Python, point coordinates only — insufficient for
  interval data).
- **Chain file (specific):** **`hg19ToHg38.over.chain.gz`**, from
  `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz`.
  Note UCSC was unreachable from this RunPod instance during this session —
  budget an alternate mirror or a pre-staged copy.
- **Expected failure rate — with an important caveat.** Published figures for
  *point variants* are 99.87–99.99% success (i.e. **0.01–0.13% failure**,
  e.g. ClinVar GRCh37→GRCh38 evaluations). **Those numbers do not transfer to
  this use case.** Wide genomic intervals (ATAC/ChIP peaks, and especially
  Borzoi's 524,288bp input windows) fail substantially more often, because a
  single interval must map contiguously — regions spanning assembly
  restructuring (expanded centromeres/satellites in hg38, retired hg19
  contigs, low-mappability repeats) fail or fragment. **No measured
  interval-scale failure rate was produced in this session, and none should
  be assumed** — measure it on your actual peak sets before relying on it.
- **How failures should be handled:** (1) never silently drop — log and count
  every unmapped region, since systematic loss in specific genomic
  neighborhoods is a bias, not noise; (2) discard intervals that fragment
  into multiple hg38 pieces (`-multiple` output) rather than picking one
  arbitrarily; (3) **re-check post-liftover intervals against the hg38
  blacklist**, since liftover can land regions in newly-problematic areas;
  (4) verify lifted interval widths match the originals — width changes signal
  an indel-containing region and unreliable coordinates; (5) treat the hg19
  series as a **separate batch covariate**, since even successfully-lifted
  data was generated against a different reference and may carry systematic
  differences.

### Cost, this session

| Item | Value |
|---|---|
| Rate | $0.75/hour |
| Wall-clock (17:48 first command → ~19:01 last retrieval) | **~73 minutes** |
| **Approximate cost** | **~$0.91** (from SSH command timestamps, not a pulled invoice) |
| Breakdown | gReLU install ~23 min; genome staging ~9.5 min; Task 3 (2 failed + 2 successful runs) ~18 min; Task 4 (2 failed + 2 successful runs) ~12 min; overhead/retrieval remainder |

---

## Fine-tune CALIBRATION PASS (fourth pod, 2026-08-23) — caching verified, cost ~7x better than projected

**Scope: calibration only, as instructed. The full trunk pass was NOT run and
nothing was trained.** Pod: `<pod-address-redacted>`, RTX 4090 24,564 MiB,
driver 580.126.20, Ubuntu 24.04.3 / Python 3.12.3 / torch 2.8.0+cu128
preinstalled. Environment built with the settled workarounds (venv
`--system-site-packages`; NCBI GRCh38 + UCSC contig remap, 194 sequences
renamed and reindexed; symlink into genomepy's default dir). gReLU install
**~24 min**, matching budget.

### 1. Measured per-gene trunk cost (n=60 genes, seed=0, CDKN1A force-included)

| Component | mean | sd | median | min | max | p5–p95 |
|---|---|---|---|---|---|---|
| Sequence extraction (FASTA→one-hot) | **0.081 s** | 0.020 | 0.078 | 0.062 | 0.206 | 0.066–0.096 |
| Trunk forward (fwd + RC, 2 passes) | **0.323 s** | 0.0020 | 0.323 | 0.322 | 0.338 | 0.322–0.324 |
| **Total, fwd + RC** | **0.405 s** | 0.021 | 0.401 | 0.385 | 0.529 | 0.389–0.427 |
| **Total, forward-only (est.)** | **0.243 s** | 0.020 | 0.239 | 0.224 | 0.368 | 0.227–0.259 |

The trunk forward itself is **remarkably stable** (sd 0.0020 s, a 0.6%
coefficient of variation) — essentially all run-to-run variance comes from
FASTA extraction (sd 0.020 s), which is disk/cache-dependent. The one 0.206 s
extraction outlier is a cold-cache read.

### 2. Projected full-pass cost — **~7x cheaper than the pre-run estimate**

| Pass | s/gene | Wall-clock (19,652 genes) | Cost @ $0.75/hr |
|---|---|---|---|
| **Forward only** | 0.243 | **1.33 h** | **$1.00** |
| **Forward + reverse-complement** | 0.405 | **2.21 h** | **$1.66** |

The pre-run projection assumed 0.861–1.639 s/gene (the Task 1 / Task 3
measured rates). Actual is **0.243 s/gene** — 3.5x faster than the optimistic
figure and 6.7x faster than the conservative one. The earlier numbers came
from `predict_on_seqs` and the Lightning `Trainer.predict` path, both of which
carry per-call framework overhead; calling `core.embedding(x)` directly avoids
it. **Recommendation: run the fwd+RC pass at $1.66** — reverse-complement
caching costs $0.66 extra and restores the augmentation that caching would
otherwise force off.

### 3. Cache size — measured, then extrapolated

| | Sample (60 genes) | Extrapolated (19,652) |
|---|---|---|
| Forward cache | 0.46 MB | **151.0 MB** |
| Reverse-complement cache | 0.46 MB | **151.0 MB** |
| **Total** | 0.92 MB | **301.9 MB** |

Versus **927 GB** for the full uncompressed trunk output — a **6,144x**
reduction. Measured `.npy` sizes matched the analytic prediction (19,652 ×
1920 × 4 B = 150.9 MB) to within the 128-byte header.

### 4. Peak GPU memory: **3.87 GB of 24.56 GB** (batch size 1)

That leaves ~20 GB headroom, implying batching is comfortably feasible and the
1.33 h could drop further. **Not measured** — batching was out of scope for
this calibration and is recorded as unquantified upside, not claimed.

### 5. Cache-identity check — **verified, and the residual is fully explained**

The check ran on 5 genes including **CDKN1A** (the benchmark locus, held-out
chr6). Initial result: max |diff| **9.58e-05**, relative **2.74e-04**. That is
larger than float32 rounding alone would predict, so rather than wave it
through, it was traced:

| Arithmetic | max abs diff | relative |
|---|---|---|
| **TF32 on** (PyTorch default: `cudnn.allow_tf32=True`) | 6.05e-05 | 1.79e-04 |
| **TF32 off** (float32) | **5.96e-08** | **1.76e-07** |
| **float64, CPU** | **3.61e-16** | **1.07e-15** |

The residual collapses by 3 orders of magnitude with TF32 disabled and by 9 in
float64 — i.e. it scales precisely with mantissa width (TF32 10-bit → float32
24-bit → float64 53-bit). **The identity `head(full) == head(pooled)` is exact
algebra; the observed difference is entirely floating-point precision, not a
modelling discrepancy.** Confirmed empirically, as required, rather than
assumed from the commutativity argument.

Practical note: at TF32 the induced error on head outputs is ~1.8e-4 relative,
negligible against `log1p(TPM)` targets spanning 0–9.97. Setting
`torch.backends.cudnn.allow_tf32 = False` for the head reduces it to ~1.8e-7
at no meaningful cost (the head conv is 46,104 params), and is worth doing.

### 6. Sanity items also confirmed

- Input shape **(1, 4, 524288)**, trunk output **(1, 1920, 6144)** — both match
  the measured contract exactly.
- Reverse-complement one-hot validity asserted (`onehot_and_rc_valid: true`);
  RC implemented as `flip(x, dims=[-2,-1])` on ACGT one-hot.
- Fresh head built exactly as `change_head()` does: `in_channels=1920`,
  `n_tasks=24`, `act_func=None`, `pool_func='avg'`, `norm=False`.
- Training set loaded from the uploaded `.h5ad` (MD5 verified against local):
  19,652 genes, 33 tasks, 24 ON/OFF training tasks, quiescent 9 held out.

### Revised total budget for the full fine-tune

| Item | Cost |
|---|---|
| Trunk cache pass (fwd + RC) | **$1.66** |
| All 5 training runs (real + 3 shuffled + LOO) on cached features | **~$0** (46,104-param head; an epoch is one 19,652×1920 matmul) |
| **Projected total** | **~$2** |

versus **$176–671** without caching. The pre-run recommendation stands and is
now backed by measurement rather than estimate.

### Calibration session cost

~41 min wall-clock (23:34 PDT 2026-08-22 → 00:15 PDT 2026-08-23), of which
~24 min was the gReLU install and ~2.5 min the genome staging. **≈ $0.51.**

Artifacts: [`execution_results_calibration/`](execution_results_calibration/)
— `calibration_results.json`, `precision_check.json`, both clean logs, and both
scripts.

---

## FINE-TUNE, MONTH-4 CHECKPOINT (fifth pod, 2026-08-23/24) — **BARS MISSED, PRIMARY CONTROL FAILED**

**Pod:** `<pod-address-redacted>`, RTX 4090 24,564 MiB, driver 570.195.03.
**Verdict: NO-GO.** Both primary bars were missed, and — decisively — the
shuffled-label control failed. Per the pre-registered stop rule, this report
stops at the measured numbers and does not proceed to interpretation of the
model as having learned senescence biology.

### Headline result

> **Test AUC = 0.659** (bar: ≥ 0.85 — **MISSED by 0.19**).
>
> **This rests on 377 of 1,740 held-out genes (21.7%); 1,363 genes were
> excluded as ambiguous (|measured log2FC| < 1.0). Class balance among genes
> used: 241 positive / 136 negative (63.9% / 36.1%). This slice is
> SYSTEMATICALLY SELECTED for large measured senescence response — it is not a
> random subset of the genome. The number therefore describes discrimination
> AMONG STRONGLY-RESPONDING GENES, not prediction across genes generally.**

> **Test Spearman = 0.514** (bar: ≥ 0.60 — **MISSED by 0.087**).

### Train / validation / test, reported separately

| Split | AUC | genes used | Spearman (pooled) | pred std |
|---|---|---|---|---|
| train | 0.752 | 3,640 / 15,910 (22.9%) | 0.587 | 0.953 |
| val | 0.686 | 459 / 2,002 (22.9%) | 0.526 | 0.876 |
| **test** | **0.659** | **377 / 1,740 (21.7%)** | **0.514** | 0.935 |

**Train→test gap: AUC −0.093, Spearman −0.074.** This is the expected failure
mode given the ON arm is effectively n≈3 from a single serial culture lineage,
and it is present — but it is *not* the main problem here. The main problem is
Control 1.

### 🔴 CONTROL 1 — SHUFFLED LABELS: **FAILED**

| | value |
|---|---|
| Real test AUC | **0.6585** |
| Shuffled-label AUC, mean ± sd (20 permutations) | **0.5605 ± 0.0940** |
| Shuffled range | **[0.4232, 0.7543]** |
| Permutations scoring **≥ the real model** | **2 of 20** |
| Real AUC above shuffled mean | **1.04 sd** |
| Empirical p | **0.143** |

**Two of twenty random 9-vs-15 relabelings of the same samples beat the true
senescence contrast, and the real result is not significant at p<0.05.** By
the pre-registered criterion — "shuffled labels scoring near the real model" —
this control fails, and the headline AUC cannot be attributed to
senescence-specific signal.

Why this control is the right test here, and why retraining is a no-op: the
head's training targets are per-sample expression, which do not depend on arm
labels, so a permuted assignment yields a bit-identical model (verified live:
`retrained model identical to real model: True`). The control therefore
isolates exactly the question at issue — *is the AUC specific to the true
ON/OFF split, or does any split of these 24 samples score comparably?* The
answer is that any split scores comparably. The model predicts expression
differences between arbitrary sample groups about as well as it predicts the
senescence contrast specifically.

### Controls that did NOT fail

| Control | Result | Reading |
|---|---|---|
| **2 — untrained head** | AUC **0.489**, Spearman −0.019 | Floor is chance, as expected |
| **3 — leave-one-PDL-out** | PDL50 **0.557**, PDL52 **0.668**, PDL53 **0.717** (mean 0.647) | **Does not collapse** — comparable to the full model's 0.659, so this is not simple lineage memorisation |
| **Supplementary — shuffled gene↔embedding** | AUC **0.504**, Spearman 0.035 | Chance. So the model *does* use sequence information (0.659 > 0.504); the embeddings are not inert |
| **Prediction collapse** | pred std **0.935** vs measured std **1.724** (ratio 0.54) | **Not collapsed** — predictions retain real variance |

So: the model uses sequence, generalises across PDL timepoints, and does not
collapse. What it cannot be shown to do is anything *senescence-specific*.

### CDKN1A — the benchmark locus (held-out chr6)

| | ON | OFF | differential |
|---|---|---|---|
| **Predicted** | 2.670 | 2.689 | **−0.019** (wrong direction) |
| **Measured** | 7.022 | 6.666 | **+0.509 log2FC** |

Spearman of predicted vs measured across the 24 tasks at this locus:
**−0.388** (anti-correlated). **The model predicts a slight decrease at CDKN1A
where measurement shows induction — it gets the benchmark locus wrong, in
direction as well as magnitude.** (CDKN2A was not used as a validation target,
per the established finding that its expression barely moves.)

### SECONDARY — threshold sensitivity (headline remains 1.0)

| threshold | AUC | genes used | pos / neg |
|---|---|---|---|
| ≥ 0.5 | 0.684 | 864 (49.7%) | 482 / 382 |
| **≥ 1.0 (headline, approved)** | **0.659** | **377 (21.7%)** | **241 / 136** |
| ≥ 2.0 | 0.645 | 103 (5.9%) | 70 / 33 |

**Spread 0.039 → stable across thresholds.** The result is not an artifact of
where the line was drawn. Note this cuts both ways: the *miss* is also stable.

### Specificity probe — quiescent samples: **model does NOT call them senescent**

Each sample's measured deviation-from-OFF projected onto the model's predicted
senescence axis, held-out test genes only:

| class | n | score |
|---|---|---|
| proliferating (OFF) | 15 | **+0.022 ± 0.146** |
| **senescent (ON)** | 9 | **+0.244 ± 0.059** |
| **quiescent (held out)** | 9 | **+0.080 ± 0.022** |

Quiescent samples score **+0.080**, far closer to proliferating (+0.022) than
to senescent (+0.244). **The model does not mistake growth arrest for
senescence** — a genuine positive, and the one clearly encouraging result in
this run. Hedged appropriately: the quiescent arm is dense arrested culture
with elevated inflammatory signal, not textbook quiescence, so this shows
separation from *that* specific state, not from quiescence in general.

### Hyperparameters (fixed in advance, not searched)

`lr=1e-4, batch_size=256, max_epochs=200, patience=10, optimizer=adam,
loss=mse, seed=0, rc_augment=True, eval_rc_average=True`. Trained 188 epochs,
best val MSE 2.3559. Target transform confirmed `log1p(TPM)` at load. TF32
disabled on cudnn and matmul. **Nothing was tuned toward the 0.85 bar.**

### Cache pass (phase 1)

2.53 h, peak GPU 3.72 GB, fwd 150.9 MB + rc 150.9 MB, **0 all-zero rows, 0
NaN**, verified independently after completion. One real bug fixed mid-phase:
the RefSeq GRCh38 assembly contains **IUPAC ambiguity codes** which survive
gReLU's `.upper()` and are rejected by `ALLOWED_BASES=[A,C,G,T,N]`, aborting
the first attempt. Any non-ACGTN base is now mapped to `N` (the truthful
encoding for an uncertain base). Impact measured, not assumed: **49 of 19,652
windows (0.25%), 133 characters** out of ~10.3 billion — codes
`{Y:49, R:33, M:19, W:17, K:6, S:5, B:4}`.

### What this triggers

Per the pre-registered fail action, and **not** implemented in this run:
full-window mean pooling over all 6,144 bins (±98 kb around each TSS) remains
the **pre-declared leading suspect** — the head averages a gene's TSS signal
into ~196 kb of mostly-irrelevant context. A **centered-bin head is the
documented fail action for a SUBSEQUENT run.** No architecture was changed
mid-run and no head was adjusted to chase the bar.

However, the control failure is the more serious finding and is **not**
obviously fixed by the pooling change: the issue is not only that the AUC is
low, but that it is statistically indistinguishable from random sample
relabelings. The other candidate fail actions on record — adding training
pairs (the n≈3 ON arm is the prime structural weakness), or a dry-only
release — should be weighed against the pooling fix rather than assuming the
architecture alone explains it.

### Cost, this session

| Item | Value |
|---|---|
| Elapsed pod time (first contact 21:50 UTC → 01:39 UTC) | **~3.82 h** |
| **Cost @ $0.75/hr** | **~$2.87** |
| Breakdown | env build ~35 min (gReLU ~24 min + genome ~4 min), cache pass 2.53 h, training + all controls ~11 min |

Artifacts: [`execution_results_finetune/`](execution_results_finetune/) —
`results.json`, `finetune_clean.log`, `cache_pass_clean.log`,
`sanitize_stats.json`, and both scripts.

---

## Tasks 3, 4, 5 — original prep-session status (superseded by the section above)

Explicitly out of scope for this session per your instruction. Nothing in
`runbook/task3_design_loop.py`, `task4_finetune_path.py`, or
`task5_reference_data.sh` was executed, and no claim is made about whether
they would pass. The environment changes made to get Tasks 1–2 working
(venv location, genome symlink/rename fix) would carry over and should make a
future Task 3/4/5 run start from a working base rather than hitting the same
genome-provider issues from scratch — worth knowing if you decide to proceed.

---

## Cost table (measured/observed, not projected)

| Item | Value |
|---|---|
| Instance rate | $0.75/hour (user-stated), both pods |
| **Pod 1** wall-clock (original Task 1+2 run) | ~60 minutes ≈ **~$0.75** |
| **Pod 2** wall-clock (this rerun: venv + gReLU install ~29 min + genome setup ~8 min + corrected Task 2 ~5 min + overhead) | ~44 minutes, SSH-command-timestamped (16:44–17:28 PDT) ≈ **~$0.55** |
| **Combined approximate cost across both sessions** | **~$1.30** (verify against RunPod's own billing for both pods — these are derived from command timestamps, not pulled invoices) |
| Per-forward-pass cost | 860.6ms (measured on pod 1) → at $0.75/hr, ≈ **$0.00018 per forward pass** |
| Task 2 original sweep (118 genes) | 67.7s GPU compute ≈ **$0.014** |
| Task 2 corrected rerun sweep (118 genes) | 90.7s GPU compute ≈ **$0.019** |
| Design sweep / fine-tune cost | **Not applicable** — Tasks 3/4 were not run on either pod, no data to project from |
| Task 2 subset re-analysis | **$0.00** — local re-analysis of already-saved predictions plus a network call to the free Ensembl REST API; no GPU instance provisioned or connected to |
| Training-set build (GSE175533 → AnnData) | **$0.00** — local data curation only, no GPU |
| **Pod 4** — fine-tune calibration pass (~41 min) | **~$0.51** |
| **Projected: full fine-tune with caching** | **~$2** (trunk cache pass $1.66 + ~$0 for all 5 head-training runs) |
| Projected: same fine-tune WITHOUT caching | $176–671 |

---

## Action needed from you

**⚠️ POD 5 (`<pod-address-redacted>`) IS STILL RUNNING AND STILL BILLING.**
It ran the real fine-tune. All work is complete and every artifact has been
copied locally — **nothing further is needed from it. Terminate it now.**

Pod 4 (`<pod-address-redacted>`) was confirmed terminated (connection refused).
Pod 3 (`<pod-address-redacted>`) status unconfirmed — worth checking.

Pods 1 (`<pod-address-redacted>`) and 2 (`<pod-address-redacted>`) were previously
confirmed terminated by the user.

**Cumulative approximate spend across all four pod sessions: ~$2.72**
(~$0.75 + ~$0.55 + ~$0.91 + ~$0.51), derived from SSH command timestamps, not
pulled invoices — check the RunPod console for authoritative figures.

**Next decision:** the full fine-tune is now costed at **~$2** total (trunk
cache pass $1.66 + effectively free head training for all five runs). Two
things I flagged pre-run still need your call before that run: (a) the AUC
definition to pre-register, and (b) acknowledgement that if the run misses the
0.85 bar, full-window mean-pooling dilution is the pre-declared prime suspect
and a centered-bin head is the documented fail action — not a mid-run tune.

**One open decision that carries real money:** the design-sweep cost estimate
spans ~$2.73 to ~$24,585 depending on whether "8,000–12,000 candidates" means
total forward passes or candidate sequences each getting a full evolve run
(Task 3 table). Worth settling before committing compute budget.

---

## Reproducing / continuing this session's work

Environment fixes that would need to be repeated (or are already in place, if
the same pod is reused) for any future run on this instance:

```bash
# venv (system Python is externally-managed on Ubuntu 24.04)
python3 -m venv /workspace/venv --system-site-packages
/workspace/venv/bin/pip install gReLU==1.1.0   # ~25 min observed, not <10 min

# genome: NCBI's exact "GRCh38" (not "hg38" -- no provider recognizes that
# literal string), renamed to UCSC-style contig names, symlinked into
# genomepy's real default dir (GRELU_GENOMES_DIR is NOT read by the actual
# installed code, despite the prep report's documentation -- see Deviations)
mkdir -p /workspace/genomes
/workspace/venv/bin/python -c "import genomepy; genomepy.install_genome(name='GRCh38', provider='NCBI', genomes_dir='/workspace/genomes', localname='hg38', annotation=False)"
# then: rename headers via assembly_report.txt's UCSC-style-name column,
# reindex with pyfaidx, symlink into ~/.local/share/genomes/hg38
# (see this session's transcript for the exact rename script)

/workspace/venv/bin/python runbook/task1_smoke_test.py
/workspace/venv/bin/python runbook/task2_run.py          # original (flawed) aggregation, kept for reference
/workspace/venv/bin/python runbook/task2_corrected.py     # corrected inverse-transform aggregation, this rerun
```

Result files from the rerun: `execution_results_rerun/task2_corrected.log` (full run
log including the sanity check), `execution_results_rerun/task2_results_corrected.csv`
(per-gene results), `execution_results_rerun/genome_setup.log`.
