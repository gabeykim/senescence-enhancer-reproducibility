# Senescence enhancer reproducibility

Code, derived results, figures and provenance records for a study of whether
senescence-associated enhancer activity is reproducible enough to model, and —
where it is — how far a sequence-to-activity model can be pushed toward
designing one.

- **Manuscript:** _(link to be added on preprint posting)_
- **Large artifacts (Zenodo):** [10.5281/zenodo.23146107](https://doi.org/10.5281/zenodo.23146107) — see
  [What is not here](#what-is-not-here)
- **Licences:** [MIT](LICENSE) for everything under `scripts/` and `env/`;
  [CC-BY 4.0](LICENSE-DATA) for `results/`, `figures/`, `audit/`,
  `audit_fixes/` and `reports/`.

This repository keeps its working reports exactly as they were written and
documents its own corrections rather than quietly fixing them: five
discrepancies found between those reports and the stored data are listed in
[reports/ERRATA.md](reports/ERRATA.md). The manuscript carries the corrected
values throughout.

---

## Read this first: every H3K27ac arm is n = 2

**Every H3K27ac dataset used here has at most two replicates per arm.** The
primary training contrast (GM21 RAS vs EV) is n = 2 vs 2; so is the held-out
validation contrast (IMR-90 SEN vs PRO). Only the GM21 ATAC data reaches
n = 4 vs 4. This is not a choice — it is what is deposited.

Everything downstream inherits that ceiling. It is why the paper reports a
measured **cross-dataset reproducibility ceiling of +0.637** and scores the
model against that rather than against 1.0, and why no result here should be
read as a population estimate. [Figure 1](figures/figure1_data_landscape.png)
exists to make this visible before anything else.

---

## The four findings

**1. Replicative senescence has no reproducible H3K27ac programme.**
Two replicative datasets agree with each other no better than a replicative
dataset agrees with an oncogene-induced one: the cleanest comparison available —
same lab, same pipeline, same genome build, no technical confound — is
Spearman **+0.054** (n = 134,336), against a within-OIS ceiling of **+0.653**
(n = 88,530). The promoter-level biology reproduces everywhere (CDKN2A falls in
all seven dataset arms), so this is not measurement failure; it is the absence
of a shared distal programme, and it is what forced the project onto OIS.

**2. A trained head on Borzoi trunk embeddings predicts the OIS response
across studies.** Training on GM21 and validating on IMR-90 — a cross-lab,
cross-tissue holdout — gives Spearman **+0.5484** on 100,070 held-out regions,
which is **86.1% of the +0.637 reproducibility ceiling**, with AUC **0.9070** on
the 26,936 evaluable regions. All three controls collapse to chance: 0 of 200
shuffled-label permutations beat it (p = 0.00498, the resolution floor), an
untrained head gives −0.0903, and shuffling the region-to-embedding assignment
takes +0.548 down to +0.025.

**3. The motif-insertion design lever saturates well below what is needed.**
Moving a sequence from the median to the 95th percentile of predicted activity
requires Δ = **+0.861**; the best condition found reaches **+0.666**, or
**0.773×**, and the dose response is non-monotonic — it peaks at 12 inserted
copies and *declines* at 18 and 20, because a 200 bp cassette tiled end to end
with consensus stops looking like an enhancer. Combining motif families buys
nothing over using more of the best single family at matched load
(+0.064, p = 0.157), while *arrangement* is worth up to +0.294, so whatever
headroom remains is in grammar and spacing rather than in motif count.

**4. The AP-1 anomaly is context dependence, not a claim about AP-1.**
Inserting the AP-1 consensus into inert background *lowers* predicted senescence
H3K27ac (Δ = −0.464), which would imply that adding AP-1 sites makes a designed
element worse — but ablating *native* AP-1 sites reverses the sign depending on
where they sit: −0.1055 in high-response regions versus +0.0658 in low-response
ones (interaction p = 0.0091), with an NF-κB positive control behaving as
expected. Genome-wide, strong AP-1 sites are **2.3× enriched** in the top
H3K27ac decile (p = 5.9 × 10⁻⁹⁴), so the insertion probe was simply only ever run
in the one context class where AP-1 is negative.

---

## What is here

```
scripts/      58 scripts, numbered in execution order; each carries a
              CONSUMES / PRODUCES header
env/          lockfiles + ENVIRONMENT.md (seven quirks that cost real time)
results/      derived tables, summary statistics and run logs, by phase
figures/      the twelve publication figures (PNG 300 dpi + vector PDF)
              and captions.md
audit/        the number audit: 72 claims checked against source files
audit_fixes/  the eight closed audit items, with their derivations
reports/      the working reports, unaltered except for redacted host
              addresses, plus ERRATA.md
data/         GSE254358 provenance evidence only (see below)
```

**The reports are a working record, superseded by the manuscript.** They are
kept because the record of what was concluded when is worth having, not because
they are correct. Five known discrepancies between a report and the stored data
are listed in [reports/ERRATA.md](reports/ERRATA.md), each with the correct
value and where the corrected version appears. Where a report and a file in
`results/` disagree, **the file is authoritative**.

### `data/GSE254358/` is evidence, not data

An accession cited across two research rounds turned out not to exist. The
directory holds the three files a fetch script wrote when NCBI returned 404:
two `.gz` files that are not gzip at all but byte-identical Apache error pages,
and a text file scraped from the same page. They are preserved byte-for-byte
with their original timestamps and checksums because the manuscript makes a
provenance claim about them. See
[data/GSE254358/README.md](data/GSE254358/README.md). Re-verified live on
2026-10-04: `esearch` returns 0, the FTP path 404s while both the parent range
directory and the adjacent accession GSE254357 return 200.

---

## What is not here

| Excluded | Why | Where to get it |
|---|---|---|
| `region_fwd.npy`, `region_rc.npy` (**818 MB each**) | Borzoi trunk embedding caches. Far over GitHub's 100 MB hard limit. | [Zenodo](https://doi.org/10.5281/zenodo.23146107), or regenerate with `scripts/81_ois_cache.py` (~1.16 h on one 24 GB GPU) |
| `cache_fwd.npy`, `cache_rc.npy` (97 MB each) | Expression-phase trunk caches | [Zenodo](https://doi.org/10.5281/zenodo.23146107) |
| `signal_matrix*.csv`, `region_responses*.csv`, `anchors*.csv` (28–66 MB) | Large intermediate matrices, fully re-derivable | [Zenodo](https://doi.org/10.5281/zenodo.23146107), or rerun steps 51/61/70 |
| Raw GEO/ENCODE data (9.4 GB) | Public and re-fetchable; redistribution adds nothing | `scripts/01, 10, 20, 22, 40, 50, 60` fetch it |
| `*.h5ad` training sets (14 MB) | Re-derivable from the fetch + build scripts | [Zenodo](https://doi.org/10.5281/zenodo.23146107), or rerun steps 43/71 |

The Zenodo record holds **24 files, 2.35 GB**, with a SHA-256 manifest and its
own README describing what each file is, which figure it supports, and what it
cannot be used for. `10.5281/zenodo.23146107` is the version DOI and always
resolves to that exact deposit; the concept DOI `10.5281/zenodo.23146106`
resolves to the newest version instead.

**Git LFS is deliberately not used.** Its quota and bandwidth limits break
anonymous cloning, which would defeat the point of a reproducibility
repository. Large artifacts go to Zenodo instead.

The trained head itself (`results/ois_model/head_ois.pt`, 10 KB) **is** here —
it is a single 1920→1 convolution and needs no special handling.

---

## Reproducing each result

### What needs a GPU and what does not

| Steps | Hardware | Wall-clock |
|---|---|---|
| `01`–`71` — fetch, QC, region building, trainset assembly | **CPU only** | hours, mostly download |
| `80`–`84` — env, trunk cache, head training, probe | **One 24 GB GPU** | ~3.5 h (cache 1.16 h, probe 1.7 h) |
| `91`–`94` — design gates | **One 24 GB GPU** | ~0.85 h |
| `93`, `95`–`99b` — scan, audit, figures | **CPU only** | minutes |

Head *training* is effectively free once the cache exists — an epoch is one
`111671 × 1920` matmul. The GPU cost is almost entirely the trunk forward
passes.

### Finding 1 — the reproducibility matrix (no GPU)

```bash
bash scripts/50_fetch_enhancer_signal.sh      # bigWigs + peaks
bash scripts/60_fetch_replicative_ceiling.sh
python scripts/51_build_enhancer_regions.py   # 237,824 anchors, hg19 -> hg38
python scripts/61_replicative_ceiling.py      # hg19 -> hg18, quantify each arm
python scripts/62_replicative_ceiling_analysis.py
python scripts/63_shared_denominator_null.py
```
Check against `results/replicative_ceiling/replicative_correlations.csv` and
`results/enhancer_cross_mechanism/cross_mechanism_correlations.csv`.

### Finding 2 — the OIS model (GPU for 81–83)

```bash
python scripts/70_build_ois_regions.py
python scripts/71_ois_trainset.py             # -> 111,671-region AnnData
bash   scripts/80_setup_pod.sh                # GPU; see env/ENVIRONMENT.md
python scripts/81_ois_cache.py                # GPU, ~1.16 h -> the 818 MB caches
python scripts/84_verify_cache_artifacts.py   # shape / finite / fwd != rc
python scripts/82_ois_train.py                # GPU; head + 200 permutations
python scripts/83_ois_probe.py                # GPU, ~1.7 h
```
Check against `results/ois_model/results_train.json`. To skip the GPU entirely,
take the caches from Zenodo and run `82`/`83` only — or skip even that and use
the stored `pred_rest.npy` + `head_ois.pt`.

### Findings 3 and 4 — the design gates (GPU for 91, 92, 94)

```bash
python scripts/93_gate2_scan.py               # CPU, ~30 s
python scripts/91_gate1_additivity.py         # GPU
python scripts/92_gate1b_ceiling.py           # GPU
python scripts/94_gate2_ap1.py                # GPU
python scripts/95_gate2_finish.py             # CPU
python scripts/96_gate2_h3.py                 # CPU
python scripts/97_recompute_additivity_contrast.py   # CPU
```
Check against `results/design_gates/`. Findings 3 and 4 depend only on
`head_ois.pt` and the Borzoi trunk, both small or publicly downloadable.

### The figures (no GPU, no model)

```bash
pip install -r env/local-analysis.lock.txt
python scripts/99_figures_main.py             # Figures 1-7
python scripts/99b_figures_supp.py            # Figures S1-S5
```
Two self-checks run inside these and will fail loudly if the stored data has
drifted: the ROC AUC recomputed from `pred_rest.npy` must equal the stored
0.9070, and S1's same-direction count must equal the stored 12,429.

---

## Figure → script → source data

| Figure | Script | Source data |
|---|---|---|
| 1 Data landscape | `99_figures_main.py::fig1` | `figures/dataset_table.csv`, built from `results/replicative_ceiling/step2_stats.json`, `scripts/60_fetch_replicative_ceiling.sh`, `results/ois_trainset/trainset_obs.csv`, `results/ois_trainset/STEP1_LOG.txt` |
| 2 Reproducibility matrix | `99_figures_main.py::fig2` | `results/replicative_ceiling/replicative_correlations.csv`; `results/enhancer_cross_mechanism/cross_mechanism_correlations.csv` |
| 3 Threshold sensitivity | `99_figures_main.py::fig3` | same two files, all four thresholds |
| 4 CDKN2A promoter vs distal | `99_figures_main.py::fig4` | `results/replicative_ceiling/locus_sanity_checks.csv` |
| 5 OIS model performance | `99_figures_main.py::fig5` | `results/ois_model/results_train.json`, `pred_rest.npy`, `idx_rest.npy`, `measured_IMR90_SEN_vs_PRO_rest.npy` |
| 6 Design ceiling | `99_figures_main.py::fig6` | `results/design_gates/gate1b_dose_curve.csv`, `gate1_spacing.csv`, `gate1_results.json`; `audit_fixes/additivity_contrast_per_seed.csv` |
| 7 AP-1 context dependence | `99_figures_main.py::fig7` | `results/ois_model/results_probe.json`; `results/design_gates/gate2_ap1_results.json`, `gate2_h3_frequency_by_decile.csv` |
| S1 ATAC vs H3K27ac | `99b_figures_supp.py::s1` | `results/enhancer_cross_mechanism/step2_stats.json`, `gm21_atac_vs_k27ac.csv` |
| S2 Batch structure | `99b_figures_supp.py::s2` | `results/ois_trainset/batch_pca_stats.csv`, `STEP2_LOG.txt` |
| S3 Shared-control artifact | `99b_figures_supp.py::s3` | `results/replicative_ceiling/shared_denominator_null.txt`, `replicative_correlations.csv` |
| S4 Exhaustive 3v3 null | `99b_figures_supp.py::s4` | `results/expression_phase/results_perm200.json` |
| S5 Motif dose response | `99b_figures_supp.py::s5` | `results/ois_model/motif_insertion.csv`, `results_probe.json` |

**Every source file in this table is in this repository**, so all twelve figures
regenerate from a clone with no Zenodo download and no GPU. Verified: the twelve
PNGs regenerate byte-identical. Four of the entries are small faithful extracts of
larger deposited artifacts (`trainset_obs.csv`, `measured_IMR90_SEN_vs_PRO_*.npy`,
`gm21_atac_vs_k27ac.csv`); each script falls back to the full artifact when it is
present, so nothing diverges.

Full stand-alone captions, with n and dispersion type for every panel, are in
[figures/captions.md](figures/captions.md). Every error bar in the set is SEM
and is labelled as such on the panel and in its caption.

---

## Environment

Two lockfiles, because two environments were used:

- [`env/gpu-environment.lock.txt`](env/gpu-environment.lock.txt) — the GPU
  environment that produced every model result. **Read the honesty note in
  it:** the installer ran `pip -q`, so no transitive freeze was ever captured
  and the hosts are gone. gReLU 1.1.0, torch 2.8.0+cu128 and CUDA 12.8 are
  evidenced by the run logs; the rest is marked as evidenced, inferred or
  not captured. Nothing is invented.
- [`env/local-analysis.lock.txt`](env/local-analysis.lock.txt) — the CPU
  environment for the audit and all twelve figures. These versions **are**
  exact.

[`env/ENVIRONMENT.md`](env/ENVIRONMENT.md) documents seven things that cost
real time, so nobody rediscovers them: PEP 668 blocking system pip on Ubuntu
24.04 and the `--system-site-packages` venv it forces; UCSC being unreachable
from some cloud providers and the NCBI GRCh38 fallback; the 194-contig
RefSeq→UCSC rename that fallback requires; gReLU 1.1.0 ignoring
`GRELU_GENOMES_DIR` despite the docs; `checkpoint=False` raising through an
identity comparison; `Aggregate` needing `length_aggfunc` and
`return_seqs="none"` raising `KeyError`; and RefSeq GRCh38 carrying IUPAC
codes that survive uppercasing and are rejected by the allowed-base set.

---

## Provenance and auditing

This repository carries its own audit trail, which is unusual and deliberate.

- [`audit/NUMBER_AUDIT.md`](audit/NUMBER_AUDIT.md) — 72 claimed numbers checked
  against the file each came from: 60 matched, 4 mismatched, 1 not found, 4
  verified but misattributed or mislabelled.
- [`audit/accession_log_merged.csv`](audit/accession_log_merged.csv) — all 44
  accessions the project touched, each re-verified live against NCBI eutils or
  the ENCODE API on 2026-10-04, including the one that does not exist.
- [`audit/n_consistency_check.txt`](audit/n_consistency_check.txt) — every
  reported n re-added; `evaluable + ambiguous == total` holds in all eight AUC
  reports.
- [`audit_fixes/`](audit_fixes/) — the eight closed items, including two where
  the audit's own explanation turned out to be wrong and was replaced by a
  reproduced one.
- [`reports/ERRATA.md`](reports/ERRATA.md) — the five known report/data
  discrepancies.

## Citing

Please cite the manuscript _(link pending; not yet posted)_ and the Zenodo
deposition [10.5281/zenodo.23146107](https://doi.org/10.5281/zenodo.23146107)
for the large artifacts. The underlying GEO and ENCODE datasets remain under
their own terms and should be cited directly; all 44 are listed in
`audit/accession_log_merged.csv`.
