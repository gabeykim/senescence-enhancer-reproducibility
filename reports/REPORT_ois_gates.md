# OIS enhancer head — Gate 1 (motif additivity) and Gate 2 (the AP-1 anomaly)

Run date 2026-10-04 · RTX 3090 · gReLU 1.1.0 / Borzoi `human_rep0` trunk + trained
`head_ois.pt` · no re-caching, no retraining.

**Headline.** Gate 1: combinations do **not** reach unity and do **not** beat the best
single motif at matched load. The ceiling is **0.773× of the +0.861 requirement**, and the
dose response **peaks at 12 copies and then reverses**. Gate 2: the AP-1 anomaly is
**H1, context dependence** — the sign of AP-1's contribution flips with the activity of the
surrounding region (interaction p = 0.009), and genome-wide AP-1 sites are strongly
*enriched* where senescence H3K27ac is high. **H3 is refuted in the opposite direction.**

---

## 0. Artifact verification (before any compute)

| Check | Result |
|---|---|
| `region_fwd.npy` / `region_rc.npy` shape, dtype | `(111671, 1920)` float32, both |
| Non-finite values | 0 / 0 |
| All-zero rows | 0 / 0 (matches `sanitize_stats.json: all_zero_rows = 0`) |
| Rows where fwd ≠ rc | 111,671 / 111,671; mean abs(fwd − rc) = 0.0999 |
| Max abs embedding value | 45.76 |
| `head_ois.pt` | loads; `conv.weight [1,1920,1]`, `bias [1]`; all finite |

The cache was reused as-is. The 1.16 h of caching was not repeated. Environment rebuilt from
`setup_pod.sh` in ~8 min (NCBI GRCh38, 194 contigs renamed to UCSC names and reindexed,
symlinked into genomepy's default dir; TF32 off on both matmul and cudnn).

---

## 1. GATE 1 — motif additivity

### 1.1 How load was balanced (read this before the numbers)

The thing that would otherwise make combinations look good for trivial reasons is unequal
insertion load. Three controls:

1. **Total copy number is held constant** within each load level *K* ∈ {6, 12, 18, 20}.
   A single-motif arm places *K* copies of one motif; a pair places *K*/2 of each; the triple
   places *K*/3 of each. Every arm at a given *K* places exactly *K* motif instances.
2. **The slot grid is identical** across every arm at a given *K*: slot *j* occupies
   `span = 200 // K` bases starting at `j*span`, motif centred in its slot. Only the
   *identity* of the motif in each slot changes between arms. Positions are therefore not a
   confound.
3. **Equal copies is not equal bases.** NF-κB is 10 bp, C/EBPβ 8 bp, ETS 7 bp, so mixtures
   containing the shorter motifs insert *fewer* bases than an NF-κB-only arm at the same copy
   number (at *K* = 12: 120 bp for NF-κB alone vs 102 bp for NF-κB+ETS). **This biases against
   combinations, not toward them.** Inserted bases are reported per condition in
   `gate1_additivity_summary.csv`.

**Paired design.** One random 200 bp background per seed, reused by every arm at that seed, so
each arm's Δ is computed against its own seed's empty cassette. This removes
background-to-background variance (an improvement on the prior single-motif probe, which drew
an independent background per arm). Insertion methodology is otherwise identical to that probe
— same neutral context `chr1:265,888` (design_anchor), same 200 bp cassette, same background
composition (A .295 / C .205 / G .205 / T .295), same scoring (head over the trunk embedding
mean-pooled across the central 1 kb). **n = 10 background seeds per arm.**

### 1.2 Adjusted deltas at matched load (K = 12)

adjusted = real − scrambled, paired within seed. Required to move median → 95th percentile:
**+0.8613**.

| Condition | Arrangement | bp | raw Δ | scrambled Δ | **adjusted Δ** | ±SEM | p | frac of required |
|---|---|---|---|---|---|---|---|---|
| NF-κB + ETS | clustered | 102 | +0.7403 | +0.0746 | **+0.6657** | 0.0487 | 2.5e-07 | **0.773×** |
| NF-κB | single | 120 | +0.7420 | +0.1402 | **+0.6018** | 0.0490 | 6.3e-07 | 0.699× |
| C/EBPβ | single | 96 | +0.6911 | +0.1074 | +0.5837 | 0.0359 | 5.6e-08 | 0.678× |
| NF-κB + ETS | interleaved | 102 | +0.5584 | +0.0747 | +0.4838 | 0.0349 | 2.2e-07 | 0.562× |
| C/EBPβ + ETS | clustered | 90 | +0.4515 | +0.0299 | +0.4216 | 0.0653 | 1.2e-04 | 0.490× |
| NF-κB+C/EBPβ+ETS | clustered | 100 | +0.4717 | +0.0829 | +0.3889 | 0.0744 | 5.5e-04 | 0.452× |
| NF-κB + C/EBPβ | clustered | 108 | +0.4481 | +0.1054 | +0.3427 | 0.0729 | 1.1e-03 | 0.398× |
| C/EBPβ + ETS | interleaved | 90 | +0.2898 | −0.0026 | +0.2924 | 0.0838 | 6.8e-03 | 0.339× |
| NF-κB+C/EBPβ+ETS | interleaved | 100 | +0.3575 | +0.1113 | +0.2462 | 0.0405 | 1.8e-04 | 0.286× |
| ETS | single | 84 | +0.1562 | −0.0585 | +0.2147 | 0.0333 | 1.2e-04 | 0.249× |
| NF-κB + C/EBPβ | interleaved | 108 | +0.2023 | +0.1536 | +0.0487 | 0.0413 | 2.7e-01 | 0.057× |

Full table including *K* = 6 in `gate1_additivity_summary.csv`.

### 1.3 Does any combination reach or approach unity? No.

- **Best condition overall: 0.773× of the requirement.** Shortfall **+0.1956**.
- **The best combination does not beat the best single motif at matched load.**
  NF-κB+ETS clustered vs NF-κB alone, both at *K* = 12, paired across the same 10 backgrounds:
  **+0.0639 ± 0.0414, t = +1.54, p = 0.157, n = 10.** Not significant.
- So the honest statement is not merely "combinations fall short of unity" but **"combining
  motif families buys nothing measurable over simply using the best single family at the same
  insertion load."**

### 1.4 Additive, sub-additive, or saturating? Sub-additive, and then reversing.

Against a copy-share-scaled expectation from the singles (`gate1_additivity_vs_expected.csv`),
observed/expected spans **0.082–1.631, median 0.846, sub-additive in 9 of 16 combinations**.
That ratio is noisy and arrangement-dependent, so the dose curve is the more reliable read.

**The dose curve settles it (Gate 1b, n = 10 seeds per point):**

| Condition | K=6 | K=12 | K=18 | K=20 |
|---|---|---|---|---|
| NF-κB + ETS (clustered) | +0.3919 | **+0.6657** | +0.5444 | +0.3620 |
| NF-κB (single) | +0.4196 | **+0.6018** | +0.5686 | +0.3884 |
| C/EBPβ (single) | +0.1877 | **+0.5837** | +0.3585 | +0.2576 |
| NF-κB+C/EBPβ+ETS (clustered) | +0.3033 | **+0.3889** | +0.1455 | +0.2753 |
| ETS (single) | **+0.2277** | +0.2147 | +0.1541 | +0.1022 |

**Every condition peaks at K = 12 and declines beyond it.** The marginal gain per added copy
turns negative after K = 12 for all five (e.g. NF-κB+ETS: +0.046/copy over 6→12, then
−0.020/copy over 12→18 and −0.091/copy over 18→20). *K* = 20 is the hard ceiling for a 200 bp
cassette with 10 bp motifs — the cassette is then 100% motif, and the response collapses.

**This is saturation with an interior optimum, not a plateau.** The implication for directed
evolution is concrete: *the remaining headroom is not in adding more motif copies.* Simply
packing more consensus sites into a 200 bp element makes it worse past ~12 copies
(~50–60% motif occupancy). Reaching +0.861 requires something other than motif dosage —
longer elements, or sequence features the model values that a consensus-tiling design does not
produce.

### 1.5 Spacing matters, and clustered always wins

Interleaved vs clustered, paired within seed, adjusted for the matched scrambled arrangement:

| Condition | K | interleaved | clustered | difference | p |
|---|---|---|---|---|---|
| NF-κB + C/EBPβ | 12 | +0.0487 | +0.3427 | **−0.2940** | **0.00046** |
| NF-κB + ETS | 12 | +0.4838 | +0.6657 | **−0.1820** | **0.0087** |
| C/EBPβ + ETS | 12 | +0.2924 | +0.4216 | **−0.1292** | **0.047** |
| NF-κB+C/EBPβ+ETS | 12 | +0.2462 | +0.3889 | −0.1427 | 0.105 |
| NF-κB+C/EBPβ+ETS | 6 | +0.1566 | +0.3033 | −0.1467 | 0.166 |
| NF-κB + C/EBPβ | 6 | +0.2608 | +0.3169 | −0.0561 | 0.351 |
| C/EBPβ + ETS | 6 | +0.0748 | +0.1261 | −0.0514 | 0.289 |
| NF-κB + ETS | 6 | +0.3634 | +0.3919 | −0.0285 | 0.456 |

**Clustered beats interleaved in 8 of 8 comparisons**, significantly in 3 (all at *K* = 12).
Homotypic clustering — like motifs adjacent — is worth up to +0.29, which is larger than the
entire benefit of adding a second motif family. For NF-κB+C/EBPβ at *K* = 12 the arrangement
is the difference between a null result (+0.049, p = 0.27) and a clear one (+0.343, p = 0.001).
**Spacing is a bigger lever than composition** in this model, and any design work should treat
it as a first-class parameter.

---

## 2. GATE 2 — the AP-1 anomaly

### 2.1 Verdict

**The evidence supports H1 (context dependence). H3 is refuted, and in the opposite direction
to its prediction. H2 is not supported as a crude artifact, though it does show the random-
background insertion assay is a poor instrument for AP-1.**

### 2.2 H1 — context dependence: SUPPORTED

**H1a — saved ISM importance (n = 18 high-response held-out regions).** Attribution enrichment
at PWM-matched positions (rel ≥ 0.85), in-motif / out-motif mean |Δ|:

| Motif | n regions | mean enrichment | SD | p vs 1.0 |
|---|---|---|---|---|
| **AP-1 (FOS::JUN)** | 16 | **3.609** | 2.155 | **0.0002** |
| ETS | 17 | 2.740 | 1.962 | 0.0021 |
| NF-κB | 5 | 1.915 | 2.357 | 0.434 |
| C/EBPβ | 15 | 1.248 | 0.627 | 0.148 |

AP-1 carries the **highest** model attribution of all four motifs in real high-response
regions. **Caveat, stated plainly: these arrays store mean |Δ| over the three alternative
bases, so they are unsigned.** They establish that the model attends strongly to AP-1 in
context; they cannot give the direction. That is why H1b was added.

**H1b — signed ablation of native sites (the test that decides H1).** For real regions carrying
a strong native AP-1 match (rel ≥ 0.95, within ±100 bp of the design anchor), the matched bases
were replaced by a mononucleotide shuffle of themselves (length and composition preserved) and
the signed change in prediction read off. The motif-specific effect subtracts a same-length
shuffle at a motif-free position in the same window, removing any generic
"shuffling 10 bp perturbs the prediction" component.

| Set | n | ablate motif | ablate control | **motif-specific** | ±SEM | p |
|---|---|---|---|---|---|---|
| AP-1, **high**-response regions | 24 | −0.1133 | −0.0078 | **−0.1055** | 0.0463 | **0.032** |
| AP-1, **low**-response regions | 24 | +0.1065 | +0.0408 | **+0.0658** | 0.0425 | 0.135 |
| NF-κB, high-response (positive control) | 24 | −0.1221 | +0.0183 | **−0.1404** | 0.0330 | **0.0003** |

- **The sign of AP-1's contribution flips with context.** High vs low interaction:
  **difference −0.1713, Welch t = −2.73, p = 0.0091; Mann-Whitney p = 0.0076.**
- It is **continuous, not just a two-bin effect**: Spearman(region's measured GM21 H3K27ac
  response, AP-1-specific ablation effect) = **−0.3635, p = 0.011, n = 48**. The more active
  the region, the more positively AP-1 contributes.
- The **NF-κB positive control behaves exactly as it should** (−0.1404, p = 0.0003), which
  validates the ablation assay: it can detect a motif that contributes positively.
- Negative effect = ablating the motif *lowers* the prediction = the motif was contributing
  *positively* in that real context.

**H1c — is the negative insertion result just one odd background?** AP-1 insertion repeated in
four independent neutral contexts on four chromosomes (adjusted Δ at 8 copies, n = 10 seeds):

| Context | adjusted Δ | ±SEM | p |
|---|---|---|---|
| chr10:860,711 | −0.6158 | 0.0960 | 0.0001 |
| chr1:265,888 *(the prior probe's background)* | −0.5556 | 0.1449 | 0.0040 |
| chr11:7,513,386 | −0.4192 | 0.0820 | 0.0006 |
| chr12:5,468,210 | +0.0658 | 0.0999 | 0.526 |

Negative in three of four, so **the prior result was not an artifact of one background** — it
replicates in two further contexts. Note all four contexts are *low-activity* regions by
construction (|measured response| < 0.05), which is precisely the regime where H1b says AP-1
contributes negatively. H1c is therefore **consistent with H1**, not a counterexample: the
insertion probe has only ever been run in the context class where AP-1 is expected to be
negative.

### 2.3 H2 — insertion artifact: NOT supported as a crude artifact

Variant panel at 8 copies in the original context, adjusted against the scrambled TRE control,
n = 10 seeds each. Selected rows (full table in `gate2_h2_variants.csv`):

| Variant | Sequence | AP-1 PWM rel | adjusted Δ | ±SEM | p |
|---|---|---|---|---|---|
| TRE heptamer | **TGACTCA** | 0.983 | **−0.5556** | 0.1449 | 0.0040 |
| m4 C→G | **TGAGTCA** | 0.975 | **−0.5153** | 0.1433 | 0.0058 |
| m6 C→T | TGACTTA | 0.878 | −0.0225 | 0.0167 | 0.212 |
| m2 G→A | TAACTCA | 0.865 | −0.0738 | 0.0219 | 0.0082 |
| CRE octamer | TGACGTCA | 0.827 | **+0.1163** | 0.0245 | 0.0011 |
| m5 T→G | TGACGCA | 0.950 | **+0.3610** | 0.0501 | 0.0001 |
| m7 A→G | TGACTCG | 0.878 | **+0.4049** | 0.0474 | <0.0001 |

- **Only two variants are strongly negative, and they are exactly the two canonical TRE
  heptamers: TGACTCA and TGAGTCA** — i.e. the AP-1 consensus TGA(C/G)TCA. Every one of the
  other 21 variants is near zero or positive.
- **The effect does not grade with PWM affinity.** Across all 23 variants Spearman
  ρ = +0.004 (p = 0.986); Pearson r = −0.44 (p = 0.036) **driven entirely by the two extreme
  points**. Restricted to the 21 single-base mutants: **ρ = +0.13, p = 0.575** — no
  relationship. m5 T→G scores rel 0.950 on the PWM yet gives **+0.361**.
- **The insertion does not create confounding sites.** Mean per cassette at rel ≥ 0.90:
  NF-κB **0.00**, C/EBPβ 0.54, ETS 0.10. The one notable exception is informative — m5 T→G
  (TGACGCA) creates **3.0** C/EBPβ sites per cassette, which plausibly explains its large
  positive value.

**Interpretation, carefully.** By the stated H2 criterion — "varies sharply with small sequence
changes in ways that do not track AP-1 binding affinity" — the literal criterion is met. But
the *pattern* of that sharpness argues against a generic artifact: the two negative variants are
not arbitrary 7-mers, they are precisely the two canonical TRE forms, and the CRE octamer
TGACGTCA (bound by CREB/ATF rather than classic AP-1) behaves differently (+0.116). A model
that distinguishes TRE from CRE and from every single-base deviation is exhibiting **sharper
sequence recognition than MA0099.3 encodes**, not insensitivity to AP-1. The right conclusion
is that **the PWM is the wrong yardstick here**, and that the insertion assay measures
something real about the exact TRE — in a context class where that effect happens to be
negative (per H1b).

### 2.4 H3 — assay-specific biology: REFUTED, in the opposite direction

Matched within-study contrast in GM21 (GSE205898): H3K27ac RAS vs EV (**n = 2 vs 2**) and ATAC
RAS vs EV (**n = 4 vs 4**) over the **same 111,671 regions** and the **same perturbation**, so
the two readouts differ only in assay. Quantitative ATAC came from
`region_responses_ois.csv::GM21_ATAC_RAS_vs_EV` (the h5ad carries only the boolean
`passes_atac_filter`); the join was verified at corr = 1.000000 against the h5ad H3K27ac
values. Motif calls are PWM matches at rel ≥ 0.95 in a 300 bp window on the design anchor.
Overall corr(H3K27ac, ATAC) = +0.595 Pearson / +0.545 Spearman.

| Motif | n with | n without | ΔH3K27ac | p | ΔATAC | p | z-gap |
|---|---|---|---|---|---|---|---|
| **AP-1** | 16,909 | 94,762 | **+0.2061** | 2e-283 | +0.0549 | 9.8e-34 | **−0.1724** |
| NF-κB | 2,089 | 109,582 | +0.3194 | 4.5e-71 | +0.1869 | 1.5e-42 | −0.0842 |
| C/EBPβ | 3,002 | 108,669 | −0.0419 | 0.011 | +0.0360 | 5.3e-04 | +0.1197 |
| ETS | 4,978 | 106,693 | −0.1526 | 2.2e-61 | −0.0534 | 2.0e-12 | +0.1047 |

z-gap = z(ATAC) − z(H3K27ac), motif-present minus motif-absent; both readouts z-scored across
regions first, because raw ATAC and raw H3K27ac log2FCs are not on a common scale.

**H3 predicts AP-1 regions should show elevated ATAC with flat-or-negative H3K27ac — a positive
z-gap, and the largest of the four. Every part of that fails:**

- AP-1 regions have **strongly elevated H3K27ac** (+0.2061, p = 2e-283), not flat or negative.
- AP-1's z-gap is **negative (−0.1724) and the most negative of the four motifs** — AP-1 is the
  most *acetylation*-shifted motif here, the exact opposite of the prediction.
- **The specificity control matters and it passes the wrong way:** C/EBPβ and ETS have
  *positive* z-gaps. If anything is accessibility-shifted relative to acetylation in this data,
  it is those two, not AP-1.
- **Sharpest form of the test:** restricted to the 27,918 regions in the **top ATAC quartile**,
  AP-1 presence still predicts **higher** H3K27ac (+0.0482, p = 4.9e-05, n = 4,657 with vs
  23,261 without). The H3 mechanism would require the opposite sign.

### 2.5 AP-1 motif frequency in high- vs low-response real regions

Fraction of regions with ≥1 strong (rel ≥ 0.95) site, by decile of measured GM21 H3K27ac
response — n = 11,167–11,168 per decile:

| Decile | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| mean response | −1.17 | −0.71 | −0.48 | −0.29 | −0.11 | +0.09 | +0.30 | +0.56 | +0.89 | +1.48 |
| **AP-1** | .086 | .093 | .109 | .128 | .155 | .173 | .199 | .197 | .196 | .179 |
| NF-κB | .010 | .012 | .014 | .013 | .018 | .016 | .020 | .019 | .029 | .038 |
| C/EBPβ | .032 | .028 | .024 | .028 | .025 | .028 | .027 | .027 | .024 | .026 |
| ETS | .044 | .070 | .062 | .056 | .045 | .038 | .035 | .034 | .028 | .034 |

Top vs bottom decile (n = 11,167 vs 11,168):

| Motif | top | bottom | ratio | OR | p |
|---|---|---|---|---|---|
| **AP-1** | 0.1786 | 0.0861 | **2.07×** | **2.306** | 5.9e-94 |
| NF-κB | 0.0376 | 0.0098 | 3.82× | 3.929 | 8.0e-45 |
| C/EBPβ | 0.0258 | 0.0323 | 0.80× | 0.793 | 4.1e-03 |
| ETS | 0.0339 | 0.0441 | 0.77× | 0.762 | 1.1e-04 |

**AP-1 sites are strongly and monotonically enriched where senescence H3K27ac is high** — a
direct contradiction of what the insertion probe implies. The enrichment is **steeper for
H3K27ac than for ATAC** (AP-1 OR 2.31 for H3K27ac vs 1.45 for ATAC), again the opposite of H3.
The same pattern holds in held-out IMR90 (OR 1.85, p = 5.7e-56).

**A discordance worth flagging for the preprint:** C/EBPβ and ETS are *depleted* in
high-H3K27ac regions (OR 0.79, 0.76) yet both give *positive* insertion deltas in Gate 1. So
the insertion-vs-observational mismatch is **not unique to AP-1** — it affects three of four
motifs, in both directions. This is a general caution about reading insertion probes as claims
about motif biology.

---

## 3. What this means for design work

1. **The motif-dosage lever is exhausted at 0.773× of requirement.** It peaks at ~12 copies in
   a 200 bp cassette and reverses beyond that. Directed evolution over "how many consensus
   sites to insert" has essentially no headroom left.
2. **Combining families adds nothing beyond the best single family** at matched load
   (p = 0.157). Compositional diversity is not the missing ingredient.
3. **Arrangement is a real and under-explored lever.** Clustered beat interleaved 8/8, by up to
   +0.29 — larger than the gain from adding a second family. If directed evolution has room, it
   is in grammar and spacing, not in motif count or identity.
4. **Do not penalise AP-1 sites in a designed element.** The prior conclusion ("per this model
   adding AP-1 sites would make a designed element worse") holds only in low-activity
   backgrounds. In real high-response regions AP-1 contributes positively (p = 0.032), the
   model attends to it more than to any other motif tested (enrichment 3.61, p = 0.0002), and
   AP-1 sites are 2.07× enriched in the top H3K27ac decile (p = 5.9e-94).
5. **The random-background insertion probe is the weak instrument here.** It was only ever run
   in near-zero-response contexts. Its sign for AP-1 is a property of that context class, not a
   claim about AP-1. Any future probe should insert into contexts spanning the activity range.

## 4. Limitations

- **Replication depth in the underlying data is n = 2 vs 2 for H3K27ac in both studies**
  (GM21 RAS vs EV, IMR90 SEN vs PRO/QUI). ATAC is n = 4 vs 4. Every H3K27ac number inherits
  the n = 2 ceiling; the cross-dataset reproducibility ceiling of +0.637 is the relevant frame.
- Gate 1 / Gate 1b / H1c / H2 use **n = 10 background seeds** per arm; H1b uses **n = 24
  regions** per stratum; H1a uses **n = 18 regions** (16 with an AP-1 match). H3 uses the full
  111,671 regions.
- H1a's importance arrays are **unsigned**; direction comes only from H1b.
- H1b selects on a strong PWM match (rel ≥ 0.95) near the anchor; regions with weak or
  non-consensus AP-1 sites are not represented.
- H3's ATAC and H3K27ac come from different assays with different replicate counts and dynamic
  ranges; the z-scoring addresses scale but not all assay-specific bias. The four-motif control
  is what makes the AP-1 comparison interpretable, and it is reported above.
- All results are predictions from one model with one trunk (`human_rep0`). Nothing here is a
  wet-lab measurement.

## 5. Files

Scripts `gate_common.py`, `gate1_additivity.py`, `gate1b_ceiling.py`, `gate2_scan.py`,
`gate2_ap1.py`, `gate2_finish.py`, `gate2_h3.py`, `make_figures.py`, `verify_artifacts.py`.

Results `artifact_verification.json`, `gate1_*.csv/json`, `gate1b_*.csv/json`,
`gate2_motif_calls.csv` (111,671 × 30), `gate2_h1a_ism_enrichment.csv`,
`gate2_h1b_ablation.csv`, `gate2_h1c_contexts.csv`, `gate2_h1c_h2_raw.csv`,
`gate2_h2_variants.csv`, `gate2_h3_*.csv`, `gate2_ap1_results.json`, `gate2_h3_results.json`,
plus `GATE*_LOG.txt`.

Figures `fig1_gate1_additivity.png`, `fig2_gate1b_dose_ceiling.png`, `fig3_h1b_ablation.png`,
`fig4_h1c_h2.png`, `fig5_h3_atac_vs_k27ac.png`.

**Known issue in the record:** `gate2_ap1.py` raised `TypeError` in its final print block
(`r.sem` resolved to the pandas `Series.sem` *method* rather than the column). Every CSV had
already been written; `gate2_finish.py` recomputes the H2 summary and the H1b interaction test
from those CSVs. No GPU work was lost or repeated.
