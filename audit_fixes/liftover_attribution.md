# Corrected liftover attribution

Three builds in this project perform a liftover. Each has its own figures. The
99.92% / 99.96% pair belongs to the **enhancer cross-mechanism test**, not to the OIS
training set.

| Build | Regions | Liftover | Midpoint success | Failures | Span integrity |
|---|---|---|---|---|---|
| **Enhancer cross-mechanism test** | 237,824 anchors | hg19 → hg38 | **237,628 / 237,824 = 99.92%** | 196 (0.08%) | both ends lift **4,996 / 5,000 = 99.92%**; of those, span within 50 bp of 1 kb **4,994 / 4,996 = 99.96%** |
| **Replicative ceiling test** | the same 237,824 anchors, reused | hg19 → **hg18** | **237,512 / 237,824 = 99.87%** | 312 (0.13%) | both ends lift 4,995 / 5,000; span within 50 bp **100.00%** |
| **OIS enhancer training set** | 112,099 K27ac summits | hg19 → hg38 | **112,078 / 112,099 = 99.98%** | 21 (0.02%) | **no span-integrity check was run** |

Sources:
- enhancer test — `output/enhancer_test/harmonization_stats.json`
  (`liftover_midpoint_success_pct 99.91758611410118`, `pct_span_preserved 99.95996797437951`);
  `output/enhancer_test/STEP1_LOG.txt:30-31`;
  `output/enhancer_test/REPORT_enhancer_cross_mechanism.md:83-85`
- replicative ceiling — `output/replicative_ceiling_test/quantification_stats.json`
  (`liftover_hg19_to_hg18_pct 99.8688105489774`, `liftover_hg19_to_hg18_failed 312`,
  `hg18_span_integrity_pct 100.0`); `output/replicative_ceiling_test/STEP1_LOG.txt`
- OIS trainset — `output/ois_enhancer_trainset/region_stats.json`
  (`liftover_pct 99.98126655902372`, `liftover_failed 21`);
  `output/ois_enhancer_trainset/STEP1_LOG.txt` §3;
  `REPORT_ois_enhancer_trainset.md:42`

## For the Methods section

1. **The training set's liftover is 99.98%, 21 failed.** Those 21 are the same 21 that appear
   in the loss breakdown (21 liftover + 154 window-before-chromosome-start + 247
   window-past-chromosome-end + 6 ambiguous = 428 = 112,099 − 111,671 = 0.38%).
2. **There is no 99.96%-equivalent for the OIS build** — a span-integrity check was never run
   on it. Either state that explicitly or run one. Do not carry the enhancer-test figure across.
3. **99.96% is conditional.** It is 4,994 of the 4,996 anchors whose *both ends* lifted, not
   99.96% of all anchors. The unconditional span figure for that build is 4,994 / 5,000 = 99.88%.

## Other Methods figures checked for the same error

| Figure | Correct build | Cross-contaminated? |
|---|---|---|
| 112,099 consensus regions → 111,671 surviving | OIS trainset | No |
| Loss breakdown 247 / 154 / 21 / 6, 0.38% | OIS trainset (`STEP2_LOG.txt` §5 and §7) | No |
| 237,824 anchors | Enhancer test; **reused verbatim** by the replicative ceiling test | No — but say "the same anchor set" where it appears under replicative results |
| Cross-dataset reproducibility +0.637 (n = 83,716, thr 0.5) | OIS trainset | No |
| Within-OIS ceiling +0.653 (n = 88,530, thr 0.5) | Enhancer test | No |
| **"35% show no H3K27ac change; 93.24% directional agreement"** | **Enhancer test** | **YES** — the claim inventory files it under "training set build". The OIS trainset's own ATAC concordance is an inverted statistic over a different region set: **25,082 / 57,169 = 43.9%** of differential-**H3K27ac** regions are also differential in ATAC, **40.8%** agreeing in direction (`ois_enhancer_trainset/STEP2_LOG.txt` §4). The two condition on opposite assays and must not be interchanged. |
| **Consensus-region criterion** | differs by build | State it. Enhancer test: 1 kb **midpoint**-centred anchors from a **combined ATAC + K27ac** consensus, ≥2 sample support, 500 bp dedup. OIS trainset: 1 kb **summit**-centred regions from **K27ac peaks only**, summit supported by ≥2 of 4 samples, 250 bp summit cluster. Already documented in `ois_enhancer_trainset/STEP2_LOG.txt` §2 as the reason +0.653 and +0.637 are not expected to be equal. |

Verified 2026-10-04.
