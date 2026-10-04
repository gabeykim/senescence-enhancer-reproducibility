# Provenance note — directional agreement between accessibility and H3K27ac change

**Authoritative value: 12,429 of 13,330 co-differential anchors agree in direction = 93.24%.**
**Authoritative file: `output/enhancer_test/step2_stats.json`.**

## The discrepancy

| Source | Count | Percentage |
|---|---|---|
| `output/enhancer_test/step2_stats.json` (`n_same_direction`) | 12,429 | 93.2408% of co-differential |
| `output/enhancer_test/STEP2_LOG.txt:81` (same run) | 12,429 (as "60.6% agree in direction") | — |
| `output/enhancer_test/REPORT_enhancer_cross_mechanism.md:225` | **12,432** | stated **93.3%** |
| Claim inventory / manuscript draft | — | stated **93.0%** |

## How it was resolved

Recomputed from the stored per-anchor table `output/enhancer_test/region_responses.csv`, using
the criterion in `scripts/52_enhancer_cross_mechanism.py:220-228`:

```
da   = both[both["gm21_atac"].abs() >= t]                    # differential in ATAC
conc = da["gm21_k27ac"].abs() >= t                           # also differential in H3K27ac
same = conc & (np.sign(da["gm21_k27ac"]) == np.sign(da["gm21_atac"]))
```

At the primary threshold t = 0.5, over the 54,498 anchors with both assays measured:

```
n_diff_atac        = 20,515
n_also_diff_k27ac  = 13,330   (65.0%  -> 35.0% show no H3K27ac change)
n_same_direction   = 12,429   (60.58% of all diff-ATAC; 93.2408% of the co-differential)
```

This reproduces `step2_stats.json` exactly.

## It is not a tie or boundary difference

Both candidate mechanisms were tested directly and both are ruled out:

- **Ties.** Exact-zero signs among the differential-ATAC rows: **0**. `np.sign` never returns 0
  here, so there is no tie for a different pass to resolve differently.
- **Boundary inclusion.** Re-running with strict `>` instead of `>=` changes **nothing**:
  n_diff 20,515, co-differential 13,330, same-direction 12,429 either way.

A third possibility — a different filter or a different anchor subset — would have to change
`n_diff_atac` or `n_also_diff_k27ac` as well, and the report quotes both of those at the
machine-written values (20,515 and 13,330).

## Conclusion

**12,432 is not reproducible from the stored data by any variant tested, and the generating
script never computes it.** The script also never computes 93.3% — only `frac_same_direction`
(60.6% of all differential-ATAC anchors) is written out. The "% of the co-differential" framing
exists only in the report, derived by hand from the larger count.

This is a report-side transcription error (12,429 → 12,432) that propagated into a derived
percentage. The claim inventory's 93.0% matches neither and appears to be a further rounding
drift.

**Use 12,429 / 13,330 = 93.24%.** The machine-written JSON and the run log agree with each
other and with a fresh recomputation from the per-anchor table.

## Secondary observation

At the stricter threshold |log2FC| >= 1.0 the agreement is **3,373 / 3,376 = 99.91%**
(n_diff_atac 5,339). That is a markedly stronger statement than the 0.5 threshold supports and
is already in `step2_stats.json`.

Verified 2026-10-04.
