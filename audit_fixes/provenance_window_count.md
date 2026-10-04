# Provenance note — 8,432 versus 8,459 distinct 524 kb windows

**Use 8,459 windows for 111,671 regions = 13.201 regions per window.**
**Authoritative file: `output/ois_enhancer_run/ois_cache/sanitize_stats.json`.**

## The two numbers

| Artifact | Regions | Windows | Ratio |
|---|---|---|---|
| `ois_enhancer_trainset/step2_stats.json`, `STEP2_LOG.txt` §6, `REPORT_ois_enhancer_trainset.md:14,20,279,368`, `window_covering_per_chrom.csv` (sums exactly) | 111,677 | **8,432** | 13.244 |
| `ois_enhancer_run/ois_cache/sanitize_stats.json`, corroborated by `windows.csv` (8,459 rows) and `region_window_map.csv` (111,671 rows) | 111,671 | **8,459** | 13.201 |

## The number audit's explanation was wrong

The audit attributed the gap to the covering being computed before 6 segmental-duplication
regions were dropped, with the cache re-deriving it. That cannot be the mechanism: **under a
fixed algorithm, removing regions can only hold the window count equal or lower.** Here the
count *rose* by 27 while the region count *fell* by 6.

## What was actually tested

Both covering algorithms were re-implemented from their source and run on the **same** region
set — the 111,671 cached regions, centres taken from `region_window_map.csv`:

| Algorithm | Source | Windows on the same 111,671 regions |
|---|---|---|
| Point-covering: a region counts as covered if its position falls in the predicted span | `scripts/71_ois_trainset.py:233-238` | **8,432** |
| Width-aware: the region's full 1 kb must lie inside the predicted span | `output/ois_enhancer_run/scripts/ois_cache.py:93` | **8,459** |

Each reproduces its reported number **exactly**, from an identical region set. The 6 dropped
regions contribute nothing to the difference.

## The real cause

The planning script advanced the greedy covering with

```python
i = np.searchsorted(pos, pos[i] + PRED_SPAN, side="left")   # 71_ois_trainset.py
```

treating each region as a dimensionless point: any region whose position fell inside the
196,608 bp predicted span was counted as served.

The executed cache requires the whole region to fit:

```python
if c - HALF_REG >= s_lo and c + HALF_REG <= s_hi:           # ois_cache.py, HALF_REG = 500
```

which shrinks the usable span from 196,608 bp to **195,608 bp** and costs 27 additional windows
genome-wide.

**Chromosome-end clamping is not involved.** `ois_cache.py:86` clamps `wstart` into
`[0, chrom_size − SEQ_LEN]`, but the width-aware run above was *unclamped* and already gives
exactly 8,459.

## Status of each number

- **8,432** is a *planning estimate*. It under-counts because it ignores region width. It set
  the GPU budget and appears in the trainset report in four places.
- **8,459** is what was executed and cached. Use it.

The "13.2× compression" headline holds either way (13.244 vs 13.201), but the two figures are
not the same quantity and must not be mixed — in particular, do not pair "111,671 regions" with
"8,432 windows".

Verified 2026-10-04.
