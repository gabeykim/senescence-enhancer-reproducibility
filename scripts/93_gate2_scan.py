#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 93_gate2_scan.py
#
# Genome-wide PWM scan of all 111,671 regions at four stringencies. CPU, ~30 s.
#
# CONSUMES: ois_enhancer_trainset.h5ad, region_responses_ois.csv, hg38
# PRODUCES: results/gate2_motif_calls.csv (20 MB), gate2_scan_meta.json
# ---------------------------------------------------------------------------
"""
GATE 2, step 1 (CPU) -- motif calls for every region in the trainset.

Fetches a WIN bp window centred on design_anchor_hg38 (the ATAC summit; falls back
to region_center_hg38 where no anchor exists) for all 111,671 regions and scans it
with the four JASPAR PWMs used by the prior probe. Writes a per-region table of
match counts and best relative scores, which H1 (native-site ablation) and H3
(ATAC vs H3K27ac) both consume.

Scanning is vectorised across regions: for each PWM column j we gather the column's
base scores for every region at once, so the whole genome-wide scan is L gathers per
motif per strand rather than a Python loop over regions.
"""
import json
import sys
import time
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/ois_gates/scripts")
from gate_common import (JASPAR, OUT, RESPONSES, TRAINSET, build_pwms, fetch,
                         revcomp)

WIN = 300                 # +-150 bp around the anchor: the element core
THRESHOLDS = [0.80, 0.85, 0.90, 0.95]   # 0.80 is the prior probe's; it calls a
                                        # motif in >90% of regions, so stricter
                                        # cuts are recorded for the H3 contrast
REL_THRESH = 0.80
CHUNK = 20_000
_IDXB = {b: i for i, b in enumerate("ACGT")}


def main():
    t0 = time.time()
    outdir = Path(OUT) / "results"
    outdir.mkdir(parents=True, exist_ok=True)

    ad = anndata.read_h5ad(TRAINSET)
    var = ad.var.reset_index(drop=True)
    tasks = list(ad.obs.index)
    X = np.asarray(ad.X, np.float32)
    n = len(var)
    print(f"{n} regions, tasks {tasks}", flush=True)

    anchor = np.where(var["design_anchor_hg38"].to_numpy() > 0,
                      var["design_anchor_hg38"].to_numpy(),
                      var["region_center_hg38"].to_numpy())
    half = WIN // 2

    # ---- fetch all windows -----------------------------------------------------
    seqs = []
    for i in range(0, n, CHUNK):
        sl = slice(i, min(i + CHUNK, n))
        iv = pd.DataFrame({"chrom": var["chrom"].astype(str).to_numpy()[sl],
                           "start": anchor[sl] - half,
                           "end": anchor[sl] - half + WIN,
                           "strand": "+"})
        import grelu.sequence.format
        s = grelu.sequence.format.convert_input_type(iv, output_type="strings",
                                                     genome="hg38")
        seqs.extend(list(s))
        print(f"  fetched {min(i+CHUNK,n)}/{n}  ({time.time()-t0:.0f}s)", flush=True)
    assert len(seqs) == n
    lens = np.array([len(s) for s in seqs])
    print(f"window lengths: min {lens.min()} max {lens.max()} "
          f"(expected {WIN}); n != {WIN}: {(lens != WIN).sum()}", flush=True)

    # ---- encode to an (n, WIN) int8 matrix, -1 for N/other ---------------------
    idx = np.full((n, WIN), -1, dtype=np.int8)
    for i, s in enumerate(seqs):
        s = s.upper()[:WIN]
        a = np.frombuffer(s.encode(), dtype=np.uint8)
        col = np.full(WIN, -1, dtype=np.int8)
        for b, v in _IDXB.items():
            col[:len(a)][a == ord(b)] = v
        idx[i] = col
    n_amb = int((idx < 0).sum())
    print(f"ambiguous/N bases across all windows: {n_amb} "
          f"({100*n_amb/(n*WIN):.4f}%)", flush=True)

    pwms = build_pwms(list(JASPAR))
    for nm, p in pwms.items():
        print(f"  PWM {nm}: {p['id']} ({p['name']}) length {p['length']}", flush=True)

    # ---- vectorised scan -------------------------------------------------------
    def scan_all(w):
        """w: (4, L) log-odds. Returns (n, n_pos) best-of-both-strands rel scores."""
        L = w.shape[1]
        npos = WIN - L + 1
        wrc = w[::-1, ::-1].copy()                       # reverse complement PWM
        best = np.full((n, npos), -np.inf, dtype=np.float32)
        for mat in (w, wrc):
            sc = np.zeros((n, npos), dtype=np.float32)
            bad = np.zeros((n, npos), dtype=bool)
            for j in range(L):
                sub = idx[:, j:j + npos]
                bad |= sub < 0
                lut = np.append(mat[:, j].astype(np.float32), 0.0)   # index 4 -> pad
                sc += lut[np.where(sub < 0, 4, sub)]
            sc[bad] = -np.inf
            np.maximum(best, sc, out=best)
        return best

    rows = {}
    hitrates = {}
    for nm, p in pwms.items():
        sc = scan_all(p["w"])
        rel = (sc - p["mn"]) / (p["mx"] - p["mn"])
        rel[~np.isfinite(sc)] = -np.inf
        finite = np.where(np.isfinite(rel), rel, -np.inf)
        rows[f"{nm}_best_rel"] = finite.max(1).astype(np.float32)
        rows[f"{nm}_best_pos"] = finite.argmax(1).astype(np.int32)
        hitrates[nm] = {}
        for th in THRESHOLDS:
            hit = rel >= th
            tag = f"{nm}_count_{int(th*100)}"
            rows[tag] = hit.sum(1).astype(np.int32)
            frac = float((hit.sum(1) > 0).mean())
            hitrates[nm][str(th)] = dict(regions_with_hit=int((hit.sum(1) > 0).sum()),
                                         frac=frac,
                                         mean_hits=float(hit.sum(1).mean()))
            print(f"  {nm} @{th:.2f}: regions with >=1 hit "
                  f"{int((hit.sum(1)>0).sum())} ({100*frac:.2f}%)  "
                  f"mean hits/region {hit.sum(1).mean():.3f}", flush=True)
            del hit
        rows[f"{nm}_count"] = rows[f"{nm}_count_80"]
        print(f"    ({time.time()-t0:.0f}s)", flush=True)
        del sc, rel, finite

    out = pd.DataFrame(rows)
    out.insert(0, "region_index", np.arange(n))
    out.insert(1, "chrom", var["chrom"].astype(str).to_numpy())
    out.insert(2, "anchor", anchor)
    out["passes_atac_filter"] = var["passes_atac_filter"].to_numpy()
    out["design_anchor_dist"] = var["design_anchor_dist"].to_numpy()
    out["split"] = var["split"].astype(str).to_numpy()
    for t in tasks:
        out[t] = X[tasks.index(t)]

    # matched quantitative ATAC lives only in the companion CSV
    resp = pd.read_csv(RESPONSES)
    print(f"responses csv: {resp.shape}, cols {list(resp.columns)}", flush=True)
    key_a = var["chrom"].astype(str) + ":" + var["summit_hg19"].astype(str)
    resp_key = resp["chrom"].astype(str) + ":" + resp["summit_hg19"].astype(str)
    m = resp.set_index(resp_key)
    m = m[~m.index.duplicated(keep="first")]
    aligned = m.reindex(key_a.to_numpy())
    out["GM21_ATAC_RAS_vs_EV"] = aligned["GM21_ATAC_RAS_vs_EV"].to_numpy()
    out["atac_summit_dist"] = aligned["atac_summit_dist"].to_numpy()
    matched = int(np.isfinite(out["GM21_ATAC_RAS_vs_EV"]).sum())
    print(f"regions with matched quantitative ATAC: {matched}/{n}", flush=True)

    # sanity: the h5ad H3K27ac values should agree with the CSV for the same key
    chk = np.corrcoef(np.nan_to_num(aligned["GM21_RAS_vs_EV"].to_numpy()),
                      np.nan_to_num(out["GM21_RAS_vs_EV"].to_numpy()))[0, 1]
    print(f"join sanity: corr(h5ad GM21, csv GM21) = {chk:.6f}", flush=True)

    out.to_csv(outdir / "gate2_motif_calls.csv", index=False)
    json.dump(dict(window=WIN, rel_thresh=REL_THRESH, thresholds=THRESHOLDS,
                   hit_rates=hitrates, n_regions=int(n),
                   ambiguous_bases=n_amb, matched_atac=matched,
                   join_corr_check=float(chk), seconds=time.time() - t0,
                   motif_sources={k: dict(id=v["id"], name=v["name"],
                                          length=v["length"]) for k, v in pwms.items()},
                   regions_with_hit={k: int((rows[f"{k}_count"] > 0).sum()) for k in pwms}),
              open(outdir / "gate2_scan_meta.json", "w"), indent=2)
    print(f"SCAN_DONE in {(time.time()-t0)/60:.1f} min -> gate2_motif_calls.csv",
          flush=True)


if __name__ == "__main__":
    main()
