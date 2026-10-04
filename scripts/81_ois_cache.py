#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 81_ois_cache.py
#
# Cache Borzoi trunk embeddings for all 111,671 regions over 8,459 windows. GPU STEP, ~1.16 h.
#
# CONSUMES: ois_enhancer_trainset.h5ad, /workspace/genomes/hg38/
# PRODUCES: ois_cache/{region_fwd.npy,region_rc.npy} (818 MB each, ZENODO), region_window_map.csv, windows.csv, preflight.json, sanitize_stats.json
# ---------------------------------------------------------------------------
"""
STEP 1 -- trunk cache for the OIS enhancer set, WINDOW-SHARED.

111,671 regions collapse to ~8,432 distinct 524 kb windows because Borzoi predicts a
196,608 bp span (6144 bins x 32 bp) from each input, and enhancers cluster. One forward
pass therefore serves ~13 regions. The cache stored is PER REGION (mean-pooled over that
region's bins), computed from PER WINDOW forward passes.

COORDINATE SYSTEM: caching uses region_center_hg38 -- the H3K27ac measurement coordinate,
which is where the target was measured. design_anchor_hg38 (the ATAC summit, median 245 bp
away) is used only by the generation probe, where the question is what sequence a TF reads.

BIN ARITHMETIC (the part that can silently go wrong):
  input is SEQ_LEN=524,288 bp; the predicted span is the central PRED_SPAN=196,608 bp,
  i.e. input coordinates [163840, 360448). For a window starting at genomic W and a region
  spanning [C-500, C+500):
      b0 = floor((C - 500 - W - 163840) / 32)
      b1 = ceil ((C + 500 - W - 163840) / 32)
  and 0 <= b0 < b1 <= 6144 must hold.
  For the REVERSE COMPLEMENT pass the output bins are reversed, so the same region occupies
  RC bins [6144 - b1, 6144 - b0).

Preflight aborts the run if the caching identity or the bin arithmetic fails.
"""
import json
import sys
import time
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import torch

import grelu.resources
import grelu.sequence.format
from grelu.model.heads import ConvHead

torch.backends.cudnn.allow_tf32 = False
torch.backends.cuda.matmul.allow_tf32 = False

SEQ_LEN = 524_288
PRED_SPAN = 196_608
N_BINS = 6144
BIN = 32
SPAN_OFF = (SEQ_LEN - PRED_SPAN) // 2      # 163840
HALF_REG = 500
C = 1920

TRAINSET = "/workspace/output/ois_enhancer_trainset/ois_enhancer_trainset.h5ad"
OUT = Path("/workspace/ois_cache")
OUT.mkdir(parents=True, exist_ok=True)
CHUNK = 200

ad = anndata.read_h5ad(TRAINSET)
var = ad.var.copy()
var["idx"] = np.arange(len(var))
n_reg = len(var)
print(f"regions {n_reg:,}  tasks {ad.shape[0]}  TF32 cudnn={torch.backends.cudnn.allow_tf32} "
      f"matmul={torch.backends.cuda.matmul.allow_tf32}", flush=True)

sizes = {}
for line in open("/workspace/genomes/hg38/hg38.fa.sizes"):
    a, b = line.split()
    sizes[a] = int(b)


# =====================================================================================
# WINDOW ASSIGNMENT (greedy covering, then clamp to chromosome, then re-verify)
# =====================================================================================
def assign_windows():
    wins, assign = [], np.full(n_reg, -1, np.int64)
    for chrom, g in var.groupby("chrom", sort=False):
        csz = sizes.get(chrom)
        if csz is None:
            continue
        g = g.sort_values("region_center_hg38")
        centers = g["region_center_hg38"].to_numpy()
        idxs = g["idx"].to_numpy()
        pending = list(range(len(centers)))
        while pending:
            i0 = pending[0]
            span_start = centers[i0] - HALF_REG
            wstart = span_start - SPAN_OFF
            wstart = max(0, min(wstart, csz - SEQ_LEN))     # clamp into the chromosome
            s_lo = wstart + SPAN_OFF
            s_hi = s_lo + PRED_SPAN
            wid = len(wins)
            took = []
            for k in pending:
                c = centers[k]
                if c - HALF_REG >= s_lo and c + HALF_REG <= s_hi:
                    assign[idxs[k]] = wid
                    took.append(k)
            if not took:                                     # region cannot be placed at all
                print(f"  WARNING unplaceable region idx={idxs[i0]} {chrom}:{centers[i0]}",
                      flush=True)
                pending.pop(0)
                continue
            wins.append((chrom, int(wstart)))
            tk = set(took)
            pending = [k for k in pending if k not in tk]
    return wins, assign


t0 = time.time()
windows, assign = assign_windows()
n_win = len(windows)
placed = int((assign >= 0).sum())
print(f"windows {n_win:,} for {placed:,}/{n_reg:,} regions  "
      f"({placed/max(n_win,1):.1f} regions/window)  [{time.time()-t0:.0f}s]", flush=True)

wdf = pd.DataFrame(windows, columns=["chrom", "wstart"])
wdf["wend"] = wdf["wstart"] + SEQ_LEN
var["window"] = assign

# bin ranges
b0 = np.full(n_reg, -1, np.int64)
b1 = np.full(n_reg, -1, np.int64)
for i in range(n_reg):
    w = assign[i]
    if w < 0:
        continue
    W = wdf["wstart"].iloc[w]
    cc = int(var["region_center_hg38"].iloc[i])
    b0[i] = (cc - HALF_REG - W - SPAN_OFF) // BIN
    b1[i] = -((-(cc + HALF_REG - W - SPAN_OFF)) // BIN)      # ceil div
var["bin0"], var["bin1"] = b0, b1
ok = (assign >= 0) & (b0 >= 0) & (b1 <= N_BINS) & (b1 > b0)
print(f"regions with valid bin range: {int(ok.sum()):,}/{n_reg:,}", flush=True)
if int(ok.sum()) != placed:
    print("BIN RANGE INCONSISTENT WITH ASSIGNMENT -- aborting", flush=True)
    sys.exit(1)

# ---- bin-arithmetic verification: two regions in one window must resolve to
# ---- different bin ranges consistent with their genomic offset
print("\n=== PREFLIGHT A: bin arithmetic ===", flush=True)
vv = var[ok]
multi = vv.groupby("window").filter(lambda g: len(g) >= 2).groupby("window")
checked, bad = 0, 0
for wid, g in list(multi)[:400]:
    g = g.sort_values("region_center_hg38")
    cs = g["region_center_hg38"].to_numpy()
    bs = g["bin0"].to_numpy()
    for a in range(len(g) - 1):
        d_gen = cs[a + 1] - cs[a]
        d_bin = (bs[a + 1] - bs[a]) * BIN
        checked += 1
        if abs(d_gen - d_bin) > BIN:                # within one bin of quantisation
            bad += 1
            if bad <= 3:
                print(f"  MISMATCH window {wid}: dgen={d_gen} dbin={d_bin}", flush=True)
        if bs[a + 1] == bs[a] and d_gen > BIN:
            bad += 1
print(f"  checked {checked:,} adjacent pairs sharing a window; inconsistent: {bad}",
      flush=True)
if bad:
    print("PREFLIGHT A FAILED", flush=True)
    sys.exit(1)
ex = list(multi)[0][1].sort_values("region_center_hg38").head(3)
for _, r in ex.iterrows():
    print(f"  eg {r['chrom']}:{int(r['region_center_hg38']):,} -> window {int(r['window'])} "
          f"bins [{int(r['bin0'])},{int(r['bin1'])})", flush=True)
print("PREFLIGHT A PASSED", flush=True)

# =====================================================================================
# MODEL
# =====================================================================================
model = grelu.resources.load_model(repo_id="Genentech/borzoi-model",
                                   filename="human_rep0.ckpt")
model.eval()
core = model.model.to("cuda")
print("model loaded", flush=True)

ALLOWED = set(grelu.sequence.format.ALLOWED_BASES)
_SAN = {c: "N" for c in map(chr, range(256)) if c.upper() not in ALLOWED}
san = {"windows_sanitized": 0, "chars_replaced": 0, "chars_seen": {}, "bases_scanned": 0}


def get_window_onehot(chrom, wstart):
    iv = pd.DataFrame({"chrom": [chrom], "start": [int(wstart)],
                       "end": [int(wstart) + SEQ_LEN], "strand": ["+"]})
    s = grelu.sequence.format.convert_input_type(iv, output_type="strings", genome="hg38")
    s = s[0] if isinstance(s, list) else s
    san["bases_scanned"] += len(s)
    extra = set(s) - ALLOWED
    if extra:
        san["windows_sanitized"] += 1
        for c in extra:
            k = s.count(c)
            san["chars_replaced"] += k
            san["chars_seen"][c] = san["chars_seen"].get(c, 0) + k
        s = s.translate(str.maketrans(_SAN))
    x = grelu.sequence.format.convert_input_type([s], output_type="one_hot")
    if not torch.is_tensor(x):
        x = torch.as_tensor(np.asarray(x))
    x = x.float()
    return x if x.dim() == 3 else x.unsqueeze(0)


def rc(x):
    return torch.flip(x, dims=[-2, -1])


# ---- PREFLIGHT B: caching identity ------------------------------------------------
print("\n=== PREFLIGHT B: caching identity head(full)[bins] == head(pooled bins) ===",
      flush=True)
probe = ConvHead(n_tasks=8, in_channels=C, act_func=None, pool_func="avg",
                 norm=False).to("cuda").eval()
devs = []
with torch.no_grad():
    for wid in [0, n_win // 3, 2 * n_win // 3, n_win - 1]:
        ch, ws = wdf["chrom"].iloc[wid], wdf["wstart"].iloc[wid]
        x = get_window_onehot(ch, ws).cuda()
        trunk = core.embedding(x)                     # (1, 1920, 6144)
        sub = var[ok & (var["window"] == wid)]
        if not len(sub):
            continue
        r = sub.iloc[0]
        lo, hi = int(r["bin0"]), int(r["bin1"])
        full = probe(trunk[:, :, lo:hi]).squeeze()            # head over the region's bins
        pooled = trunk[:, :, lo:hi].mean(dim=-1, keepdim=True)
        cached = probe(pooled).squeeze()
        den = full.abs().max().item()
        rel = (full - cached).abs().max().item() / (den if den > 0 else 1.0)
        devs.append(rel)
        print(f"  window {wid:5d} {ch} bins[{lo},{hi}) rel dev {rel:.3e}", flush=True)
worst = max(devs)
TOL = 1e-4
print(f"  worst {worst:.3e}  tolerance {TOL:.0e}", flush=True)
if not worst < TOL:
    print("PREFLIGHT B FAILED -- aborting", flush=True)
    sys.exit(1)
print("PREFLIGHT B PASSED", flush=True)
json.dump({"bin_arithmetic_pairs_checked": checked, "bin_arithmetic_bad": bad,
           "identity_devs": devs, "identity_worst": worst, "tolerance": TOL,
           "verdict": "PASS"}, open(OUT / "preflight.json", "w"), indent=2)
del probe
san.update({"windows_sanitized": 0, "chars_replaced": 0, "chars_seen": {}, "bases_scanned": 0})

# =====================================================================================
# CACHE
# =====================================================================================
fwd_p, rc_p = OUT / "region_fwd.npy", OUT / "region_rc.npy"
done_p = OUT / "chunks_done.txt"
mode = "r+" if fwd_p.exists() else "w+"
cf = np.lib.format.open_memmap(fwd_p, mode=mode, dtype=np.float32, shape=(n_reg, C))
cr = np.lib.format.open_memmap(rc_p, mode=mode, dtype=np.float32, shape=(n_reg, C))
done = set()
if done_p.exists():
    done = {int(x) for x in done_p.read_text().split() if x.strip()}
    print(f"resuming: {len(done)} chunks done", flush=True)

by_win = {}
for i in np.where(ok)[0]:
    by_win.setdefault(int(var["window"].iloc[i]), []).append(i)

t0 = time.time()
for c0 in range(0, n_win, CHUNK):
    cid = c0 // CHUNK
    if cid in done:
        continue
    c1 = min(c0 + CHUNK, n_win)
    with torch.no_grad():
        for wid in range(c0, c1):
            regs = by_win.get(wid)
            if not regs:
                continue
            x = get_window_onehot(wdf["chrom"].iloc[wid], wdf["wstart"].iloc[wid]).cuda()
            tf = core.embedding(x)[0]                    # (1920, 6144)
            tr = core.embedding(rc(x))[0]
            for i in regs:
                lo, hi = int(var["bin0"].iloc[i]), int(var["bin1"].iloc[i])
                cf[i] = tf[:, lo:hi].mean(dim=-1).cpu().numpy()
                cr[i] = tr[:, N_BINS - hi:N_BINS - lo].mean(dim=-1).cpu().numpy()
            del x, tf, tr
    cf.flush(); cr.flush()
    done.add(cid)
    done_p.write_text(" ".join(str(d) for d in sorted(done)))
    nd = len(done) * CHUNK
    el = time.time() - t0
    eta = (n_win - min(nd, n_win)) * (el / max(nd, 1)) / 3600
    print(f"  chunk {cid:3d}  windows {c1:6,}/{n_win:,}  elapsed {el/60:6.1f} min  "
          f"ETA {eta:5.2f} h", flush=True)

peak = torch.cuda.max_memory_allocated() / 1e9
tot = (time.time() - t0) / 3600
print(f"\nCACHE COMPLETE {tot:.2f} h  peak GPU {peak:.2f} GB", flush=True)
z = int((np.abs(np.asarray(cf)).sum(axis=1) == 0).sum())
print(f"  all-zero rows fwd {z} (expect {n_reg - placed})", flush=True)
print(f"  NaNs fwd {int(np.isnan(np.asarray(cf)).sum())}", flush=True)
print(f"  IUPAC: {san['windows_sanitized']} windows, {san['chars_replaced']} chars of "
      f"{san['bases_scanned']:,} scanned; codes {san['chars_seen']}", flush=True)
san.update(dict(cache_hours=tot, peak_gpu_gb=peak, n_windows=n_win, n_regions_cached=placed,
                all_zero_rows=z))
json.dump(san, open(OUT / "sanitize_stats.json", "w"), indent=2)
var[["chrom", "region_center_hg38", "design_anchor_hg38", "window", "bin0", "bin1",
     "split", "passes_atac_filter"]].to_csv(OUT / "region_window_map.csv")
wdf.to_csv(OUT / "windows.csv", index=False)
print("CACHE_DONE", flush=True)
