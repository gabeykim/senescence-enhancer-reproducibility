#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 90_gate_common.py
#
# Shared model loading, scoring and PWM utilities for the design gates. Imported, not run.
#
# CONSUMES: results/head_ois.pt, Borzoi trunk (HuggingFace)
# PRODUCES: none (module)
# ---------------------------------------------------------------------------
"""
Shared loading / scoring for the two follow-up gates.

Identical scoring contract to the prior probe (ois_probe.py):
  prediction = trained ConvHead applied to the Borzoi trunk embedding,
  mean-pooled over the central POOL_W=1000 bp of the predicted span.
Coordinate system is design_anchor_hg38 (ATAC summit), as in the probe.
"""
import numpy as np
import pandas as pd
import torch

import grelu.resources
import grelu.sequence.format
from grelu.model.heads import ConvHead

# the pod's known-good numerics settings
torch.backends.cudnn.allow_tf32 = False
torch.backends.cuda.matmul.allow_tf32 = False

DEV = "cuda"
RUN = "/workspace/ois_enhancer_run"
OUT = "/workspace/ois_gates"
TRAINSET = "/workspace/ois_enhancer_trainset/ois_enhancer_trainset.h5ad"
RESPONSES = "/workspace/ois_enhancer_trainset/region_responses_ois.csv"
GENOME_SIZES = "/workspace/genomes/hg38/hg38.fa.sizes"

SEQ_LEN = 524_288
PRED_SPAN = 196_608
BIN = 32
POOL_W = 1000
CASSETTE = 200

CB0 = (PRED_SPAN // 2 - POOL_W // 2) // BIN
CB1 = -((-(PRED_SPAN // 2 + POOL_W // 2)) // BIN)

# consensus strings exactly as used by the prior probe, so results are comparable
MOTIF_SEQ = {"NFKB_RELA": "GGGACTTTCC", "CEBPB": "TTGCGCAA",
             "AP1_FOSJUN": "TGACTCA", "ETS1": "AGGAAGT"}
JASPAR = {"NFKB_RELA": "MA0107.1", "CEBPB": "MA0466.2",
          "AP1_FOSJUN": "MA0099.3", "ETS1": "MA0098.3"}
BG = {"A": .295, "C": .205, "G": .205, "T": .295}

ALLOWED = set(grelu.sequence.format.ALLOWED_BASES)
_SAN = {c: "N" for c in map(chr, range(256)) if c.upper() not in ALLOWED}


def onehot(s):
    if set(s) - ALLOWED:
        s = s.translate(str.maketrans(_SAN))
    x = grelu.sequence.format.convert_input_type([s], output_type="one_hot")
    if not torch.is_tensor(x):
        x = torch.as_tensor(np.asarray(x))
    x = x.float()
    return x if x.dim() == 3 else x.unsqueeze(0)


def load_model():
    head = ConvHead(n_tasks=1, in_channels=1920, act_func=None, pool_func="avg",
                    norm=False).to(DEV)
    head.load_state_dict(torch.load(f"{RUN}/results/head_ois.pt", map_location=DEV))
    head.eval()
    model = grelu.resources.load_model(repo_id="Genentech/borzoi-model",
                                       filename="human_rep0.ckpt")
    model.eval()
    core = model.model.to(DEV)
    return head, core


def pick_batch(core):
    for b in (8, 6, 4, 3, 2, 1):
        try:
            x = torch.zeros(b, 4, SEQ_LEN, device=DEV)
            with torch.no_grad():
                core.embedding(x)
            del x
            torch.cuda.empty_cache()
            return b
        except RuntimeError:
            torch.cuda.empty_cache()
    return 1


def make_scorer(head, core, batch):
    def score(seqs):
        out = []
        for b0 in range(0, len(seqs), batch):
            x = torch.cat([onehot(s) for s in seqs[b0:b0 + batch]]).to(DEV)
            with torch.no_grad():
                t = core.embedding(x)[:, :, CB0:CB1].mean(dim=-1, keepdim=True)
                out.append(head(t).squeeze(-1).squeeze(-1).cpu().numpy())
            del x
        return np.concatenate(out)
    return score


def fetch(chrom, start, end):
    iv = pd.DataFrame({"chrom": [chrom], "start": [int(start)], "end": [int(end)],
                       "strand": ["+"]})
    s = grelu.sequence.format.convert_input_type(iv, output_type="strings", genome="hg38")
    s = s[0] if isinstance(s, list) else s
    return s.translate(str.maketrans(_SAN)) if (set(s) - ALLOWED) else s


def load_sizes():
    sizes = {}
    for line in open(GENOME_SIZES):
        a, b = line.split()
        sizes[a] = int(b)
    return sizes


def window_on(chrom, pos, sizes):
    """524 kb window centred on pos. Returns (seq, window_start, offset_of_pos_in_seq)."""
    ws = int(pos) - SEQ_LEN // 2
    ws = max(0, min(ws, sizes[chrom] - SEQ_LEN))
    return fetch(chrom, ws, ws + SEQ_LEN), ws, int(pos) - ws


def scramble(m, seed):
    a = list(m)
    np.random.default_rng(seed).shuffle(a)
    return "".join(a)


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


# ---------------------------------------------------------------- PWM utilities
def fetch_pfm(mid):
    import json
    import subprocess
    for ep in ["https://jaspar.elixir.no/api/v1/matrix/{}/?format=json",
               "https://jaspar.genereg.net/api/v1/matrix/{}/?format=json"]:
        try:
            r = subprocess.run(["curl", "-sS", "-f", "--max-time", "40", ep.format(mid)],
                               capture_output=True, text=True, timeout=60)
            if r.returncode == 0 and r.stdout.strip().startswith("{"):
                j = json.loads(r.stdout)
                if j.get("pfm"):
                    return j
        except Exception:
            continue
    return None


def build_pwms(names=None):
    """name -> dict(w=log-odds (4,L), mn=min score, mx=max score, id=, length=)"""
    names = names or list(JASPAR)
    pwms = {}
    for nm in names:
        j = fetch_pfm(JASPAR[nm])
        if j is None:
            raise RuntimeError(f"JASPAR fetch failed for {nm} {JASPAR[nm]}")
        pfm = np.array([j["pfm"][b] for b in "ACGT"], float) + 0.25
        pfm = pfm / pfm.sum(0, keepdims=True)
        w = np.log2(pfm / np.array([[BG[b]] for b in "ACGT"]))
        pwms[nm] = dict(w=w, mn=float(w.min(0).sum()), mx=float(w.max(0).sum()),
                        id=j.get("matrix_id", JASPAR[nm]), name=j.get("name"),
                        length=int(w.shape[1]))
    return pwms


_IDXB = {b: i for i, b in enumerate("ACGT")}


def seq_to_idx(s):
    return np.array([_IDXB.get(c, -1) for c in s.upper()], dtype=np.int64)


def scan_pwm(seq, pwm, rel_thresh=0.80, both_strands=True):
    """Return list of (start, rel_score, strand) for matches at >= rel_thresh."""
    w, mn, mx = pwm["w"], pwm["mn"], pwm["mx"]
    L = w.shape[1]
    hits = []
    for strand, s in ((("+"), seq), (("-"), revcomp(seq))) if both_strands else ((("+"), seq),):
        idx = seq_to_idx(s)
        n = len(idx) - L + 1
        if n <= 0:
            continue
        for i in range(n):
            sub = idx[i:i + L]
            if (sub < 0).any():
                continue
            sc = float(w[sub, np.arange(L)].sum())
            rel = (sc - mn) / (mx - mn)
            if rel >= rel_thresh:
                start = i if strand == "+" else len(seq) - i - L
                hits.append((start, rel, strand))
    return sorted(hits)


def scan_pwm_fast(seq, pwm, rel_thresh=0.80):
    """Vectorised version of scan_pwm; same return contract."""
    w, mn, mx = pwm["w"], pwm["mn"], pwm["mx"]
    L = w.shape[1]
    out = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        idx = seq_to_idx(s)
        n = len(idx) - L + 1
        if n <= 0:
            continue
        win = np.lib.stride_tricks.sliding_window_view(idx, L)      # (n, L)
        valid = (win >= 0).all(1)
        safe = np.where(win < 0, 0, win)
        sc = w[safe, np.arange(L)[None, :]].sum(1)
        rel = (sc - mn) / (mx - mn)
        ok = valid & (rel >= rel_thresh)
        for i in np.flatnonzero(ok):
            start = int(i) if strand == "+" else len(seq) - int(i) - L
            out.append((start, float(rel[i]), strand))
    return sorted(out)
