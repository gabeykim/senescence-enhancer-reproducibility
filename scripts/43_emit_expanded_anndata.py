#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 43_emit_expanded_anndata.py
#
# Emit the multi-study AnnData and task table.
#
# CONSUMES: output/trainset_expanded/_X*.npy
# PRODUCES: output/trainset_expanded/senescence_trainset_multistudy.h5ad, multistudy_task_table.csv
# ---------------------------------------------------------------------------
"""
Stage 2: emit the expanded multi-study AnnData, matching the existing build's contract.

Contract preserved from output/trainset/senescence_trainset.h5ad:
  .X            (n_tasks, n_intervals) float32
  .var          chrom / start / end / strand / tss / gene_type / split
  intervals     524,288 bp, TSS-centered
  splits        chr9 + chr6 held out entirely as TEST; chr8/16/20 VAL; rest TRAIN
                (disjoint chromosome lists -- gReLU's split() drops val/test chroms
                 from train, so overlapping lists silently empty the train split)

What changes: .X now carries the batch-corrected target (see uns['target_transform']),
and the uncorrected values are kept in layers so nothing is destroyed.
"""
import json
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

ROOT = Path("/Users/gabeykim/Downloads/Senescence")
OUT = ROOT / "output" / "trainset_expanded"
TAB = ROOT / "output" / "tables"

TEST_CHROMS = ["chr9", "chr6"]
VAL_CHROMS = ["chr8", "chr16", "chr20"]

LOG = []


def R(s=""):
    print(s, flush=True)
    LOG.append(s)


def main():
    st = json.load(open(OUT / "_stage1.json"))
    genes = st["genes"]
    Xraw = np.load(OUT / "_Xraw.npy")
    Xc = np.load(OUT / "_Xc.npy")
    Xcs = np.load(OUT / "_Xcs.npy")
    meta = pd.read_csv(OUT / "multistudy_sample_assignments.csv")
    iv = pd.read_csv(OUT / "_intervals.csv", index_col=0).loc[genes]

    R("=" * 86)
    R("STAGE 2 -- EMIT EXPANDED ANNDATA")
    R("=" * 86)

    # ---------------- which samples become tasks ----------------
    # train_eligible -> training tasks; test_transfer / holdout_ris -> carried in the
    # object but flagged, so a training run can select on obs['role'] without needing a
    # second file. Excluded samples are dropped entirely.
    keep = meta["role"].isin(["train_eligible", "test_transfer", "holdout_ris"]).to_numpy()
    R(f"\n  samples in pooled matrix : {len(meta)}")
    R(f"  dropped (role='excluded'): {int((~keep).sum())}")
    R(f"  retained as tasks        : {int(keep.sum())}")
    meta = meta[keep].reset_index(drop=True)
    Xraw, Xc, Xcs = Xraw[keep], Xc[keep], Xcs[keep]

    # ---------------- splits ----------------
    chrom = iv["chrom"].to_numpy()
    split = np.where(np.isin(chrom, TEST_CHROMS), "test",
             np.where(np.isin(chrom, VAL_CHROMS), "val", "train"))
    iv = iv.assign(split=split)
    R(f"\n  split sizes: " + str(pd.Series(split).value_counts().to_dict()))
    tr = set(iv.loc[iv["split"] == "train", "chrom"])
    va = set(iv.loc[iv["split"] == "val", "chrom"])
    te = set(iv.loc[iv["split"] == "test", "chrom"])
    assert not (tr & va) and not (tr & te) and not (va & te), "chromosome lists overlap"
    R(f"  chromosome lists disjoint: OK  (train {len(tr)}, val {len(va)}, test {len(te)})")

    # ---------------- assemble ----------------
    obs = meta.set_index(meta["study"] + "::" + meta["sample"].astype(str))
    obs.index.name = "task"
    assert obs.index.is_unique, "task names are not unique"

    # Collapse the induction labels to the biologically distinct classes. GSE74324's
    # "HRASV12_OIS" and GSE205692's "HRASG12V_OIS" are the SAME mechanism (oncogenic RAS)
    # under two spellings of the same allele; counting them separately would overstate
    # mechanism diversity, which is precisely the quantity this build exists to improve.
    MECH_CLASS = {"HRASV12_OIS": "OIS", "HRASG12V_OIS": "OIS", "replicative": "replicative",
                  "etoposide_TIS": "TIS", "radiation_RIS": "RIS",
                  "contact_inhibition": "quiescence"}
    obs = obs.assign(mechanism_class=obs["mechanism"].map(MECH_CLASS))
    assert obs["mechanism_class"].notna().all(), "unmapped mechanism label"

    ad = anndata.AnnData(X=np.ascontiguousarray(Xcs, dtype=np.float32),
                         obs=obs, var=iv.copy())
    ad.layers["log_native"] = np.ascontiguousarray(Xraw, dtype=np.float32)
    ad.layers["study_centered"] = np.ascontiguousarray(Xc, dtype=np.float32)

    ad.uns["contract"] = {
        "seq_len": 524_288, "bin_size": 32, "label_len": 196_608,
        "genome": "hg38", "input_shape": "(4, 524288) float32 channels-first",
        "trunk_output": "(batch, 1920, 6144)",
        "head": "ConvHead + softplus (non-negative outputs)",
        "target_per_example": "(n_tasks, 1) float32",
        "NOTE_softplus": "The existing head's softplus enforces non-negative outputs. "
                         "The corrected targets here are CENTERED and therefore signed. "
                         "A softplus head cannot represent the negative half and MUST be "
                         "changed (e.g. linear output) before training on .X. This is a "
                         "required consequence of batch correction, not an oversight.",
    }
    ad.uns["target_transform"] = (
        "X = arm-balanced per-study per-gene centering of log expression, then divided by "
        "that study's pooled SD. layers['log_native'] = uncorrected log expression "
        "(log1p(TPM) for RNA-seq studies; log2 RMA x ln2 for the GSE205692 microarray). "
        "layers['study_centered'] = centered but not scaled.")
    ad.uns["splits"] = {"train": sorted(tr), "val": sorted(va), "test": sorted(te)}
    ad.uns["splits"] = {k: np.array(v, dtype=object) for k, v in ad.uns["splits"].items()}
    ad.uns["arm_rationale"] = json.dumps(
        json.load(open(TAB / "expansion_arm_rationale.json")))
    ad.uns["quantification_by_study"] = {
        "GSE175533": "TPM (deposited)",
        "GSE74324": "RPKM (deposited) -> TPM by exact per-sample rescale",
        "GSE99028": "raw counts (deposited) -> TPM via GENCODE v44 union-exon length",
        "GSE205692": "log2 RMA microarray intensity (deposited, sva-corrected by submitters)",
    }
    # h5py cannot store lists of dicts; JSON-encode the tabular bits
    ad.uns["symbol_harmonisation"] = json.dumps(st["symbol_harmonisation"])
    ad.uns["batch_correction"] = {
        "method": "arm-balanced per-study per-gene centering + per-study scale equalisation",
        "center_info": json.dumps(st["center_info"]),
        "scale_info": json.dumps(st["scale_info"]),
        "pca_stats_file": "output/trainset_expanded/batch_pca_stats.csv",
    }
    ad.uns["provenance"] = {
        "built": pd.Timestamp.today().strftime("%Y-%m-%d"),
        "scripts": "40_fetch_expansion_series.sh; 41_parse_expansion_metadata.py; 42_build_expanded_trainset.py; 43_emit_expanded_anndata.py",
        "supersedes": "output/trainset/senescence_trainset.h5ad",
        "reason": "single-study fine-tune failed its shuffled-label control (p=0.143)",
    }

    p = OUT / "senescence_trainset_multistudy.h5ad"
    ad.write_h5ad(p)
    R(f"\n  wrote {p}")
    R(f"  shape (n_tasks, n_intervals) = {ad.shape}")
    R(f"  X dtype {ad.X.dtype}  range [{ad.X.min():.3f}, {ad.X.max():.3f}]")
    R(f"  layers: {list(ad.layers.keys())}")

    # ---------------- replication accounting ----------------
    R("\n" + "=" * 86)
    R("REPLICATION ACCOUNTING -- is the effective n materially better than n~3?")
    R("=" * 86)
    trn = ad.obs[ad.obs["role"] == "train_eligible"]
    on = trn[trn["cls"] == "senescent"]
    off = trn[trn["cls"] == "proliferating"]
    R(f"\n  TRAINING tasks: {len(trn)}  (ON {len(on)}, OFF {len(off)}, "
      f"quiescent {int((trn['cls']=='quiescent').sum())})")
    R("\n  ON arm composition:")
    R(on.groupby(["study", "cell_line", "mechanism"]).size().to_string())
    R(f"\n  independent cell lines in the ON arm : {on['cell_line'].nunique()} "
      f"({sorted(on['cell_line'].unique())})")
    R(f"  independent studies in the ON arm    : {on['study'].nunique()}")
    R(f"  DISTINCT induction mechanisms        : {on['mechanism_class'].nunique()} "
      f"({sorted(on['mechanism_class'].unique())})")
    R("    (not 3 -- GSE74324 'HRASV12' and GSE205692 'HRASG12V' are the same oncogenic-RAS")
    R("     mechanism under two spellings, so OIS is represented twice, by two cell lines.)")
    R("\n  ON-arm tasks per (cell line, mechanism class) -- the replication that matters:")
    R(on.groupby(["cell_line", "mechanism_class"]).size().to_string())
    R("\n  HELD OUT (never trained on):")
    ho = ad.obs[ad.obs["role"] != "train_eligible"]
    R(ho.groupby(["role", "study", "cls"]).size().to_string())

    ad.obs.to_csv(OUT / "multistudy_task_table.csv")
    (OUT / "EMIT_LOG.txt").write_text("\n".join(LOG))


if __name__ == "__main__":
    main()
