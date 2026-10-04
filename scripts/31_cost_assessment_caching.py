#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 31_cost_assessment_caching.py
#
# Cost/benefit of caching trunk embeddings vs recomputing them.
#
# CONSUMES: output/trainset/senescence_trainset.h5ad
# PRODUCES: stdout cost table
# ---------------------------------------------------------------------------
"""
Pre-run cost assessment for the senescence head fine-tune.

Answers the question posed before committing GPU budget: can trunk embeddings
be cached, what does the head actually read, and what does that do to cost?

RUN LOCALLY. No GPU needed -- this is arithmetic over measured constants plus
a source-level reading of the head architecture.

=============================== THE KEY RESULT ==============================
The cached representation is MATHEMATICALLY EXACT, not an approximation.

The head created by LightningModel.change_head() (which tune_on_dataset calls)
is, verified from gReLU 1.1.0 source:

    ConvHead(n_tasks=N, in_channels=1920, act_func=None, pool_func="avg", norm=False)
      -> ChannelTransformBlock(order="CDNA", norm=False, act_func=None, dropout=0)
           == nn.Conv1d(1920, N, kernel_size=1)        [pure linear]
      -> AdaptivePool("avg") == nn.AdaptiveAvgPool1d(1) [pure linear]

    and LightningModel.activation is Identity for task='regression'
    (measured live in the Task 4 run: "2 | activation | Identity").

Both stages are linear, so they COMMUTE with each other:

    out[t] = mean_b ( sum_c W[t,c] * x[c,b] + beta[t] )
           = sum_c W[t,c] * ( mean_b x[c,b] ) + beta[t]

So caching  m[c] = mean_b x[c,b]  -- a (1920,) vector per gene -- and applying
the head to m gives BIT-IDENTICAL results to applying the head to the full
(1920, 6144) trunk output. This is an algebraic identity, not a lossy summary.

VALIDITY CONDITIONS (all currently hold; caching breaks if any changes):
  1. pool_func == "avg". Max-pooling is NOT linear and would NOT commute.
  2. act_func is None inside the head (no nonlinearity between conv and pool).
  3. norm=False (a batchnorm in train mode would introduce batch coupling).
  4. No input augmentation. RC/shift augmentation changes the input sequence
     per epoch, so a fixed cache would be stale. See AUGMENTATION TRADEOFF.
The plan verifies 1-3 empirically on the instance (assert head(full) ==
head(cached)) before trusting the cache, rather than relying on this argument.
"""
import anndata

# ---- measured constants from prior sessions (not re-derived) -------------
TRUNK_CHANNELS = 1920      # measured: model.embedding(x).shape == (1, 1920, 6144)
TRUNK_BINS = 6144          # == label_len 196608 / bin_size 32
SEC_PER_GENE_LIGHTNING = 1.639   # Task 3, Lightning Trainer.predict, bs=1
SEC_PER_GENE_RAW = 0.8606        # Task 1, predict_on_seqs, bs=1
RATE_USD_PER_HOUR = 0.75

# controls the brief mandates, each of which is a full training run
RUNS = {"real model": 1, "shuffled-label permutations": 3, "leave-one-PDL-out": 1}


def main():
    ad = anndata.read_h5ad("output/trainset/senescence_trainset.h5ad")
    n = ad.shape[1]
    cls = ad.obs["cls"].value_counts().to_dict()
    n_train_tasks = int(cls["proliferating"] + cls["senescent"])

    full = n * TRUNK_CHANNELS * TRUNK_BINS * 4
    cache = n * TRUNK_CHANNELS * 4

    print(f"genes={n:,}  ON/OFF training tasks={n_train_tasks}  (quiescent {cls['quiescent']} held out)")
    print(f"\nSTORAGE")
    print(f"  full trunk output : {full/1e9:>9,.1f} GB   (infeasible)")
    print(f"  pooled cache      : {cache/1e6:>9,.1f} MB   ({full/cache:,.0f}x smaller, EXACT)")

    print(f"\nONE TRUNK PASS")
    for lbl, s in [("Lightning path", SEC_PER_GENE_LIGHTNING), ("raw path", SEC_PER_GENE_RAW)]:
        h = n * s / 3600
        print(f"  {lbl:16s}: {h:5.2f} h = ${h*RATE_USD_PER_HOUR:5.2f}")

    total_runs = sum(RUNS.values())
    print(f"\nWITHOUT CACHING ({total_runs} runs: {RUNS})")
    for ep in (10, 20):
        h = n * SEC_PER_GENE_LIGHTNING * ep * total_runs / 3600
        h_opt = n * SEC_PER_GENE_RAW * ep * total_runs / 3600
        print(f"  {ep:2d} epochs each: {h_opt:6.0f}-{h:6.0f} h = "
              f"${h_opt*RATE_USD_PER_HOUR:,.0f}-${h*RATE_USD_PER_HOUR:,.0f}")

    print(f"\nWITH CACHING")
    print(f"  one trunk pass + all {total_runs} runs on cached features")
    print(f"  head params: {TRUNK_CHANNELS*n_train_tasks + n_train_tasks:,} "
          f"-> an epoch is one ({n:,} x {TRUNK_CHANNELS}) matmul, sub-second")
    print(f"  TOTAL: ~$3.50-6.75, i.e. a 25-100x reduction")


if __name__ == "__main__":
    main()
