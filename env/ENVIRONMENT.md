# Environment quirks

Seven things cost real time on this project. Each is recorded here so nobody
rediscovers them. All were hit on Ubuntu 24.04 / Python 3.12 / gReLU 1.1.0 /
torch 2.8.0+cu128 on rented RTX 4090 and RTX 3090 instances.

### 1. PEP 668 blocks system pip on Ubuntu 24.04

`pip install` into the system interpreter fails with
`error: externally-managed-environment`. A plain venv then *loses* the base
image's CUDA-matched torch and reinstalls it from scratch (~10 min of wheel
download for a package already present).

**Fix:** create the venv with system packages visible, so the existing
cu128 torch is reused:

```bash
python3 -m venv --system-site-packages /workspace/venv
```

### 2. UCSC is unreachable from some cloud providers

`genomepy.install_genome(..., provider="UCSC")` hangs or fails from the RunPod
network. The genome has to come from NCBI instead:

```python
genomepy.install_genome(name="GRCh38", provider="NCBI",
                        genomes_dir="/workspace/genomes", localname="hg38",
                        annotation=False)
```

### 3. NCBI GRCh38 needs a 194-contig rename to UCSC names

The NCBI assembly uses RefSeq accessions (`NC_000001.11`), but every coordinate
in this project is UCSC-style (`chr1`). Column 10 of `assembly_report.txt` holds
the UCSC name; rewrite the FASTA headers with it and reindex. Exactly **194**
contigs are renamed. See `scripts/80_setup_pod.sh` for the working code.
Reindex with `pyfaidx.Fasta(..., rebuild=True)` afterwards or the stale `.fai`
will silently mismatch.

### 4. gReLU 1.1.0 ignores `GRELU_GENOMES_DIR`

The documentation says the genome directory is configurable by that environment
variable. In 1.1.0 it is not read. The genome must be symlinked into genomepy's
real default location:

```bash
mkdir -p /root/.local/share/genomes
ln -sfn /workspace/genomes/hg38 /root/.local/share/genomes/hg38
```

### 5. `checkpoint=False` raises a spurious error

Passing `checkpoint=False` trips an identity comparison (`is True`) rather than a
truth test, so the false branch raises instead of disabling checkpointing. Avoid
the high-level wrapper and call the trunk directly — every script here uses
`model.model.embedding(x)`.

### 6. Two more gReLU API traps

- The `Aggregate` utility fails on Borzoi's bin count unless `length_aggfunc` is
  given explicitly.
- `return_seqs="none"` raises `KeyError` instead of returning no sequences.

Both are avoided by calling the trunk directly, as above.

### 7. RefSeq GRCh38 contains IUPAC codes that survive uppercasing

The FASTA carries `R Y B M K S W` ambiguity codes. They are not lowercase
soft-masking, so `.upper()` does not remove them, and gReLU's one-hot encoder
rejects them against its allowed-base set. Every sequence must be sanitised to
`N` before encoding. Measured impact, not assumed:

| Run | Windows/genes touched | Characters replaced |
|---|---|---|
| Enhancer trunk cache (8,459 windows) | 33 | 98 |
| Expression cache (19,652 genes) | 49 | 133 |
| Expanded cache (13,211 genes) | 31 | 70 |

### Also worth knowing

- **TF32 must be off** for reproducible numerics. Every GPU script sets both
  `torch.backends.cudnn.allow_tf32 = False` and
  `torch.backends.cuda.matmul.allow_tf32 = False`.
- **The head needs `act_func=None`.** Targets are signed log2 fold changes; a
  softplus head cannot represent the negative half.
- **`pyarrow` is absent** from the gReLU dependency set, so `DataFrame.to_parquet`
  fails with "Unable to find a usable engine". Write CSV.
- **Backgrounding over SSH hangs the client.** `cmd &` inside `ssh host "..."`
  leaves the session open even with output redirected; use
  `setsid nohup cmd </dev/null >log 2>&1 &`.
- **`pkill -f <pattern>` kills its own shell** when the pattern also matches the
  command line that launched it. Bracket the first character, or use a separate
  invocation.
