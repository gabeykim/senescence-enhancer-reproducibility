#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 32_calibration_setup.sh
#
# GPU env + NCBI GRCh38 staging for the expression calibration pass.
#
# CONSUMES: NCBI (network)
# PRODUCES: /workspace/genomes/hg38/{hg38.fa,hg38.fa.sizes}
# ---------------------------------------------------------------------------
# Environment setup for the calibration pass. Reuses the proven workarounds from
# prior sessions -- these are settled, not rediscovered:
#   - PEP 668 blocks system pip on Ubuntu 24.04 -> venv with --system-site-packages
#     so the RunPod-preinstalled torch (2.8.0+cu128) is reused rather than refetched
#   - UCSC unreachable from RunPod, and neither "hg38" (Ensembl calls it GRCh38.p14,
#     NCBI calls it GRCh38) matches the literal string gReLU searches for
#     -> explicit provider='NCBI', name='GRCh38', localname='hg38'
#   - NCBI ships RefSeq-style contig names (>1 NC_000001.11) -> remap to UCSC style
#     (chr1) via the assembly_report.txt UCSC-style-name column, then reindex
#   - gReLU 1.1.0 IGNORES GRELU_GENOMES_DIR -> symlink into genomepy's real default
#   - gReLU install takes ~25 min
set -euo pipefail

echo "=== [1/4] venv (reusing preinstalled torch) ==="
python3 -m venv /workspace/venv --system-site-packages
/workspace/venv/bin/python -c 'import torch; print("torch", torch.__version__, "cuda", torch.cuda.is_available())'

echo "=== [2/4] gReLU 1.1.0 (~25 min, run under nohup by the caller) ==="
/workspace/venv/bin/pip install --no-input gReLU==1.1.0

echo "=== [3/4] genome: NCBI GRCh38 -> UCSC contig names -> symlink ==="
mkdir -p /workspace/genomes
/workspace/venv/bin/python -c "
import genomepy
genomepy.install_genome(name='GRCh38', provider='NCBI',
                        genomes_dir='/workspace/genomes', localname='hg38',
                        annotation=False)
print('GENOME DOWNLOAD DONE')
"

/workspace/venv/bin/python3 -c "
mapping = {}
with open('/workspace/genomes/hg38/assembly_report.txt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if cols[9] != 'na':
            mapping[cols[0]] = cols[9]
n = 0
with open('/workspace/genomes/hg38/hg38.fa') as fin, \
     open('/workspace/genomes/hg38/hg38.chr.fa', 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            sid = line[1:].split()[0]
            fout.write('>' + mapping[sid] + '\n' if sid in mapping else line)
            n += sid in mapping
        else:
            fout.write(line)
print('renamed', n)
"
mv /workspace/genomes/hg38/hg38.fa /workspace/genomes/hg38/hg38.fa.ncbi_names
mv /workspace/genomes/hg38/hg38.chr.fa /workspace/genomes/hg38/hg38.fa
rm -f /workspace/genomes/hg38/hg38.fa.fai /workspace/genomes/hg38/hg38.fa.sizes
cd /workspace/genomes/hg38 && /workspace/venv/bin/python3 -c "
import pyfaidx
fa = pyfaidx.Fasta('hg38.fa', rebuild=True)
with open('hg38.fa.sizes','w') as f:
    for k in fa.keys():
        f.write(f'{k}\t{len(fa[k])}\n')
print('reindexed', len(fa.keys()), 'sequences')
"
rm -f /workspace/genomes/hg38/hg38.fa.ncbi_names
mkdir -p /root/.local/share/genomes
ln -sfn /workspace/genomes/hg38 /root/.local/share/genomes/hg38

echo "=== [4/4] verify ==="
/workspace/venv/bin/python3 -c "
import genomepy, pandas as pd, grelu.sequence.format
print('installed:', genomepy.list_installed_genomes())
iv = pd.DataFrame({'chrom':['chr1'],'start':[1000000],'end':[1000100],'strand':['+']})
s = grelu.sequence.format.convert_input_type(iv, output_type='strings', genome='hg38')
print('extraction OK, len =', len(s[0]))
"
grep -P '^chr1\t' /workspace/genomes/hg38/hg38.fa.sizes
echo "SETUP COMPLETE"
