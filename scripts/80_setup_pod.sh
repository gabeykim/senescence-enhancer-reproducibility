#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 80_setup_pod.sh
#
# GPU env build: venv, gReLU 1.1.0, NCBI GRCh38 with the 194-contig UCSC rename. GPU STEP.
#
# CONSUMES: NCBI (network)
# PRODUCES: /workspace/venv, /workspace/genomes/hg38/{hg38.fa,hg38.fa.sizes,hg38.fa.fai}
# ---------------------------------------------------------------------------
# Full env + genome setup. Known workarounds baked in:
#  - PEP 668 on this image -> venv with --system-site-packages (reuses the system torch,
#    which is already cu128 and sees the 4090; a fresh torch install would waste ~10 min)
#  - gReLU 1.1.0 IGNORES GRELU_GENOMES_DIR -> genome must be symlinked into genomepy's
#    real default dir, /root/.local/share/genomes
#  - UCSC hgdownload is unreachable from RunPod -> install NCBI GRCh38 and rename its
#    RefSeq contigs to UCSC names using assembly_report.txt column 10, then reindex
set -euo pipefail
cd /workspace

python3 -m venv --system-site-packages /workspace/venv
/workspace/venv/bin/pip -q install --upgrade pip
/workspace/venv/bin/pip -q install gReLU==1.1.0
echo "GRELU_INSTALLED"
/workspace/venv/bin/python -c "import grelu, torch; print('grelu', grelu.__version__, '| torch', torch.__version__, '| cuda', torch.cuda.is_available())"

mkdir -p /workspace/genomes
/workspace/venv/bin/python -c "
import genomepy
genomepy.install_genome(name='GRCh38', provider='NCBI',
                        genomes_dir='/workspace/genomes', localname='hg38',
                        annotation=False)
print('GENOME_DOWNLOAD_DONE')
"
/workspace/venv/bin/python -c "
mapping = {}
with open('/workspace/genomes/hg38/assembly_report.txt') as f:
    for line in f:
        if line.startswith('#'): continue
        c = line.rstrip('\n').split('\t')
        if len(c) > 9 and c[9] != 'na':
            mapping[c[0]] = c[9]
n = 0
with open('/workspace/genomes/hg38/hg38.fa') as fin, open('/workspace/genomes/hg38/hg38.chr.fa','w') as fout:
    for line in fin:
        if line.startswith('>'):
            sid = line[1:].split()[0]
            if sid in mapping:
                fout.write('>' + mapping[sid] + '\n'); n += 1
            else:
                fout.write(line)
        else:
            fout.write(line)
print('renamed', n)
"
mv /workspace/genomes/hg38/hg38.fa /workspace/genomes/hg38/hg38.fa.ncbi
mv /workspace/genomes/hg38/hg38.chr.fa /workspace/genomes/hg38/hg38.fa
rm -f /workspace/genomes/hg38/hg38.fa.fai /workspace/genomes/hg38/hg38.fa.sizes
cd /workspace/genomes/hg38
/workspace/venv/bin/python -c "
import pyfaidx
fa = pyfaidx.Fasta('hg38.fa', rebuild=True)
with open('hg38.fa.sizes','w') as f:
    for k in fa.keys(): f.write(f'{k}\t{len(fa[k])}\n')
print('reindexed', len(fa.keys()))
"
rm -f /workspace/genomes/hg38/hg38.fa.ncbi
mkdir -p /root/.local/share/genomes
ln -sfn /workspace/genomes/hg38 /root/.local/share/genomes/hg38

/workspace/venv/bin/python -c "
import genomepy, pandas as pd, grelu.sequence.format
print('installed:', genomepy.list_installed_genomes())
iv = pd.DataFrame({'chrom':['chr1'],'start':[1000000],'end':[1000100],'strand':['+']})
s = grelu.sequence.format.convert_input_type(iv, output_type='strings', genome='hg38')
print('extraction OK, len =', len(s[0]))
"
grep -P '^chr1\t' /workspace/genomes/hg38/hg38.fa.sizes
echo SETUP_COMPLETE
