#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 10_fetch_gse175533.sh
#
# Fetch GSE175533 (replicative senescence, Hayflick) metadata + supplements.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE175533/
# ---------------------------------------------------------------------------
# Fetch GSE175533 (Calico Hayflick-limit multi-omics) metadata + the processed
# supplementary files needed for this audit. Deliberately SKIPS:
#   - GSE175533_RAW.tar (14.8 GB of per-sample bigWigs -- not needed for the
#     peak-atlas/TPM-table QC this audit runs, and far too large to be a
#     reasonable download for a metadata+processed-file audit)
#   - GSE175533_sceasy_hay.h5ad.gz (single-cell RNA-seq object -- out of scope;
#     the task's RNA arm is the bulk TPM table)
# Both are noted, not silently dropped -- see the fetch log this script prints.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/data/GSE175533"
mkdir -p "$OUT"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175533"

echo ">> SOFT + matrix"
curl -sL -o "$OUT/GSE175533_family.soft.gz" "$BASE/soft/GSE175533_family.soft.gz"
curl -sL -o "$OUT/GSE175533_series_matrix.txt.gz" "$BASE/matrix/GSE175533_series_matrix.txt.gz"
gunzip -kf "$OUT/GSE175533_family.soft.gz" "$OUT/GSE175533_series_matrix.txt.gz"

echo ">> Supplementary file listing"
curl -sL "$BASE/suppl/filelist.txt" -o "$OUT/filelist.txt"
cat "$OUT/filelist.txt"

echo ">> Downloading processed files needed for this audit (peak atlas, sig-peak table, TPM table)"
for f in GSE175533_atlas.bed.gz GSE175533_atac_peaks_sig.xlsx GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx; do
  echo "   -> $f"
  curl -sL -o "$OUT/$f" "$BASE/suppl/$f"
done
gunzip -kf "$OUT/GSE175533_atlas.bed.gz"

echo ">> Skipped (documented, not downloaded): GSE175533_RAW.tar (14.8GB bigWigs), GSE175533_sceasy_hay.h5ad.gz (scRNA object, out of scope)"
ls -la "$OUT"
