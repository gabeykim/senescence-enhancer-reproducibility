#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 60_fetch_replicative_ceiling.sh
#
# Fetch the replicative/IR/OIS bigWigs for the ceiling test.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE146585/bw/, data/GSE106146/bw/
# ---------------------------------------------------------------------------
# Fetch H3K27ac bigWigs for the within-replicative ceiling test.
# Per-sample/series files only -- no RAW.tar archives.
set -uo pipefail
ROOT="/Users/gabeykim/Downloads/Senescence"; cd "$ROOT"

gsm_get () {   # $1 GSM, $2 filename, $3 destdir
  local pre="${1:0:$(( ${#1} - 3 ))}nnn" dest="$3/$2"
  mkdir -p "$3"; [ -s "$dest" ] && { echo "  cached  $2"; return 0; }
  curl -sS -f --retry 3 --max-time 3600 -o "$dest" \
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/${pre}/${1}/suppl/${2}" \
    && echo "  ok      $2 ($(du -h "$dest"|cut -f1))" || { echo "  FAILED $2"; rm -f "$dest"; }
}
ser_get () {   # $1 series, $2 filename, $3 destdir
  local pre; pre=$(echo "$1" | sed -E 's/[0-9]{1,3}$/nnn/')
  local dest="$3/$2"
  mkdir -p "$3"; [ -s "$dest" ] && { echo "  cached  $2"; return 0; }
  curl -sS -f --retry 3 --max-time 3600 -o "$dest" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${pre}/${1}/suppl/${2}" \
    && echo "  ok      $2 ($(du -h "$dest"|cut -f1))" || { echo "  FAILED $2"; rm -f "$dest"; }
}

S="_STARdef_unique_tbp1.ucsc.bigWig"
echo "=== GSE146585 / GSE146567  IMR90 REPLICATIVE H3K27ac (n=2 per arm, hg18) ==="
gsm_get GSM4395479 "GSM4395479_TSD311_IMR90_Young_anti-H3K27ac_Rep1${S}"      data/GSE146585/bw
gsm_get GSM4395480 "GSM4395480_TSD312_IMR90_Young_anti-H3K27ac_Rep2${S}"      data/GSE146585/bw
gsm_get GSM4395481 "GSM4395481_TSD313_IMR90_Senescent_anti-H3K27ac_Rep1${S}"  data/GSE146585/bw
gsm_get GSM4395482 "GSM4395482_TSD314_IMR90_Senescent_anti-H3K27ac_Rep2${S}"  data/GSE146585/bw

echo "=== GSE146585 / GSE146559  BJ REPLICATIVE H3K27ac (n=1 per arm, hg18) ==="
gsm_get GSM4395086 "GSM4395086_Tom161_BJ_Young_anti-H3K27ac${S}"              data/GSE146585/bw
gsm_get GSM4395091 "GSM4395091_Tom164_BJ_Senescent_anti-H3K27ac${S}"          data/GSE146585/bw

echo "=== GSE146585 / GSE146563  IMR90 OIS H3K27ac (n=1 per arm, hg18) -- same lab/pipeline as above ==="
gsm_get GSM4395451 "GSM4395451_TSC801_IMR90_GFP_Control_anti-H3K27ac${S}"     data/GSE146585/bw
gsm_get GSM4395452 "GSM4395452_TSC802_IMR90_RAS_OIS_anti-H3K27ac${S}"         data/GSE146585/bw

echo "=== GSE106146  fibroblast REPLICATIVE + IR H3K27ac (hg19) ==="
ser_get GSE106146 "GSE106146_P2-H3K27ac.bigWig"   data/GSE106146/bw   # proliferating (n=1)
ser_get GSE106146 "GSE106146_S1-H3K27ac.bigWig"   data/GSE106146/bw   # senescent    (n=1)
ser_get GSE106146 "GSM2830447_IR-H3K27ac.Rep1.bigWig" data/GSE106146/bw
ser_get GSE106146 "GSM2830448_IR-H3K27ac.Rep2.bigWig" data/GSE106146/bw

echo; du -sh data/GSE146585/bw data/GSE106146/bw 2>/dev/null; echo FETCH_COMPLETE
