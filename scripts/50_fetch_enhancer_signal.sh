#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 50_fetch_enhancer_signal.sh
#
# Fetch bigWigs and peaks for the enhancer cross-mechanism test.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE205898/, GSE206402/, GSE74238/ bw+peaks
# ---------------------------------------------------------------------------
# Fetch ONLY the per-sample signal/peak files needed for the cross-mechanism enhancer test.
#
# GEO serves per-sample supplementary files under /geo/samples/GSM<prefix>nnn/<GSM>/suppl/,
# so the multi-GB RAW.tar archives (GSE205898 9.0 GB, GSE206402 8.2 GB, GSE74238 3.8 GB)
# are NOT downloaded -- only the individual bigWigs/narrowPeaks actually used.
#
# curl, not python urllib: this Mac's python.org build has no configured CA bundle.
set -uo pipefail
ROOT="/Users/gabeykim/Downloads/Senescence"
cd "$ROOT"

gsm_url () {           # $1 = GSM id, $2 = filename
  local pre="${1:0:$(( ${#1} - 3 ))}nnn"
  echo "https://ftp.ncbi.nlm.nih.gov/geo/samples/${pre}/${1}/suppl/${2}"
}

get () {               # $1 = GSM, $2 = filename, $3 = dest dir
  local dest="$3/$2"
  mkdir -p "$3"
  if [ -s "$dest" ]; then echo "  cached  $2"; return 0; fi
  local url; url="$(gsm_url "$1" "$2")"
  if curl -sS -f --retry 3 --max-time 3600 -o "$dest" "$url"; then
    echo "  ok      $2  ($(du -h "$dest" | cut -f1))"
  else
    echo "  FAILED  $2   <- $url"; rm -f "$dest"; return 1
  fi
}

echo "=== GSE74238  IMR90 OIS H3K27ac bigWigs (no peaks were ever deposited) ==="
get GSM1915113 "GSM1915113_H3K27ac_proliferating.bigWig"           data/GSE74238/bw
get GSM1915115 "GSM1915115_H3K27ac_senescent..bigWig"              data/GSE74238/bw
get GSM2098176 "GSM2098176_H3K27ac_proliferating_replicate.bigWig" data/GSE74238/bw
get GSM2098178 "GSM2098178_H3K27ac_senescent_replicate.bigWig"     data/GSE74238/bw

echo "=== GSE205898  GM21 OIS H3K27ac bigWigs + narrowPeaks ==="
get GSM6235098 "GSM6235098_EV_K27ac_1_dedup_blacklisted_merge.bw"       data/GSE205898/bw
get GSM6235099 "GSM6235099_EV_K27ac_2_dedup_blacklisted_merge.bw"       data/GSE205898/bw
get GSM6235110 "GSM6235110_RAS_D18_K27ac_1_dedup_blacklisted_merge.bw"  data/GSE205898/bw
get GSM6235111 "GSM6235111_RAS_D18_K27ac_2_dedup_blacklisted_merge.bw"  data/GSE205898/bw
get GSM6235098 "GSM6235098_EV_K27ac_Rep1_peaks.narrowPeak.gz"           data/GSE205898/peaks
get GSM6235099 "GSM6235099_EV_K27ac_Rep2_peaks.narrowPeak.gz"           data/GSE205898/peaks
get GSM6235110 "GSM6235110_RAS_D18_K27ac_Rep1_peaks.narrowPeak.gz"      data/GSE205898/peaks
get GSM6235111 "GSM6235111_RAS_D18_K27ac_Rep2_peaks.narrowPeak.gz"      data/GSE205898/peaks

echo "=== GSE206402  GM21 OIS ATAC bigWigs (narrowPeaks already local) ==="
get GSM6253126 "GSM6253126_pBABE_D8_REP1.bw"   data/GSE206402/bw
get GSM6253127 "GSM6253127_pBABE_D8_REP2.bw"   data/GSE206402/bw
get GSM6253122 "GSM6253122_pBABE_D32_REP1.bw"  data/GSE206402/bw
get GSM6253123 "GSM6253123_pBABE_D32_REP2.bw"  data/GSE206402/bw
get GSM6253132 "GSM6253132_RAS_D18_REP1.bw"    data/GSE206402/bw
get GSM6253133 "GSM6253133_RAS_D18_REP2.bw"    data/GSE206402/bw
get GSM6253134 "GSM6253134_RAS_D23_REP1.bw"    data/GSE206402/bw
get GSM6253135 "GSM6253135_RAS_D23_REP2.bw"    data/GSE206402/bw

echo
echo "TOTAL DOWNLOADED:"
du -sh data/GSE74238/bw data/GSE205898/bw data/GSE206402/bw 2>/dev/null
echo FETCH_COMPLETE
