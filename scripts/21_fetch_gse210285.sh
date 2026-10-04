#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 21_fetch_gse210285.sh
#
# Fetch GSE210285 (multi-omics senescence) metadata.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE210285/
# ---------------------------------------------------------------------------
# Fetches SOFT / series-matrix / supplementary listing for GSE210285 directly
# from the NCBI GEO FTP mirror. No SRA/raw-read access. No HTML scraping of
# GEO's web UI -- only the plain-text SOFT/matrix files and a directory
# listing (parsed for filenames only) are touched.
#
# GSE210285 has TWO platforms (GPL20301 Illumina HiSeq 4000, GPL23227
# BGISEQ-500), so GEO publishes two per-platform series-matrix files
# instead of one -- both are fetched below.
#
# Supplementary data on this series is bigWig signal tracks only (4 files,
# one per GSM) bundled into one GSE210285_RAW.tar (~221 MiB / 232,181,760
# bytes per filelist.txt). Per task instructions, bigWigs are NOT downloaded
# (noted as existence + size only) -- filelist.txt (334 bytes) is fetched in
# full since it is small and gives exact per-file sizes/types without
# pulling the tar itself.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/data/GSE210285"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210285"
mkdir -p "$DATA"

echo ">> GSE210285"

echo "   family soft ..."
curl -sL -o "$DATA/GSE210285_family.soft.gz" "$BASE/soft/GSE210285_family.soft.gz"
gunzip -kf "$DATA/GSE210285_family.soft.gz"

echo "   series matrix (per-platform: GPL20301, GPL23227) ..."
curl -sL -o "$DATA/GSE210285-GPL20301_series_matrix.txt.gz" \
  "$BASE/matrix/GSE210285-GPL20301_series_matrix.txt.gz"
curl -sL -o "$DATA/GSE210285-GPL23227_series_matrix.txt.gz" \
  "$BASE/matrix/GSE210285-GPL23227_series_matrix.txt.gz"
gunzip -kf "$DATA/GSE210285-GPL20301_series_matrix.txt.gz"
gunzip -kf "$DATA/GSE210285-GPL23227_series_matrix.txt.gz"

echo "   suppl directory listing (filenames only, no scraping of rendered HTML content) ..."
curl -sL "$BASE/suppl/" \
  | grep -oE 'href="[^"/][^"]*"' \
  | sed -E 's/href="([^"]*)"/\1/' \
  | grep -vE '^https?://' \
  > "$DATA/suppl_listing.txt"
cat "$DATA/suppl_listing.txt"

echo "   suppl/filelist.txt (small, exact per-file sizes/types) ..."
curl -sL -o "$DATA/filelist.txt" "$BASE/suppl/filelist.txt"
cat "$DATA/filelist.txt"

echo "   NOTE: GSE210285_RAW.tar (~221 MiB) and the 4 constituent .bw files"
echo "   are intentionally NOT downloaded -- bigWig signal tracks, not peak"
echo "   calls, and flagged as 'very large' under task constraints. Existence"
echo "   and exact sizes are captured above via filelist.txt."

# sanity: are these real gzip files, or 404 HTML pages saved with a .gz name?
for f in "$DATA/GSE210285_family.soft.gz" \
         "$DATA/GSE210285-GPL20301_series_matrix.txt.gz" \
         "$DATA/GSE210285-GPL23227_series_matrix.txt.gz"; do
  if [ -f "$f" ] && ! file "$f" | grep -q gzip; then
    echo "   !! $f is NOT gzip (likely 404/HTML) -- accession or path may not exist"
  fi
done

echo "Done. Inspect data/GSE210285/ for results."
