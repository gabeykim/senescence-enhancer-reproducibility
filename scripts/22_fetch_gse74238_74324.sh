#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 22_fetch_gse74238_74324.sh
#
# Fetch GSE74238 (OIS H3K27ac) and GSE74324 (OIS RNA) metadata.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE74238/, data/GSE74324/
# ---------------------------------------------------------------------------
# Fetches SOFT family file + series matrix + supplementary-file listing/data for
# GSE74238 (ChIP-seq: H3K27ac/BRD4, IMR90 OIS) and GSE74324 (matched RNA-seq,
# 72 samples) -- Tasdemir et al., Cancer Discovery 2016, doi:10.1158/2159-8290.CD-16-0217.
#
# No SRA/FASTQ access. No HTML scraping for metadata (SOFT/matrix text only --
# the suppl/ directory LISTING is scraped for filenames+sizes, same pattern as
# scripts/01_fetch_geo_metadata.sh and scripts/20_fetch_gse206402_family.sh).
#
# GSE74238 supplementary data is one big GSE74238_RAW.tar (3.79 GB, bundling 14
# per-sample bigWig ChIP-signal tracks, 25 MB-780 MB each) plus filelist.txt.
# Per task constraints (no large bigWig/bam downloads), we download filelist.txt
# (small, gives per-file names+sizes+types) but do NOT download the tar or any
# bigWig. GSE74324 has one processed matrix, GSE74324_ALL_samples_rpkm.txt.gz
# (~8.6 MB gzipped) -- small enough to download in full. GSE74324 has no
# filelist.txt on FTP (only GSE74238 does) -- fetch_suppl() below detects and
# discards the 404 HTML page that request returns, rather than saving it as if
# it were real data.
set -uo pipefail  # no -e: keep going if one accession/file has a hiccup

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/data"

is_gzip () { file "$1" 2>/dev/null | grep -q gzip; }

fetch_soft_and_matrix () {
  local acc="$1" range="$2"
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${range}/${acc}"
  local outdir="$DATA/$acc"
  mkdir -p "$outdir"
  echo ">> $acc  (soft + matrix)"

  echo "   soft ..."
  curl -sL -o "$outdir/${acc}_family.soft.gz" "$base/soft/${acc}_family.soft.gz" || true
  echo "   matrix ..."
  curl -sL -o "$outdir/${acc}_series_matrix.txt.gz" "$base/matrix/${acc}_series_matrix.txt.gz" || true

  for f in "$outdir/${acc}_family.soft.gz" "$outdir/${acc}_series_matrix.txt.gz"; do
    if [ -f "$f" ] && ! is_gzip "$f"; then
      echo "   !! $f is NOT gzip (likely 404/HTML) -- accession or path may not exist"
    fi
  done

  # Decompress alongside the .gz (keep both) -- scripts/23_parse_gse74238_74324.py reads
  # the plain-text .soft file directly, same convention as scripts/02_parse_metadata.py.
  gunzip -k -f "$outdir/${acc}_family.soft.gz" 2>/dev/null || true
  gunzip -k -f "$outdir/${acc}_series_matrix.txt.gz" 2>/dev/null || true
}

fetch_suppl () {
  # Directory listing (name + size + type, from GEO's own filelist.txt where
  # present -- more reliable than scraping the Apache HTML index for sizes).
  # Small files (< ~20 MB, e.g. processed matrices/DE tables) are downloaded in
  # full; bigWig/bam/tar archives of raw tracks are recorded (name + size) but
  # NOT downloaded, per task constraints.
  local acc="$1" range="$2"
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${range}/${acc}"
  local outdir="$DATA/$acc"
  mkdir -p "$outdir"
  echo ">> $acc  (suppl listing + selective download)"

  # GEO's own machine-readable per-file manifest (name, timestamp, size in bytes, type),
  # where one exists. Not every series has one (GSE74324 doesn't) -- a missing manifest
  # 404s to an HTML error page, which we detect and discard rather than mistake for data.
  curl -sL -o "$outdir/filelist.txt" "$base/suppl/filelist.txt" || true
  if [ -f "$outdir/filelist.txt" ] && grep -qi "<html" "$outdir/filelist.txt"; then
    echo "   (no filelist.txt for $acc -- 404 HTML page discarded)"
    rm -f "$outdir/filelist.txt"
  fi

  # href-scrape of the suppl/ index purely to get top-level archive/file NAMES
  # (filelist.txt already gives sizes/types for files inside GSE74238_RAW.tar,
  # but the top-level suppl/ listing is what actually exists to fetch).
  curl -sL "$base/suppl/" \
    | grep -oE 'href="[^"/][^"]*"' \
    | sed -E 's/href="([^"]*)"/\1/' \
    | grep -vE '^https?://' \
    > "$outdir/suppl_listing.txt" || true

  while read -r fname; do
    [ -z "$fname" ] && continue
    [ "$fname" = "filelist.txt" ] && continue  # already fetched above
    lower=$(echo "$fname" | tr '[:upper:]' '[:lower:]')
    case "$lower" in
      *.bw|*.bigwig|*.bam|*.bai|*raw*.tar|*.tar)
        size=$(curl -sIL "$base/suppl/$fname" | grep -i '^content-length' | tr -d '\r' | awk '{print $2}')
        echo "   SKIP (large/raw track, not downloaded): $fname (Content-Length: ${size:-unknown} bytes)"
        continue
        ;;
    esac
    echo "   -> $fname"
    curl -sL -o "$outdir/$fname" "$base/suppl/$fname" || true
  done < "$outdir/suppl_listing.txt"
}

# ---- GSE74238: ChIP-seq (H3K27ac, BRD4; proliferating/quiescent/senescent + a
# DMSO/Etoposide pair) ----
fetch_soft_and_matrix GSE74238 GSE74nnn
fetch_suppl GSE74238 GSE74nnn

# ---- GSE74324: matched RNA-seq, 72 samples ----
fetch_soft_and_matrix GSE74324 GSE74nnn
fetch_suppl GSE74324 GSE74nnn

echo "Done. Inspect data/GSE74238/, data/GSE74324/."
echo "NOTE: GSE74238_RAW.tar (3.79 GB; 14 bigWig ChIP-signal tracks, 25 MB-780 MB"
echo "each per filelist.txt) was deliberately NOT downloaded -- see filelist.txt"
echo "for per-track names/sizes/types."
