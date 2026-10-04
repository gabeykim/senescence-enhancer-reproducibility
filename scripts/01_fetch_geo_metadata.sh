#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 01_fetch_geo_metadata.sh
#
# Fetch SOFT + series-matrix metadata for the initial accession set.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/<ACC>/*_family.soft.gz, *_series_matrix.txt.gz, filelist.txt
# ---------------------------------------------------------------------------
# Fetches SOFT / series-matrix / supplementary files for the two accessions under audit
# directly from the NCBI GEO FTP mirror. No SRA/raw-read access. No HTML scraping.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$ROOT/data"

fetch_series () {
  local acc="$1" range="$2"
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${range}/${acc}"
  local outdir="$DATA/$acc"
  mkdir -p "$outdir"
  echo ">> $acc"

  echo "   soft ..."
  curl -sL -o "$outdir/${acc}_family.soft.gz" "$base/soft/${acc}_family.soft.gz" || true
  echo "   matrix ..."
  curl -sL -o "$outdir/${acc}_series_matrix.txt.gz" "$base/matrix/${acc}_series_matrix.txt.gz" || true

  echo "   suppl listing ..."
  # keep only relative filenames (real supplementary files); drop absolute URLs
  # like the boilerplate hhs.gov footer link that also matches href="...".
  curl -sL "$base/suppl/" \
    | grep -oE 'href="[^"/][^"]*"' \
    | sed -E 's/href="([^"]*)"/\1/' \
    | grep -vE '^https?://' \
    > "$outdir/suppl_listing.txt" || true
  while read -r fname; do
    [ -z "$fname" ] && continue
    echo "     -> $fname"
    curl -sL -o "$outdir/$fname" "$base/suppl/$fname"
  done < "$outdir/suppl_listing.txt"

  # sanity: are these real gzip files, or 404 HTML pages saved with a .gz name?
  for f in "$outdir/${acc}_family.soft.gz" "$outdir/${acc}_series_matrix.txt.gz"; do
    if [ -f "$f" ] && ! file "$f" | grep -q gzip; then
      echo "   !! $f is NOT gzip (likely 404/HTML) — accession or path may not exist"
    fi
  done
}

# GSE220545 -- RNA-seq, reported as matched series. Confirmed present on GEO FTP.
fetch_series GSE220545 GSE220nnn

# GSE254358 -- ATAC-seq, senescence. Included for completeness / reproducibility.
# NOTE: as of 2026-08-20 this accession does NOT exist in the GEO database
# (esearch db=gds term=GSE254358[ACCN] -> count=0; FTP path -> HTTP 404;
# GSE254nnn range folder itself resolves fine, so this is not a path-construction bug).
# We still attempt the fetch below so a future re-run auto-detects if it appears.
fetch_series GSE254358 GSE254nnn

echo "Done. Inspect data/<ACC>/ for results; check for 404-HTML masquerading as .gz above."
