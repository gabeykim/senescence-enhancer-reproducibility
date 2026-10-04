#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 40_fetch_expansion_series.sh
#
# Fetch the additional series for the multi-study expansion.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE205692/, GSE220545/, GSE99028/
# ---------------------------------------------------------------------------
# Fetch GEO metadata (SOFT + series matrix) and enumerate supplementary files for the
# expansion series. Processed supplementary files only -- no SRA, no HTML metadata scraping
# (the suppl/ directory listing is scraped ONLY to enumerate filenames; all metadata is
# parsed from SOFT/series-matrix text).
#
# NOTE: this Mac's python has no CA bundle (SSLCertVerificationError from urllib), so every
# network call in this project shells out to curl. Do not "fix" this by switching to urllib.
set -uo pipefail
ROOT="/Users/gabeykim/Downloads/Senescence"
cd "$ROOT"

fetch_series () {
  local acc="$1"
  local pre
  pre=$(echo "$acc" | sed -E 's/[0-9]{1,3}$/nnn/')
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${pre}/${acc}"
  local dir="data/${acc}"
  mkdir -p "$dir"
  echo "=== ${acc} ==="

  # family SOFT
  if [ ! -s "${dir}/${acc}_family.soft" ]; then
    curl -sS -f -o "${dir}/${acc}_family.soft.gz" "${base}/soft/${acc}_family.soft.gz" \
      && gunzip -kf "${dir}/${acc}_family.soft.gz" \
      && echo "  soft OK ($(wc -l < "${dir}/${acc}_family.soft") lines)" \
      || echo "  soft MISSING"
  else
    echo "  soft cached ($(wc -l < "${dir}/${acc}_family.soft") lines)"
  fi

  # series matrix (may be split per-platform)
  if ! ls "${dir}"/*series_matrix.txt >/dev/null 2>&1; then
    curl -s "${base}/matrix/" \
      | grep -oE "${acc}[A-Za-z0-9_.-]*series_matrix\.txt\.gz" | sort -u \
      | while read -r m; do
          curl -sS -f -o "${dir}/${m}" "${base}/matrix/${m}" && gunzip -kf "${dir}/${m}"
        done
    echo "  matrix: $(ls "${dir}"/*series_matrix.txt 2>/dev/null | wc -l) file(s)"
  else
    echo "  matrix cached: $(ls "${dir}"/*series_matrix.txt 2>/dev/null | wc -l) file(s)"
  fi

  # supplementary listing. Filter out absolute URLs -- a prior bug in this project scraped an
  # hhs.gov footer link as if it were a filename, and the resulting curl -o failure (exit 23)
  # killed the whole fetch script before it reached the next accession.
  curl -s "${base}/suppl/" \
    | grep -oE 'href="[^"]+"' | sed 's/href="//; s/"$//' \
    | grep -vE '^https?://' | grep -vE '^/|^\?|^$' \
    | sort -u > "${dir}/suppl_listing.txt"
  echo "  suppl files:"
  sed 's/^/    /' "${dir}/suppl_listing.txt"

  # GEO's own manifest, when present (authoritative sizes; 404s return an HTML page, so
  # verify it is really TSV before keeping it)
  curl -sS -f -o "${dir}/filelist.txt" "${base}/suppl/filelist.txt" 2>/dev/null
  if [ -s "${dir}/filelist.txt" ] && head -1 "${dir}/filelist.txt" | grep -qi 'html'; then
    rm -f "${dir}/filelist.txt"
  fi
  echo
}

for acc in GSE99028 GSE206493 GSE205692 GSE210285; do
  fetch_series "$acc"
done
echo "FETCH_COMPLETE"
