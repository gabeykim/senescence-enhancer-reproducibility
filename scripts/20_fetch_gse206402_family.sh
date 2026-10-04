#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 20_fetch_gse206402_family.sh
#
# Fetch GSE206402 (OIS ATAC) metadata.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE206402/
# ---------------------------------------------------------------------------
# Fetches metadata (+ supplementary files, where reasonable) for the GSE206402
# SuperSeries family (Martinez-Zamudio et al., Cell Genomics 2023,
# doi:10.1016/j.xgen.2023.100293 -- escape from OIS in GM21 fibroblasts via ER:RAS).
#
# Primary focus: GSE206402 (ATAC-seq subseries) -- full SOFT family file, series
# matrix, AND every supplementary file listed on the FTP suppl/ directory (peak
# calls / count matrices are downloaded in full; anything that looks like a
# bigWig/bam-scale file is recorded by name+size from the directory listing but
# NOT downloaded, per task constraints -- no raw/large alignment tracks).
#
# Secondary (context only): GSE206496 (parent SuperSeries), GSE205898
# (H3K27ac/H3K4me1 ChIP-seq subseries), GSE206493 (RNA-seq subseries) -- series
# matrix only (lighter than the full family.soft, sufficient to report design).
#
# No SRA/FASTQ access. No HTML scraping for metadata (only the suppl/ directory
# LISTING is scraped, same pattern as scripts/01_fetch_geo_metadata.sh, to
# enumerate filenames -- the metadata itself always comes from SOFT/matrix text).
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
}

fetch_matrix_only () {
  # lighter fetch for the SuperSeries family members we only need design context from
  local acc="$1" range="$2"
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${range}/${acc}"
  local outdir="$DATA/$acc"
  mkdir -p "$outdir"
  echo ">> $acc  (matrix only)"
  curl -sL -o "$outdir/${acc}_series_matrix.txt.gz" "$base/matrix/${acc}_series_matrix.txt.gz" || true
  local f="$outdir/${acc}_series_matrix.txt.gz"
  if [ -f "$f" ] && ! is_gzip "$f"; then
    echo "   !! $f is NOT gzip (likely 404/HTML) -- accession or path may not exist"
  fi
}

fetch_multiplatform_matrix () {
  # GSE206496 is a SuperSeries spanning 3 platforms -- GEO does not publish a
  # single flat <acc>_series_matrix.txt.gz for it (confirmed 404 as of
  # 2026-08-20); instead one matrix per platform exists under matrix/. We only
  # need the two small ones (ATAC+ChIP GSM range, and RNA-seq GSM range) for
  # design context here; the microarray platform's matrix (7.4M, GSE205692
  # subseries, out of this audit's ATAC-seq scope) is deliberately skipped --
  # see data/GSE206496/NOTE_matrix_structure.txt.
  local acc="$1" range="$2"
  shift 2
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${range}/${acc}"
  local outdir="$DATA/$acc"
  mkdir -p "$outdir"
  echo ">> $acc  (per-platform matrix, selected platforms only)"
  for plat in "$@"; do
    local fname="${acc}-${plat}_series_matrix.txt.gz"
    curl -sL -o "$outdir/$fname" "$base/matrix/$fname" || true
    if [ -f "$outdir/$fname" ] && ! is_gzip "$outdir/$fname"; then
      echo "   !! $outdir/$fname is NOT gzip (likely 404/HTML)"
    else
      echo "   -> $fname"
    fi
  done
}

fetch_suppl () {
  # Full suppl/ directory listing (name + size, from the FTP HTML index) for
  # GSE206402 specifically. Small files (peaks/counts/DE tables) are downloaded
  # in full; files that look bigWig/bam/large-alignment-track scale are recorded
  # (name + size) but NOT downloaded -- per task constraints.
  local acc="$1" range="$2"
  local base="https://ftp.ncbi.nlm.nih.gov/geo/series/${range}/${acc}"
  local outdir="$DATA/$acc"
  mkdir -p "$outdir"
  echo ">> $acc  (suppl listing + selective download)"

  # Raw directory index (Apache-style listing: name + date + size columns).
  curl -sL "$base/suppl/" > "$outdir/suppl_index.html" || true

  # filenames only (same href-scrape pattern as 01_fetch_geo_metadata.sh)
  grep -oE 'href="[^"/][^"]*"' "$outdir/suppl_index.html" \
    | sed -E 's/href="([^"]*)"/\1/' \
    | grep -vE '^https?://' \
    > "$outdir/suppl_listing.txt" || true

  # name + size table, parsed straight out of the Apache listing text (last two
  # whitespace-separated tokens on the line containing the href are date/time
  # and size in this listing style -- verified by eyeballing suppl_index.html).
  python3 - "$outdir/suppl_index.html" "$outdir/suppl_listing_with_sizes.txt" <<'PYEOF'
import re, sys
html_path, out_path = sys.argv[1], sys.argv[2]
text = open(html_path, encoding="utf-8", errors="replace").read()
rows = []
for line in text.splitlines():
    m = re.search(r'href="([^"/][^"]*)"', line)
    if not m:
        continue
    fname = m.group(1)
    if fname.startswith("http"):
        continue
    # trailing size token, e.g. "...</a>   2023-04-05 12:34  15M" or similar
    size_m = re.search(r'(\d+(?:\.\d+)?[KMG]?)\s*$', line.strip())
    size = size_m.group(1) if size_m else "?"
    rows.append((fname, size))
with open(out_path, "w") as fh:
    for fname, size in rows:
        fh.write(f"{fname}\t{size}\n")
print(f"wrote {out_path} ({len(rows)} entries)", file=sys.stderr)
PYEOF

  while read -r fname; do
    [ -z "$fname" ] && continue
    lower=$(echo "$fname" | tr '[:upper:]' '[:lower:]')
    case "$lower" in
      *.bw|*.bigwig|*.bam|*.bai|*raw*.tar)
        echo "   SKIP (large track, not downloaded): $fname"
        continue
        ;;
    esac
    echo "   -> $fname"
    curl -sL -o "$outdir/$fname" "$base/suppl/$fname" || true
  done < "$outdir/suppl_listing.txt"
}

# ---- Primary: GSE206402 (ATAC-seq subseries) ----
fetch_soft_and_matrix GSE206402 GSE206nnn
fetch_suppl GSE206402 GSE206nnn

# ---- Context only: SuperSeries + other assay subseries ----
# GSE206496 (parent SuperSeries) has NO flat series_matrix -- it spans 3
# platforms (GPL16791 ATAC+ChIP, GPL17586 microarray, GPL24676 RNA-seq).
# Only the two small ones are pulled here; GPL17586 (7.4M, microarray
# subseries GSE205692 -- discovered during this audit, not in the original
# task list) is out of scope for an ATAC-seq-focused audit.
fetch_multiplatform_matrix GSE206496 GSE206nnn GPL16791 GPL24676
fetch_matrix_only GSE205898 GSE205nnn   # ChIP-seq (H3K27ac/H3K4me1) subseries
fetch_matrix_only GSE206493 GSE206nnn   # RNA-seq subseries

echo "Done. Inspect data/GSE206402/, data/GSE206496/, data/GSE205898/, data/GSE206493/."
