#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 21_fetch_gse206402_peaks.py
#
# Fetch GSE206402 narrowPeak supplements.
#
# CONSUMES: GEO FTP (network)
# PRODUCES: data/GSE206402/peaks/
# ---------------------------------------------------------------------------
"""
Downloads the 25 per-GSM narrowPeak.gz files for GSE206402 individually via
GEO's per-sample supplementary mirror (https://ftp.ncbi.nlm.nih.gov/geo/samples/
GSMnnnnnn/<GSM>/suppl/), rather than pulling the single 7.6 GB GSE206402_RAW.tar
that bundles them with the (also-present, not-needed-here) per-sample bigWigs.

Source of the per-GSM filename list: data/GSE206402/filelist.txt, which GEO
publishes alongside the RAW.tar and enumerates every file inside it (name +
byte size) without requiring the tar to be fetched. This script only downloads
the *.narrowPeak.gz entries; *.bw entries are recorded (name + size) but
skipped, consistent with the audit's "note bigWigs, don't download" instruction.

Coordinate system / peak caller are NOT assumed here -- they are inspected
directly from the downloaded file headers and the SOFT record's
Sample_data_processing text in 21_parse_gse206402_metadata.py.
"""
import re
import sys
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "GSE206402"
FILELIST = DATA / "filelist.txt"


def gsm_range_dir(gsm):
    # GEO per-sample suppl mirror uses a GSMnnnnnn "range" folder, same nnn-suffix
    # pattern as the series-level GSExxxnnn folders.
    return re.sub(r"\d{3}$", "nnn", gsm)


def main():
    if not FILELIST.exists():
        sys.exit(f"missing {FILELIST} -- run scripts/20_fetch_gse206402_family.sh first")

    entries = []
    with open(FILELIST, encoding="utf-8") as fh:
        next(fh)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] != "File":
                continue
            fname, size = parts[1], parts[3]
            entries.append((fname, int(size)))

    narrowpeak = [(f, s) for f, s in entries if f.endswith(".narrowPeak.gz")]
    bigwig = [(f, s) for f, s in entries if f.endswith(".bw")]
    print(f"filelist.txt: {len(entries)} per-sample files total "
          f"({len(narrowpeak)} narrowPeak.gz, {len(bigwig)} bigWig)", file=sys.stderr)

    print("\n-- bigWig files present but NOT downloaded (name, size bytes) --", file=sys.stderr)
    for fname, size in bigwig:
        print(f"  {fname}\t{size}", file=sys.stderr)

    ok, failed = 0, []
    for fname, size in narrowpeak:
        gsm_m = re.match(r"(GSM\d+)_", fname)
        if not gsm_m:
            failed.append((fname, "no GSM prefix parsed"))
            continue
        gsm = gsm_m.group(1)
        url = (f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm_range_dir(gsm)}/"
               f"{gsm}/suppl/{fname}")
        dest = DATA / "peaks" / fname
        dest.parent.mkdir(exist_ok=True)
        r = subprocess.run(["curl", "-sL", "-o", str(dest), url])
        if r.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
            failed.append((fname, "download failed/empty"))
            continue
        # sanity: expected size vs actual (gzip is deterministic-ish but NCBI
        # timestamps can differ by re-gzip; just flag gross mismatches)
        actual = dest.stat().st_size
        if abs(actual - size) > 0.05 * size:
            print(f"  !! size mismatch for {fname}: filelist.txt says {size}, got {actual}",
                  file=sys.stderr)
        ok += 1
        print(f"  OK  {fname} ({actual} bytes) -> {dest}", file=sys.stderr)

    print(f"\nDownloaded {ok}/{len(narrowpeak)} narrowPeak.gz files to {DATA/'peaks'}", file=sys.stderr)
    if failed:
        print("FAILED:", failed, file=sys.stderr)


if __name__ == "__main__":
    main()
