#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 22_parse_gse206402_metadata.py
#
# Parse GSE206402 metadata and peak inventory.
#
# CONSUMES: data/GSE206402/*_family.soft.gz, peaks/*.narrowPeak.gz
# PRODUCES: output/tables/gse206402_sample_metadata.csv
# ---------------------------------------------------------------------------
"""
Parse GSE206402_family.soft (ATAC-seq subseries) into a tidy one-row-per-GSM
CSV. Direct SOFT parsing (no GEOparse dependency), same pattern as
scripts/02_parse_metadata.py.

Output: output/tables/gse206402_sample_metadata.csv, columns:
  gsm, title, cell_line, condition, timepoint, treatment, replicate_id,
  raw_characteristics
(plus a few extra audit columns: source_name, genotype/library info,
supplementary_file names, and a peak-file-internal-label cross-check flag
-- see PEAK_LABEL_MISMATCH note below.)

ASSUMPTION / PARSING NOTE (flagged, not silently resolved): GEO's
Sample_characteristics_ch1 for this series carries exactly three keys --
"cell line", "cell type", "treatment" -- and NO explicit per-sample
"condition" or "senescence/escape status" field. "condition" in this CSV is
therefore derived from the "treatment" characteristic (Empty vector control
vs H-RAS-G12-V overexpression) plus the parsed timepoint, NOT from any
GEO-native status vocabulary. In particular: GEO gives NO machine-readable
label distinguishing "still senescent" from "escaped" within the RAS arm's
7 timepoints (D8/14/18/23/32/45/56) -- that distinction, if it exists, lives
only in the paper's Results text (not fetched here; out of scope per the
GEO-metadata-only constraint), so `condition` for RAS-arm rows is left as
the coarse label "RAS overexpression (OIS induction/time-course; escape
sub-status NOT in GEO metadata)" rather than guessed.
"""
import re
import sys
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "GSE206402"
OUT = ROOT / "output" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

SOFT = DATA / "GSE206402_family.soft"
SOFT_GZ = DATA / "GSE206402_family.soft.gz"


def open_text(path_plain, path_gz):
    if path_plain.exists():
        return open(path_plain, encoding="utf-8", errors="replace")
    if path_gz.exists():
        return gzip.open(path_gz, mode="rt", encoding="utf-8", errors="replace")
    sys.exit(f"neither {path_plain} nor {path_gz} exists -- run scripts/20_fetch_gse206402_family.sh first")


def parse_soft(fh):
    """Return (series_dict, [sample_dict, ...]) -- same shape as 02_parse_metadata.py."""
    samples = []
    series = {}
    cur = None
    cur_type = None
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            if cur is not None and cur_type == "SAMPLE":
                samples.append(cur)
            cur = {"characteristics": []}
            cur_type = "SAMPLE"
            continue
        if line.startswith("^SERIES"):
            cur = series
            cur_type = "SERIES"
            continue
        if line.startswith("^PLATFORM"):
            if cur is not None and cur_type == "SAMPLE":
                samples.append(cur)
            cur = None
            cur_type = "PLATFORM"
            continue
        if not line.startswith("!") or cur is None:
            continue
        m = re.match(r"!(\w+)\s*=\s*(.*)", line)
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        if key == "Sample_characteristics_ch1":
            cur["characteristics"].append(val)
        elif key in cur:
            cur.setdefault(key + "_extra", []).append(val)
        else:
            cur[key] = val
    if cur is not None and cur_type == "SAMPLE":
        samples.append(cur)
    return series, samples


def chars_to_dict(char_list):
    d = {}
    for c in char_list:
        if ":" in c:
            k, v = c.split(":", 1)
            d[k.strip().lower()] = v.strip()
        else:
            d[c.strip().lower()] = ""
    return d


def infer_timepoint(title, chars_text):
    # "day 32 post-transduction" / "day 8 post-transduction" etc.
    m = re.search(r"day\s*(\d+)\s*post-transduction", title, re.I)
    if m:
        return f"day {m.group(1)}"
    m = re.search(r"day\s*(\d+)", chars_text, re.I)
    if m:
        return f"day {m.group(1)}"
    return "UNRESOLVED -- no day/timepoint pattern matched in title or characteristics"


def infer_replicate(title):
    m = re.search(r"rep\s*(\d+)", title, re.I)
    return m.group(1) if m else "UNRESOLVED -- no 'repN' pattern in title"


def infer_condition(treatment_raw):
    t = treatment_raw.lower()
    if "empty vector" in t:
        return "empty-vector control (proliferating/non-RAS baseline at matched timepoint)"
    if "ras" in t and "overexpression" in t:
        return ("RAS overexpression (OIS induction/time-course; escape sub-status "
                "NOT resolvable from GEO metadata -- see script docstring)")
    return f"UNRESOLVED -- unrecognized treatment string: {treatment_raw!r}"


def load_peak_label_mismatches():
    """
    Cross-check: the narrowPeak.gz files downloaded by 21_fetch_gse206402_peaks.py
    carry an internal peak-name prefix (column 4, e.g. 'RAS_D18_1') that in
    several cases DISAGREES with the GEO filename/title day or replicate number
    (e.g. GSM6253135 filename says RAS_D23_Rep2, internal peak names say
    RAS_D22_1). This does not change what's written to the metadata CSV (GEO's
    own title/characteristics are authoritative for the CSV), but it is
    recorded per-GSM here so the flag travels with the row instead of being
    buried in a separate log only.
    """
    peaks_dir = DATA / "peaks"
    mism = {}
    if not peaks_dir.exists():
        return mism
    import gzip as _gzip
    for f in sorted(peaks_dir.glob("*.narrowPeak.gz")):
        gsm_m = re.match(r"(GSM\d+)_", f.name)
        if not gsm_m:
            continue
        gsm = gsm_m.group(1)
        # expected label derived from the filename itself, e.g.
        # GSM6253135_RAS_D23_Rep2_peaks.narrowPeak.gz -> RAS_D23_Rep2
        fname_label_m = re.match(r"GSM\d+_(.+)_peaks\.narrowPeak\.gz$", f.name)
        fname_label = fname_label_m.group(1) if fname_label_m else "?"
        try:
            with _gzip.open(f, "rt") as fh:
                first = fh.readline()
            internal = first.split("\t")[3] if first else ""
            internal_prefix = re.sub(r"_peak_\d+$", "", internal)
        except Exception as e:
            internal_prefix = f"READ_ERROR:{e}"
        # normalize case/underscore-only comparison; flag if the day number or
        # rep number visibly differs, OR if the internal label isn't even a
        # short sample tag (some files carry a full lab filesystem path instead).
        is_path_leak = "/" in internal_prefix
        # compare digits after 'D' in filename label vs internal label
        fname_day = re.search(r"[Dd](\d+)", fname_label)
        internal_day = re.search(r"[Dd](\d+)", internal_prefix)
        day_mismatch = (fname_day and internal_day and fname_day.group(1) != internal_day.group(1))
        fname_rep = re.search(r"[Rr]ep(\d+)", fname_label)
        internal_rep = re.search(r"[Rr]ep(\d+)", internal_prefix)
        rep_mismatch = (fname_rep and internal_rep and fname_rep.group(1) != internal_rep.group(1))
        flag = ""
        if is_path_leak:
            flag = f"internal peak name is a raw lab filesystem path, not a sample tag: {internal_prefix!r}"
        elif day_mismatch:
            flag = f"internal peak-name day ({internal_day.group(0)}) != filename day ({fname_day.group(0)})"
        elif rep_mismatch:
            flag = f"internal peak-name rep ({internal_rep.group(0)}) != filename rep ({fname_rep.group(0)})"
        mism[gsm] = flag
    return mism


def main():
    with open_text(SOFT, SOFT_GZ) as fh:
        series, samples = parse_soft(fh)

    if not samples:
        sys.exit("no ^SAMPLE blocks parsed -- check data/GSE206402/GSE206402_family.soft(.gz)")

    peak_flags = load_peak_label_mismatches()

    rows = []
    for s in samples:
        chars = chars_to_dict(s.get("characteristics", []))
        title = s.get("Sample_title", "")
        gsm = s.get("Sample_geo_accession", "")
        chars_text = " ; ".join(s.get("characteristics", []))
        treatment_raw = chars.get("treatment", "")
        row = {
            "gsm": gsm,
            "title": title,
            "cell_line": chars.get("cell line") or "UNRESOLVED -- no 'cell line' key",
            "condition": infer_condition(treatment_raw),
            "timepoint": infer_timepoint(title, chars_text),
            "treatment": treatment_raw or "UNRESOLVED -- no 'treatment' key",
            "replicate_id": infer_replicate(title),
            "raw_characteristics": chars_text,
            # extra audit columns (not in the minimum spec, kept for traceability)
            "cell_type": chars.get("cell type", ""),
            "source_name": s.get("Sample_source_name_ch1", ""),
            "library_strategy": s.get("Sample_library_strategy", ""),
            "instrument": s.get("Sample_instrument_model", ""),
            "genome_build_in_data_processing": (
                "hg19" if "hg19" in (
                    s.get("Sample_data_processing", "")
                    + " ".join(s.get("Sample_data_processing_extra", []))
                ) else "NOT FOUND in Sample_data_processing text"
            ),
            "supplementary_file_bw": s.get("Sample_supplementary_file_1", ""),
            "supplementary_file_narrowpeak": s.get("Sample_supplementary_file_2", ""),
            "peak_file_internal_label_flag": peak_flags.get(gsm, ""),
        }
        rows.append(row)

    fieldnames = list(rows[0].keys())
    out_csv = OUT / "gse206402_sample_metadata.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    n_flagged = sum(1 for r in rows if r["peak_file_internal_label_flag"])
    print(f"parsed {len(rows)} GSM rows -> {out_csv}", file=sys.stderr)
    print(f"{n_flagged}/{len(rows)} rows have a peak-file internal-label mismatch flag "
          f"(see peak_file_internal_label_flag column)", file=sys.stderr)
    unresolved = [r["gsm"] for r in rows if any(
        isinstance(v, str) and v.startswith("UNRESOLVED") for v in r.values())]
    if unresolved:
        print(f"UNRESOLVED fields present for GSMs: {unresolved}", file=sys.stderr)
    else:
        print("No UNRESOLVED fields -- title/characteristics cleanly parsed for all rows.", file=sys.stderr)


if __name__ == "__main__":
    main()
