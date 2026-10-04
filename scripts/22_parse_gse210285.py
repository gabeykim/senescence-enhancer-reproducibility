#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 22_parse_gse210285.py
#
# Parse GSE210285 metadata.
#
# CONSUMES: data/GSE210285/*_family.soft.gz
# PRODUCES: output/tables/gse210285_sample_metadata.csv
# ---------------------------------------------------------------------------
"""
Parse GSE210285_family.soft into a tidy one-row-per-GSM CSV.

Direct SOFT parsing (no GEOparse dependency), same approach as
scripts/02_parse_metadata.py in the sibling GSE254358/GSE220545 audit --
SOFT is a simple line-oriented key=value format with ^ENTITY markers.

GSE210285 is small (4 GSMs) and every !Sample_characteristics_ch1 line is
quoted verbatim into raw_characteristics per the task spec ("be
exhaustive"), rather than only keeping keyword-matched fields.
"""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOFT = ROOT / "data" / "GSE210285" / "GSE210285_family.soft"
OUT = ROOT / "output" / "tables" / "gse210285_sample_metadata.csv"
OUT.parent.mkdir(parents=True, exist_ok=True)


def parse_soft(path):
    """Return (series_dict, [sample_dict, ...]) from a family.soft file."""
    samples = []
    series = {}
    cur = None
    cur_type = None
    with open(path, encoding="utf-8", errors="replace") as fh:
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
            elif key == "Sample_data_processing":
                cur.setdefault("data_processing_all", []).append(val)
            elif key == "Sample_extract_protocol_ch1":
                cur.setdefault("extract_protocol_all", []).append(val)
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


def main():
    series, samples = parse_soft(SOFT)
    assert len(samples) == 4, f"expected 4 GSMs, parsed {len(samples)} -- re-check SOFT file"

    rows = []
    for s in samples:
        chars = chars_to_dict(s.get("characteristics", []))
        title = s.get("Sample_title", "")

        # cell line: explicit "cell line: 2BS" characteristic
        cell_line = chars.get("cell line", "")

        # condition / senescence status: explicit "phenotype: Growing" /
        # "phenotype: Senescence" characteristic -- no keyword inference
        # needed, GEO gives an unambiguous field for this series.
        condition = chars.get("phenotype", "")

        # passage number / population doubling level: NOT present as a
        # characteristics key, NOT present in title, NOT present in any
        # protocol field on this series. Explicitly flagged rather than
        # left blank so it reads as "checked and absent," not "forgot to
        # extract."
        passage_or_pdl = "NOT REPORTED IN GEO METADATA (no passage/PDL characteristic, " \
                          "title, or protocol field on any of the 4 GSMs)"

        # replicate id: explicit in Sample_title ("...-rep1"/"...-rep2")
        m = re.search(r"rep(\d+)\s*$", title, re.I)
        replicate_id = m.group(1) if m else ""

        raw_char = " ; ".join(f'"{c}"' for c in s.get("characteristics", []))

        row = {
            "gsm": s.get("Sample_geo_accession", ""),
            "title": title,
            "cell_line": cell_line,
            "condition_senescence_status": condition,
            "passage_or_pdl": passage_or_pdl,
            "replicate_id": replicate_id,
            "raw_characteristics": raw_char,
            # extra context columns, not in the minimum spec but useful for
            # the design-adequacy / cross-dataset sections below
            "source_name_ch1": s.get("Sample_source_name_ch1", ""),
            "organism_ch1": s.get("Sample_organism_ch1", ""),
            "library_strategy": s.get("Sample_library_strategy", ""),
            "instrument_model": s.get("Sample_instrument_model", ""),
            "platform_id": s.get("Sample_platform_id", ""),
            "extract_protocol_ch1_all": " || ".join(s.get("extract_protocol_all", [])),
            "data_processing_all": " || ".join(s.get("data_processing_all", [])),
            "supplementary_file_1": s.get("Sample_supplementary_file_1", ""),
            "biosample_relation": next(
                (r.split(": ", 1)[1] for r in s.get("Sample_relation_extra", []) if "BioSample" in r),
                s.get("Sample_relation", "") if "BioSample" in s.get("Sample_relation", "") else ""),
            "sra_relation": next(
                (r.split(": ", 1)[1] for r in s.get("Sample_relation_extra", []) if "SRA" in r),
                s.get("Sample_relation", "") if "SRA" in s.get("Sample_relation", "") else ""),
        }
        rows.append(row)

    # sort by GSM for a stable, predictable row order
    rows.sort(key=lambda r: r["gsm"])

    fieldnames = list(rows[0].keys())
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {OUT} ({len(rows)} rows)")

    # ---- quick design-adequacy printout (also re-derived formally in 23_*) ----
    from collections import Counter
    cond_counts = Counter(r["condition_senescence_status"] for r in rows)
    print("\nCondition counts (from 'phenotype' characteristic):")
    for cond, n in cond_counts.items():
        print(f"  {cond!r}: n={n}")

    platform_counts = Counter(r["platform_id"] for r in rows)
    print("\nPlatform counts:")
    for p, n in platform_counts.items():
        print(f"  {p}: n={n}")

    print("\nCondition x platform crosstab (checking for platform/condition confound):")
    for r in rows:
        print(f"  {r['gsm']} ({r['title']}): condition={r['condition_senescence_status']!r}, "
              f"platform={r['platform_id']}, instrument={r['instrument_model']}")


if __name__ == "__main__":
    main()
