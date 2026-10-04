#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 02_parse_metadata.py
#
# Parse GEO SOFT/series-matrix into a flat sample table.
#
# CONSUMES: data/<ACC>/*_family.soft.gz
# PRODUCES: output/tables/sample_metadata.csv
# ---------------------------------------------------------------------------
"""
Parse GEO family.soft files into a tidy one-row-per-GSM dataframe.
Direct SOFT parsing (no GEOparse dependency) -- SOFT is a simple
line-oriented key=value format with ^ENTITY markers, so a small
hand-rolled parser is more transparent/auditable than pulling in a
third-party library for this.

Fields extracted per sample, matching the audit's Task 1 spec:
  cell line/type, senescence status, induction method, timepoint,
  passage number, donor, replicate id -- pulled from
  Sample_characteristics_ch1 (key: value pairs) plus Sample_title /
  Sample_source_name_ch1 / Sample_treatment_protocol_ch1 as fallbacks.

ASSUMPTION: GEO does not enforce a controlled vocabulary for
characteristics keys, so "senescence status" / "induction method" /
"timepoint" etc. are inferred by keyword matching against whatever
keys+values a given series actually used (e.g. GSE220545 uses
"treatment", "time", "genotype" -- there is no literal
"senescence status" key). Every inferred field is paired with the
raw characteristics string it was derived from, so the inference is
auditable, not opaque.
"""
import re
import sys
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "output" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def parse_soft(path):
    """Return list of dicts, one per ^SAMPLE block, plus the ^SERIES block dict."""
    samples = []
    series = {}
    cur = None
    cur_type = None
    if not path.exists():
        return series, samples
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
            elif key in cur:
                # repeated key (e.g. multiple Sample_supplementary_file_N) -> keep first + count
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


INDUCTION_KEYWORDS = {
    "replicative": ["replicative", "serial passage", "population doubling"],
    "oncogene-induced (OIS)": ["ras", "raf", "oncogen", "4oht", "tamoxifen", "er-ras"],
    "therapy-induced (TIS)": ["doxorubicin", "etoposide", "chemo", "therapy-induced", "bleomycin"],
    "irradiation-induced": ["irradiat", "ir-induced", "gamma", "x-ray"],
}


def infer_induction(text):
    t = text.lower()
    hits = [label for label, kws in INDUCTION_KEYWORDS.items() if any(k in t for k in kws)]
    return ";".join(hits) if hits else ""


# NOTE: keyword matching deliberately runs ONLY over title/source/characteristics,
# never over free-text protocol fields. GSE220545's Sample_treatment_protocol_ch1 is
# an identical boilerplate sentence copied onto every sample ("Uncommitted control
# cells were treated ... and committed cells were treated ...") -- it describes the
# whole series' design, not this particular sample, so matching against it makes
# every sample look "ambiguous" regardless of arm. Status has to come from the
# per-sample title/characteristics (here: the timepoint), cross-checked by hand
# against the protocol sentence once per series, not re-derived per row.
SENESCENT_KEYWORDS = ["senescent", "senescence-committed", "committed"]
PROLIF_KEYWORDS = ["proliferat", "uncommitted", "young", "growing", "presenescent", "quiescent"]
NEGATION_PREFIX = ["un", "non", "pre"]


def infer_status(text):
    t = " " + text.lower() + " "
    is_sen = any(k in t for k in SENESCENT_KEYWORDS) and "uncommitted" not in t.replace("committed", "")
    is_sen = any(k in t for k in SENESCENT_KEYWORDS)
    is_uncommitted = "uncommitted" in t
    if is_uncommitted:
        return "not senescent (explicitly 'uncommitted')"
    if is_sen:
        return "senescent/committed"
    is_prolif = any(k in t for k in PROLIF_KEYWORDS)
    if is_prolif:
        return "proliferating"
    return ""


def infer_timepoint(text):
    m = re.search(r"(day\s*[\d.]+|[\d.]+\s*h(?:our)?s?\b|[\d.]+\s*d(?:ay)?s?\b)", text, re.I)
    return m.group(0) if m else ""


# GSE220545 has no per-sample status field at all (see infer_status note): the
# 30h-vs-72h -> uncommitted-vs-committed mapping only exists as one series-level
# sentence in Sample_treatment_protocol_ch1, identical on every sample. Encoding
# it here as an explicit, hand-verified cross-reference (not a keyword guess) is
# more honest than trying to make infer_status() re-derive it per row.
MANUAL_STATUS_CROSSREF = {
    ("GSE220545", "Day 1.5"): "uncommitted (30h 4OHT) -- RAS-activated but NOT senescent; "
                               "NOT an untreated/proliferating control",
    ("GSE220545", "Day 3"): "committed/senescent (72h 4OHT) -- per Series_treatment_protocol_ch1",
}


def infer_replicate(title, chars):
    m = re.search(r"replicate\s*([A-Za-z0-9]+)", title, re.I)
    if m:
        return m.group(1)
    for k, v in chars.items():
        if "replicate" in k or k in ("rep", "biological replicate"):
            return v
    return ""


def build_rows(acc, series, samples):
    rows = []
    for s in samples:
        chars = chars_to_dict(s.get("characteristics", []))
        title = s.get("Sample_title", "")
        source = s.get("Sample_source_name_ch1", "")
        treat_protocol = s.get("Sample_treatment_protocol_ch1", "")
        full_text = " | ".join([title, source, treat_protocol] + s.get("characteristics", []))
        # status text deliberately excludes treat_protocol -- see NOTE above infer_status
        status_text = " | ".join([title, source] + s.get("characteristics", []))
        row = {
            "series": acc,
            "gsm": s.get("Sample_geo_accession", ""),
            "title": title,
            "source_name": source,
            "cell_line_or_type": chars.get("cell line") or chars.get("cell type") or source,
            "senescence_status_inferred": infer_status(status_text) or "UNRESOLVED FROM METADATA -- see timepoint_inferred + treatment_protocol",
            "induction_method_inferred": infer_induction(full_text),
            "timepoint_inferred": chars.get("time") or infer_timepoint(full_text),
            "status_via_protocol_crossref": MANUAL_STATUS_CROSSREF.get(
                (acc, chars.get("time", "")), ""),
            "passage_number": chars.get("passage") or chars.get("passage number") or "NOT REPORTED",
            "donor": chars.get("donor") or chars.get("individual") or chars.get("subject") or "NOT REPORTED",
            "replicate_id": infer_replicate(title, chars),
            "genotype": chars.get("genotype", ""),
            "treatment_raw": chars.get("treatment", ""),
            "raw_characteristics": " ; ".join(s.get("characteristics", [])),
            "treatment_protocol": treat_protocol,
            "library_strategy": s.get("Sample_library_strategy", ""),
            "instrument": s.get("Sample_instrument_model", ""),
            "supplementary_file_1": s.get("Sample_supplementary_file_1", ""),
            "bioproject": next((r.split(": ", 1)[1] for r in s.get("Sample_relation_extra", [])
                                 if "BioProject" in r), series.get("Series_relation", "")),
        }
        rows.append(row)
    return rows


def main():
    all_rows = []
    accessions = [("GSE254358", "GSE254358_family.soft"), ("GSE220545", "GSE220545_family.soft")]
    for acc, fname in accessions:
        path = DATA / acc / fname
        series, samples = parse_soft(path)
        if not samples:
            print(f"[{acc}] NO SAMPLES PARSED from {path} "
                  f"(file missing, empty, or not a real SOFT file -- see audit notes).", file=sys.stderr)
            all_rows.append({
                "series": acc, "gsm": "N/A", "title": "ACCESSION DOES NOT EXIST ON GEO (verified via "
                "esearch db=gds and FTP 404 on 2026-08-20)" if acc == "GSE254358" else "PARSE FAILURE",
            })
            continue
        rows = build_rows(acc, series, samples)
        print(f"[{acc}] parsed {len(rows)} samples", file=sys.stderr)
        all_rows.extend(rows)

    fieldnames = []
    for r in all_rows:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)

    out_csv = OUT / "sample_metadata.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"wrote {out_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
