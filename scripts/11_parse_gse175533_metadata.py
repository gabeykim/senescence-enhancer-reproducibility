#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 11_parse_gse175533_metadata.py
#
# Parse GSE175533 sample metadata.
#
# CONSUMES: data/GSE175533/*_family.soft.gz
# PRODUCES: output/tables/gse175533_sample_metadata.csv
# ---------------------------------------------------------------------------
"""
Task 1 for GSE175533: parse !Sample_characteristics_ch1 (and title/library
fields) into a tidy one-row-per-GSM dataframe covering all 147 samples.

Direct SOFT parsing, no GEOparse dependency (see 02_parse_metadata.py from
the prior audit in this project for the original rationale).

Design of this series (established by manual inspection before writing this
parser -- see the interactive audit transcript / REPORT.md Task 1 section for
the evidence trail):

  - Sample_title encodes everything structurally: an assay/arm prefix, a PDL
    or day or timepoint token, and a replicate letter (A/B/C, or _1/_2 for
    scRNA). Sample_library_strategy independently confirms ATAC-seq vs
    RNA-Seq per sample and agrees with the title prefix in 100% of samples
    checked -- used here as the authoritative assay-type field rather than
    inferring from title alone.
  - !Sample_characteristics_ch1 gives "population doublings: N" for the RS
    (replicative senescence) and hTERT arms only -- CD (cell density) and RIS
    (radiation-induced senescence) samples do NOT carry a PDL characteristic
    (they're staged by day post-treatment instead), so PDL is legitimately
    blank/NaN for those arms, not a parsing failure.
  - !Sample_extract_protocol_ch1 and !Sample_growth_protocol_ch1 contain
    boilerplate describing ALL arms (RS/hTERT, RIS, CD) concatenated
    identically on every single RNA-seq sample regardless of which arm it
    actually belongs to -- this is the same trap as GSE220545 in the prior
    audit (shared series-level protocol text copied onto every sample). Never
    used for per-sample condition inference here; condition comes from the
    structured Sample_title prefix instead.
"""
import re
import sys
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOFT = ROOT / "data" / "GSE175533" / "GSE175533_family.soft"
OUT = ROOT / "output" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def parse_soft(path):
    samples = []
    cur = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur is not None:
                    samples.append(cur)
                cur = {"characteristics": [], "data_processing": []}
                continue
            if line.startswith("^SERIES") or line.startswith("^PLATFORM"):
                if cur is not None:
                    samples.append(cur)
                    cur = None
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
                cur["data_processing"].append(val)
            elif key not in cur:
                cur[key] = val
    if cur is not None:
        samples.append(cur)
    return samples


def chars_to_dict(char_list):
    d = {}
    for c in char_list:
        if ":" in c:
            k, v = c.split(":", 1)
            d[k.strip().lower()] = v.strip()
    return d


# --- title parsing: arm, condition label, PDL/day, timepoint code, replicate ---
TITLE_PATTERNS = [
    # RS ATAC-seq: RS_atac_PDL20_TP1_A / RS_atac_hTERT_TP1_A
    (re.compile(r"^RS_atac_PDL(?P<pdl>\d+)_TP(?P<tp>\d+)_(?P<rep>[A-Z])$"),
     lambda m: dict(arm="RS_atac", cell_context="WI-38", pdl=int(m["pdl"]), day=None,
                    tp=f"TP{m['tp']}", replicate=m["rep"])),
    (re.compile(r"^RS_atac_hTERT_TP(?P<tp>\d+)_(?P<rep>[A-Z])$"),
     lambda m: dict(arm="RS_atac_hTERT", cell_context="WI-38 hTERT", pdl=None, day=None,
                    tp=f"TP{m['tp']}", replicate=m["rep"])),
    # RS bulk RNA: RS_PDL20_TP1_A
    (re.compile(r"^RS_PDL(?P<pdl>\d+)_TP(?P<tp>\d+)_(?P<rep>[A-Z])$"),
     lambda m: dict(arm="RS_bulkRNA", cell_context="WI-38", pdl=int(m["pdl"]), day=None,
                    tp=f"TP{m['tp']}", replicate=m["rep"])),
    # hTERT bulk RNA: hTERT_TP1_A
    (re.compile(r"^hTERT_TP(?P<tp>\d+)_(?P<rep>[A-Z])$"),
     lambda m: dict(arm="hTERT_bulkRNA", cell_context="WI-38 hTERT", pdl=None, day=None,
                    tp=f"TP{m['tp']}", replicate=m["rep"])),
    # RIS: RIS_d2_no_xray_A / RIS_d3_xray_A
    (re.compile(r"^RIS_d(?P<day>[\d.]+)_(?P<xray>no_xray|xray)_(?P<rep>[A-Z])$"),
     lambda m: dict(arm="RIS", cell_context="WI-38",
                    pdl=None, day=float(m["day"]),
                    tp=("mock-IR control" if m["xray"] == "no_xray" else "irradiated"),
                    replicate=m["rep"])),
    # CD: CD_D1.5_A
    (re.compile(r"^CD_D(?P<day>[\d.]+)_(?P<rep>[A-Z])$"),
     lambda m: dict(arm="CD", cell_context="WI-38", pdl=None, day=float(m["day"]),
                    tp=None, replicate=m["rep"])),
    # scRNA PDL: RS_sc_PDL_25_1
    (re.compile(r"^RS_sc_PDL_(?P<pdl>\d+)_(?P<rep>\d+)$"),
     lambda m: dict(arm="RS_scRNA", cell_context="WI-38", pdl=int(m["pdl"]), day=None,
                    tp=None, replicate=m["rep"])),
    # scRNA hTERT: RS_sc_hTert_t2
    (re.compile(r"^RS_sc_hTert_t(?P<tp>\d+)$"),
     lambda m: dict(arm="hTERT_scRNA", cell_context="WI-38 hTERT", pdl=None, day=None,
                    tp=f"TP{m['tp']}", replicate="1")),
]


def parse_title(title):
    for pat, fn in TITLE_PATTERNS:
        m = pat.match(title)
        if m:
            return fn(m.groupdict())
    return dict(arm="UNPARSED", cell_context=None, pdl=None, day=None, tp=None, replicate=None)


# senescence-condition labeling -- explicit, arm-by-arm, not a keyword guess.
# See REPORT Task 2(b) for the PDL cutoff justification (data-driven, from
# Task 4's marker/PCA results) -- this column intentionally does NOT assert a
# proliferating/senescent binary for RS/hTERT PDL points; that call is made
# downstream once the marker and PCA evidence is in hand, and encoded in a
# separate script (16_pdl_cutoff_and_diffacc.py) rather than hard-coded here.
ARM_NOTES = {
    "RS_atac": "WI-38 replicative senescence ATAC time course (PDL-staged) -- senescence status is a PDL continuum, not labeled binary in GEO",
    "RS_atac_hTERT": "hTERT-immortalized ATAC counterpart, matched timepoints -- not expected to senesce (telomerase-immortalized)",
    "RS_bulkRNA": "WI-38 replicative senescence bulk RNA-seq time course (PDL-staged)",
    "hTERT_bulkRNA": "hTERT-immortalized bulk RNA counterpart, matched timepoints",
    "RIS": "radiation-induced senescence (10 Gy X-ray) or its day-2 mock-IR control; RNA-seq ONLY, no ATAC exists for this arm",
    "CD": "cell-density/growth-arrest control (contact inhibition), RNA-seq ONLY, no ATAC exists for this arm",
    "RS_scRNA": "single-cell RNA-seq, WI-38 replicative senescence (not used in this audit's bulk QC)",
    "hTERT_scRNA": "single-cell RNA-seq, hTERT counterpart (not used in this audit's bulk QC)",
}


def main():
    samples = parse_soft(SOFT)
    rows = []
    for s in samples:
        title = s.get("Sample_title", "")
        parsed = parse_title(title)
        chars = chars_to_dict(s.get("characteristics", []))
        row = {
            "gsm": s.get("Sample_geo_accession", ""),
            "title": title,
            "assay": s.get("Sample_library_strategy", ""),
            "arm": parsed["arm"],
            "arm_note": ARM_NOTES.get(parsed["arm"], ""),
            "cell_line_or_context": parsed["cell_context"] or chars.get("cell line", ""),
            "pdl": parsed["pdl"] if parsed["pdl"] is not None else chars.get("population doublings", ""),
            "day_post_treatment": parsed["day"],
            "timepoint_code": parsed["tp"],
            "replicate_id": parsed["replicate"],
            "molecule": s.get("Sample_molecule_ch1", ""),
            "instrument": s.get("Sample_instrument_model", ""),
            "biosample": s.get("Sample_relation", ""),
            "raw_characteristics": " ; ".join(s.get("characteristics", [])),
        }
        rows.append(row)

    n_unparsed = sum(1 for r in rows if r["arm"] == "UNPARSED")
    print(f"Parsed {len(rows)} samples; {n_unparsed} titles did not match any known pattern "
          f"(should be 0 -- investigate if not).", file=sys.stderr)
    if n_unparsed:
        for r in rows:
            if r["arm"] == "UNPARSED":
                print("  UNPARSED:", r["gsm"], r["title"], file=sys.stderr)

    fieldnames = list(rows[0].keys())
    out_csv = OUT / "gse175533_sample_metadata.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out_csv}", file=sys.stderr)

    # quick inventory print for the audit log
    import collections
    by_arm_assay = collections.Counter((r["arm"], r["assay"]) for r in rows)
    print("\n=== Sample inventory by arm x assay ===", file=sys.stderr)
    for (arm, assay), n in sorted(by_arm_assay.items()):
        print(f"  {arm:20s} {assay:10s} n={n}", file=sys.stderr)


if __name__ == "__main__":
    main()
