#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 23_parse_gse74238_74324.py
#
# Parse GSE74238/GSE74324 metadata.
#
# CONSUMES: data/GSE74238/, data/GSE74324/ SOFT
# PRODUCES: output/tables/gse74238_sample_metadata.csv, gse74324_sample_metadata.csv
# ---------------------------------------------------------------------------
"""
Parse GSE74238 (ChIP-seq: H3K27ac/BRD4/input, IMR90 OIS) and GSE74324
(matched RNA-seq, 72 samples) family.soft files into tidy one-row-per-GSM
CSVs. Direct SOFT parsing (line-oriented ^ENTITY / !key = value format),
same approach as scripts/02_parse_metadata.py -- no GEOparse, no HTML
scraping (the suppl/ directory LISTING was scraped separately in
scripts/22_fetch_gse74238_74324.sh purely to get filenames; all metadata
here comes from the SOFT text itself).

Output:
  output/tables/gse74238_sample_metadata.csv  (22 rows)
  output/tables/gse74324_sample_metadata.csv  (72 rows)

Columns (per task spec): gsm, title, assay, cell_line, condition, treatment,
timepoint, replicate_id, raw_characteristics.

Design notes baked into the parsing logic (see inline comments):
  - GSE74238: "assay" = ChIP target (H3K27ac / BRD4 / input-no-antibody),
    parsed from the `antibody:` characteristic. "condition" = proliferating/
    quiescent/senescent/DMSO/Etoposide, parsed from the `treatment:`
    characteristic (DMSO/Etoposide are NOT part of the P/Q/S triad -- a
    separate genotoxic-stress arm, flagged in condition_group).
    replicate_id is inferred from the literal "replicate" suffix in
    Sample_title (absent -> rep 1, present -> rep 2) since GEO gives no
    explicit replicate-number field for this series.
  - GSE74324: "assay" is always RNA-seq. "condition" (P/Q/S) comes from the
    `cell condition:` characteristic. "treatment" captures BOTH axes GEO
    conflates under different characteristic keys for different sample
    groups -- shRNA identity from `viral transduction:` (e.g. "shBrd4",
    "shRen") AND, for the DMSO/JQ1 arms, the drug from `treatment:` (e.g.
    "DMSO 48h", "100nM JQ1 48h") -- both are recorded, pipe-separated, so
    the perturbation axis is fully auditable per row. replicate_id is the
    trailing _N in Sample_title (all groups are explicit triplicates).
  - Neither series has a per-sample GEO field for "timepoint" in the sense
    of hours/days post-induction; both protocols state harvest at "Day 12
    post-infection" as SERIES-LEVEL protocol text (Sample_growth_protocol_ch1,
    identical across all samples in a series), not a per-sample
    discriminating field. Recorded as such rather than invented per-row.
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
    """Return (series_dict, [sample_dict, ...]) from a family.soft file."""
    samples, series = [], {}
    cur, cur_type = None, None
    if not path.exists():
        return series, samples
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur is not None and cur_type == "SAMPLE":
                    samples.append(cur)
                cur, cur_type = {"characteristics": []}, "SAMPLE"
                continue
            if line.startswith("^SERIES"):
                cur, cur_type = series, "SERIES"
                continue
            if line.startswith("^PLATFORM"):
                if cur is not None and cur_type == "SAMPLE":
                    samples.append(cur)
                cur, cur_type = None, "PLATFORM"
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
                cur.setdefault("data_processing", []).append(val)
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


# ---------------------------------------------------------------- GSE74238 --

def build_gse74238_rows(series, samples):
    rows = []
    for s in samples:
        chars = chars_to_dict(s.get("characteristics", []))
        title = s.get("Sample_title", "")
        gsm = s.get("Sample_geo_accession", "")
        antibody = chars.get("antibody", "")
        if antibody.lower().startswith("h3k27ac"):
            assay = "H3K27ac (ChIP)"
        elif antibody.lower().startswith("brd4"):
            assay = "BRD4 (ChIP)"
        elif "none" in antibody.lower() or antibody.strip() in ("-", ""):
            assay = "input/no-antibody"
        else:
            assay = f"UNRESOLVED antibody={antibody!r}"

        treatment_raw = chars.get("treatment", "")
        t = treatment_raw.lower()
        if t in ("proliferating", "quiescent", "senescent"):
            condition = treatment_raw
            condition_group = "P/Q/S triad"
        elif t in ("dmso", "etoposide"):
            condition = treatment_raw
            condition_group = "DMSO/Etoposide genotoxic-stress arm (NOT part of P/Q/S triad)"
        else:
            condition = treatment_raw or "UNRESOLVED"
            condition_group = "UNRESOLVED"

        replicate_id = "2 (biological replicate; separate BioSample/SRX)" if "replicate" in title.lower() else "1"

        rows.append({
            "gsm": gsm,
            "title": title,
            "assay": assay,
            "cell_line": chars.get("cell line", "IMR90"),
            "condition": condition,
            "condition_group": condition_group,
            "treatment": treatment_raw,
            "timepoint": "NOT a per-sample field -- series-level protocol text says "
                         "'Cells were processed at day 12 post-infection' "
                         "(Sample_growth_protocol_ch1, identical across all 22 samples)",
            "replicate_id": replicate_id,
            "cell_type": chars.get("cell type", ""),
            "antibody_raw": antibody,
            "supplementary_file": s.get("Sample_supplementary_file_1", ""),
            "genome_build": next((dp.split(":", 1)[1].strip() for dp in s.get("data_processing", [])
                                   if dp.lower().startswith("genome_build")), ""),
            "raw_characteristics": " ; ".join(s.get("characteristics", [])),
        })
    return rows


# ---------------------------------------------------------------- GSE74324 --

PERTURBATION_LABELS = {
    "shren": "shRen (control shRNA, targets Renilla luciferase -- non-targeting control)",
    "shp53": "shp53 (p53 knockdown)",
    "shp65": "shp65 (RelA/p65 knockdown)",
    "shbrd4": "shBrd4 (BRD4 knockdown)",
    "shp16/p21": "shp16/p21 (CDKN2A+CDKN1A double knockdown)",
    "shp53/rb": "shp53/Rb (p53+RB double knockdown)",
}


def infer_perturbation(title):
    t = title.lower()
    # longest key first: "shp53" is a substring of "shp53/rb", so a naive
    # first-match-wins scan over dict-insertion order would silently
    # mislabel every shp53/Rb sample as plain shp53 -- caught by the
    # assay x condition sanity print in main() (shp53 showed n=6 instead
    # of 3 until this was sorted by key length).
    for key in sorted(PERTURBATION_LABELS, key=len, reverse=True):
        if key in t:
            return PERTURBATION_LABELS[key]
    if "dmso" in t:
        return "shRen background + DMSO vehicle (48h, control for JQ1 arm)"
    if "jq1" in t:
        return "shRen background + JQ1 100nM (48h, BRD4 bromodomain inhibitor)"
    return f"UNRESOLVED from title={title!r}"


def build_gse74324_rows(series, samples):
    rows = []
    for s in samples:
        chars = chars_to_dict(s.get("characteristics", []))
        title = s.get("Sample_title", "")
        gsm = s.get("Sample_geo_accession", "")

        condition = chars.get("cell condition", "")
        if not condition:
            tl = title.lower()
            condition = ("proliferating" if tl.startswith("proliferating") else
                         "quiescent" if tl.startswith("quiescent") else
                         "senescent" if tl.startswith("senescent") else "UNRESOLVED")

        viral = chars.get("viral transduction", "")
        drug_treatment = chars.get("treatment", "")  # "no treatment" / "DMSO 48h" / "100nM JQ1 48h"
        perturbation = infer_perturbation(title)
        is_baseline = "shren" in title.lower().split("_")[0].lower() or "shRen" in title
        is_baseline = "shren" in re.sub(r"[^a-z0-9/]", "", title.lower())
        is_clean_triad = bool(re.search(r"\bshRen_\d\b", title)) and drug_treatment.lower() == "no treatment"
        is_dmso_baseline = "dmso" in title.lower()

        m = re.search(r"_(\d+)$", title.strip())
        replicate_id = m.group(1) if m else ""

        rows.append({
            "gsm": gsm,
            "title": title,
            "assay": "RNA-seq",
            "cell_line": chars.get("cell line", "IMR90"),
            "condition": condition,
            "treatment": f"perturbation={perturbation} | drug_treatment={drug_treatment or 'NOT REPORTED'}",
            "perturbation_arm": perturbation,
            "is_clean_PQS_triad_baseline": is_clean_triad,  # shRen, "no treatment" -- closest match to GSE74238's P/Q/S
            "is_dmso_vehicle_baseline": is_dmso_baseline,   # shRen background + DMSO vehicle (2nd candidate baseline)
            "timepoint": "NOT a per-sample field -- series-level protocol text says "
                         "'RNA was extracted at Day 12 post-infection' "
                         "(Sample_growth_protocol_ch1, identical across all 72 samples); "
                         "DMSO/JQ1 arms additionally specify a 48h drug-exposure window ending at harvest",
            "replicate_id": replicate_id,
            "viral_transduction_raw": viral,
            "drug_treatment_raw": drug_treatment,
            "supplementary_file": s.get("Sample_supplementary_file_1", ""),
            "genome_build": next((dp.split(":", 1)[1].strip() for dp in s.get("data_processing", [])
                                   if dp.lower().startswith("genome_build")), ""),
            "raw_characteristics": " ; ".join(s.get("characteristics", [])),
        })
    return rows


def write_csv(rows, out_path):
    fieldnames = []
    for r in rows:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out_path} ({len(rows)} rows)", file=sys.stderr)


def main():
    series238, samples238 = parse_soft(DATA / "GSE74238" / "GSE74238_family.soft")
    if not samples238:
        sys.exit("[GSE74238] NO SAMPLES PARSED -- check data/GSE74238/GSE74238_family.soft exists "
                 "and is a real SOFT file (run scripts/22_fetch_gse74238_74324.sh first).")
    rows238 = build_gse74238_rows(series238, samples238)
    print(f"[GSE74238] parsed {len(rows238)} samples "
          f"(expected 22 per Series_sample_id count)", file=sys.stderr)
    write_csv(rows238, OUT / "gse74238_sample_metadata.csv")

    series324, samples324 = parse_soft(DATA / "GSE74324" / "GSE74324_family.soft")
    if not samples324:
        sys.exit("[GSE74324] NO SAMPLES PARSED -- check data/GSE74324/GSE74324_family.soft exists "
                 "and is a real SOFT file (run scripts/22_fetch_gse74238_74324.sh first).")
    rows324 = build_gse74324_rows(series324, samples324)
    print(f"[GSE74324] parsed {len(rows324)} samples "
          f"(expected 72 per Series_sample_id count)", file=sys.stderr)
    write_csv(rows324, OUT / "gse74324_sample_metadata.csv")

    # quick sanity summaries to stderr (full tables written to CSV)
    from collections import Counter
    print("\n[GSE74238] assay x condition counts:", file=sys.stderr)
    c = Counter((r["assay"], r["condition"]) for r in rows238)
    for k, v in sorted(c.items()):
        print(f"   {k}: {v}", file=sys.stderr)

    print("\n[GSE74324] perturbation_arm x condition counts:", file=sys.stderr)
    c2 = Counter((r["perturbation_arm"].split(" (")[0], r["condition"]) for r in rows324)
    for k, v in sorted(c2.items()):
        print(f"   {k}: {v}", file=sys.stderr)

    n_clean = sum(1 for r in rows324 if r["is_clean_PQS_triad_baseline"])
    n_dmso = sum(1 for r in rows324 if r["is_dmso_vehicle_baseline"])
    print(f"\n[GSE74324] clean shRen/'no treatment' P/Q/S baseline rows: {n_clean} / 72", file=sys.stderr)
    print(f"[GSE74324] shRen+DMSO-vehicle baseline rows (2nd candidate): {n_dmso} / 72", file=sys.stderr)


if __name__ == "__main__":
    main()
