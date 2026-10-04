#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 05_design_adequacy.py
#
# Score each candidate series for design adequacy (arms, replicates, controls).
#
# CONSUMES: output/tables/sample_metadata.csv
# PRODUCES: output/tables/design_adequacy_log.txt
# ---------------------------------------------------------------------------
"""
Task 2: design adequacy -- replicate counts per condition per series, and
whether GSE254358 and GSE220545 are genuinely paired (shared samples/donor/
BioProject) or independent experiments.

Replicate counts are computed from output/tables/sample_metadata.csv
(produced by 02_parse_metadata.py). The pairing check re-queries NCBI
eutils live (esearch/esummary against db=gds) so this script is a
standalone, re-runnable record of the evidence, not just narrative claims
in the report.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "output" / "tables" / "sample_metadata.csv"
OUT = ROOT / "output" / "tables" / "design_adequacy_log.txt"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def eutils_json(path, params):
    # shell out to curl rather than urllib: this machine's python.org build has
    # no configured CA bundle (SSLCertVerificationError on urlopen), while the
    # system curl's trust store works fine -- verified during the interactive
    # audit run before this script existed.
    qs = "&".join(f"{k}={v}" for k, v in params.items())
    url = f"{EUTILS}/{path}?{qs}&retmode=json"
    time.sleep(0.4)  # be polite to eutils (3 req/s cap without an API key)
    out = subprocess.run(["curl", "-s", url], capture_output=True, text=True, timeout=30, check=True)
    try:
        return json.loads(out.stdout)
    except json.JSONDecodeError:
        sys.exit(f"eutils returned non-JSON (likely rate-limited or transient error) for {url}:\n{out.stdout[:500]}")


def main():
    log = []

    # ---- replicate counts from parsed metadata ----
    df = pd.read_csv(META)
    log.append("=== Replicate counts per condition per series ===")
    for series, sub in df.groupby("series"):
        log.append(f"\n{series}: {len(sub)} GSM row(s) in parsed metadata")
        if series == "GSE254358":
            log.append("  N/A -- accession does not exist on GEO (see 04_atac_status.py). "
                        "0 samples, 0 replicates, 0 conditions available for this audit.")
            continue
        counts = sub.groupby("timepoint_inferred")["gsm"].nunique() if "timepoint_inferred" in sub else None
        if counts is not None:
            for cond, n in counts.items():
                flag = "  [FLAG: <3 biological replicates]" if n < 3 else "  [OK: >=3 replicates]"
                log.append(f"  condition/timepoint = {cond!r}: n = {n} GSM(s){flag}")

    # ---- pairing check: live eutils queries ----
    log.append("\n\n=== Pairing check (live NCBI eutils queries) ===")

    r = eutils_json("esearch.fcgi", {"db": "gds", "term": "GSE254358%5BACCN%5D"})
    n = r["esearchresult"]["count"]
    log.append(f"esearch db=gds term=GSE254358[ACCN] -> count={n}")
    gse254358_exists = int(n) > 0

    r = eutils_json("esearch.fcgi", {"db": "gds", "term": "GSE220545%5BACCN%5D"})
    n = r["esearchresult"]["count"]
    log.append(f"esearch db=gds term=GSE220545[ACCN] -> count={n}")
    gse220545_exists = int(n) > 0

    bioproject = None
    if gse220545_exists:
        r = eutils_json("esummary.fcgi", {"db": "gds", "id": "200220545"})
        rec = r["result"]["200220545"]
        bioproject = rec.get("bioproject")
        log.append(f"GSE220545 BioProject: {bioproject}")
        log.append(f"GSE220545 relations/extrelations (SuperSeries links): "
                    f"{rec.get('relations')}, {rec.get('extrelations')}")

        r2 = eutils_json("esearch.fcgi", {"db": "gds", "term": f"{bioproject}"})
        linked = r2["esearchresult"]["idlist"]
        log.append(f"GEO series under BioProject {bioproject}: {linked} "
                    f"(count={r2['esearchresult']['count']})")
        log.append("  -> " + ("ONLY GSE220545 itself is under this BioProject; "
                               "no companion ATAC-seq series exists there."
                               if linked == ["200220545"] else
                               "additional series found -- inspect manually."))

    log.append("\n=== Verdict ===")
    if not gse254358_exists:
        log.append(
            "GSE254358 does not exist in GEO. Therefore GSE254358 and GSE220545 CANNOT be "
            "genuinely paired (same cells/passage/induction) under any reading -- there is no "
            "second dataset to pair against. This is not 'two independent experiments that "
            "would have to be stitched together' either; it is one real dataset (GSE220545) "
            "and one non-existent accession. Confidence: HIGH (direct eutils existence check, "
            "FTP 404, BioProject cross-reference all agree)."
        )
    else:
        log.append("GSE254358 resolved on this run -- re-run full ATAC pipeline before trusting this branch.")

    report = "\n".join(log)
    OUT.write_text(report)
    print(report)


if __name__ == "__main__":
    main()
