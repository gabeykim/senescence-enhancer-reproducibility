#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# 04_atac_status.py
#
# Report which accessions actually carry usable ATAC data.
#
# CONSUMES: data/<ACC>/ filelists + metadata
# PRODUCES: stdout summary
# ---------------------------------------------------------------------------
"""
Task 3/4, ATAC-seq arm (GSE254358).

This script intentionally does no analysis. It exists so the audit has a
single reproducible entry point that documents *why* nothing was analyzed,
rather than the ATAC arm just silently being absent from the pipeline.

Verified 2026-08-20:
  - NCBI eutils esearch, db=gds, term=GSE254358[ACCN]  -> count=0
  - FTP https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254358/ -> HTTP 404
  - The GSE254nnn range folder itself resolves fine (HTTP 200), ruling out a
    path-construction bug -- the specific accession GSE254358 does not exist.
  - No GEO series is linked to GSE220545's BioProject (PRJNA910178) other than
    GSE220545 itself, and GSE220545 has no SuperSeries/relations pointing at
    a companion ATAC-seq series.
  - Checked 6 plausible single-digit-transposition typos of GSE254358
    (GSE254538, GSE245358, GSE253458, GSE254385, GSE245538, GSE253548): all
    exist as *unrelated* series (mouse myocarditis, MLL-AF9 leukemia, U87MG
    glioma, osteosarcoma FFPE -- none is ATAC-seq/senescence/IMR90). This is
    not proof no typo exists, but it rules out the simplest single-swap cases.

Conclusion: GSE254358 is not a resolvable GEO accession as given. Tasks 1/3/4
cannot be executed for the ATAC-seq arm from GEO metadata or supplementary
files, because there is no GEO record to pull them from. See the main report
for the impact on Task 2 (pairing) and Task 5 (verdict).
"""
print(__doc__)
