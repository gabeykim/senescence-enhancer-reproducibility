# ⚠ GSE254358 DOES NOT EXIST — THESE FILES ARE NOT DATA

**Do not load, parse, glob, or cite any file in this directory as a dataset.**

This directory is kept deliberately, as evidence. It holds what NCBI returns when you request
a GEO accession that does not exist. It is referenced by the provenance note in the manuscript
and by `output/number_audit/NUMBER_AUDIT.md`.

## What the accession is

`GSE254358` was the accession originally supplied as the ATAC-seq arm of a senescence
sequence-to-activity fine-tune, to be paired with `GSE220545` (RNA-seq). It is **absent from
GEO entirely** — not private, not embargoed, and not a path-construction mistake. No data from
it exists anywhere in this project, and no result depends on it.

## How that was verified

Verified independently on 2026-08-20 (`REPORT.md`, top level) and **re-verified live on
2026-10-04** during the number audit. All four checks reproduce:

| Check | Command / URL | Result |
|---|---|---|
| GEO DataSets index | `esearch db=gds term=GSE254358[ACCN]` | **`count = 0`**, empty `idlist` |
| GEO FTP series path | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254358/` | **HTTP 404** |
| Parent range directory | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/` | **HTTP 200** — the range folder is real, so the 404 is not a path bug |
| Adjacent sibling accession | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254357/` | **HTTP 200** — an accession one lower in the same folder resolves, so the 404 is specific to GSE254358 |
| BioProject link | `esearch db=bioproject term=GSE254358` | **`count = 0`** |

Additionally, `GSE220545`'s BioProject (`PRJNA910178`) was checked for a linked companion
series and has none — GSE220545 is the only GEO record under it. Six plausible
single-digit-transposition typos (`GSE254538`, `GSE245358`, `GSE253458`, `GSE254385`,
`GSE245538`, `GSE253548`) were each checked: all six resolve to **real but unrelated** series
(mouse myocarditis, MLL-AF9 leukemia, U87MG TMZ resistance, osteosarcoma FFPE — none is
ATAC-seq, senescence, or IMR90). That rules out the simple typo cases without proving no typo
exists.

## What the files actually are

All three were produced by `scripts/01_fetch_geo_metadata.sh` requesting URLs that returned 404.
The fetch did not fail loudly, so the error bodies were written to disk under data-looking names.

| File | Size | True type | md5 | sha256 |
|---|---|---|---|---|
| `GSE254358_family.soft.gz` | 990 B | **XHTML, not gzip** — an Apache `Error 404 / Object not found!` page from `ftp.ncbi.nlm.nih.gov` | `af4cfff178d8056c25c148f96e742185` | `04b10f58e3f1340fc897c7057a0001c82ebd216d8ced1a96a2f05676b5d22494` |
| `GSE254358_series_matrix.txt.gz` | 990 B | **XHTML, not gzip** — byte-identical to the file above | `af4cfff178d8056c25c148f96e742185` | `04b10f58e3f1340fc897c7057a0001c82ebd216d8ced1a96a2f05676b5d22494` |
| `suppl_listing.txt` | 66 B | plain text; two copies of the line `mailto:webadmin@ncbi.nlm.nih.gov`, scraped out of the same 404 page | `58320622d970ec45733291510b51c32a` | `1bcac4a29f61d5cc2ceec0b2fd707a27a942544e88115342eaedb70ea514a29e` |

Note the two `.gz` files are **not compressed at all** and are byte-identical to each other.
`file` reports them as `XML 1.0 document text, ASCII text`; `gzcat` fails with
`not in gzip format`. An earlier description of them as "decompressing to a mailto address" was
imprecise — there is nothing to decompress.

## Why this matters

This is the project's demonstrated provenance failure mode: a 404 body saved under a `.gz`
filename inside a `data/<accession>/` directory, which a directory glob or a naive loader would
treat as a dataset. Checksums are recorded above so these specific files can be cited as the
example without re-fetching.

Verified 2026-10-04 · see `output/audit_fixes/` (item 7) and `output/audit_fixes/accession_log_merged.csv`.
