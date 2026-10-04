# Placeholders to fill, in one pass

Fourteen locations across two files need the Zenodo DOI or the preprint link.
Nothing else in the repository contains a placeholder — `figures/captions.md`,
`reports/`, `audit/` and every script are clean, and no figure source points
outside the repository any more.

**Order matters.** Zenodo only archives GitHub releases created *after* the
repository is linked to Zenodo. Do not tag `v1.0.0` until the link exists, or
the release will not be archived.

Suggested sequence:

1. Deposit the large artifacts (staged separately; see that directory's
   `DEPOSITION_README.md` and `CHECKSUMS.sha256`) and record the concept DOI.
2. Link the GitHub repository to Zenodo in the Zenodo GitHub settings.
3. Fill the DOI placeholders below and commit.
4. Post the preprint, fill the link placeholders, commit.
5. Only then tag and publish the release, so Zenodo archives it.

---

## `README.md` — 6 DOI, 2 link

| Line | Current text | Replace with |
|---|---|---|
| 8 | `- **Manuscript:** _(link to be added on preprint posting)_` | the preprint URL |
| 9 | `- **Large artifacts (Zenodo):** _(DOI to be added on deposition — see …` | the Zenodo DOI |
| 121 | `… | Zenodo _(DOI pending)_, or regenerate with \`scripts/81_ois_cache.py\` …` | DOI |
| 122 | `… | Zenodo _(DOI pending)_ |` | DOI |
| 123 | `… | Zenodo _(DOI pending)_, or rerun steps 51/61/70 |` | DOI |
| 125 | `… | Zenodo _(DOI pending)_, or rerun steps 43/71 |` | DOI |
| 283 | `Please cite the manuscript _(link pending)_ and the Zenodo deposition` | preprint URL |
| 284 | `_(DOI pending)_ for the large artifacts. …` | DOI |

Line 124 (the raw GEO row) deliberately has no placeholder — that data is
re-fetched with the scripts, not downloaded from Zenodo.

## `CITATION.cff` — 8

| Line | Field | Current | Replace with |
|---|---|---|---|
| 49 | `version` | `0.0.0-pre-release` | the release tag, e.g. `1.0.0` |
| 50 | `date-released` | `2026-01-01` | the release date, `YYYY-MM-DD` |
| 51 | `doi` | `10.5281/zenodo.0000000` | the Zenodo **concept** DOI (resolves to the newest version) |
| 54 | `identifiers[0].value` | `10.5281/zenodo.0000000` | the Zenodo DOI (version-specific is fine here) |
| 55 | `identifiers[0].description` | `… (DOI pending)` | drop the "(DOI pending)" |
| 59 | `preferred-citation.title` | `PREPRINT TITLE PENDING` | the preprint title |
| 65 | `preferred-citation.doi` | `10.0000/preprint.pending` | the preprint DOI |
| 66 | `preferred-citation.url` | `https://example.org/preprint-link-pending` | the preprint URL |

After editing, re-check that the file still parses:

```bash
python3 -c "import yaml; yaml.safe_load(open('CITATION.cff')); print('ok')"
```

---

## Find them again

```bash
git grep -nE "DOI pending|link pending|link to be added|zenodo\.0000000|preprint\.pending|preprint-link-pending|PREPRINT TITLE PENDING|0\.0\.0-pre-release|2026-01-01"
```

Expected: 14 hits across `README.md` and `CITATION.cff`. When that command
returns nothing, this file can be deleted.
