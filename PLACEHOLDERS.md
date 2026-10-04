# Placeholders

**Zenodo DOI: DONE.** All 9 DOI placeholders were filled on 2026-10-04 with
[10.5281/zenodo.23146107](https://doi.org/10.5281/zenodo.23146107) — 6 in
`README.md`, 3 in `CITATION.cff`.

**7 placeholders remain**, in two groups that unblock at different times.
Neither can be filled yet.

---

## Group A — preprint-gated (5)

Fill when the preprint is posted. Nothing else blocks these.

| File | Line | Field / text | Replace with |
|---|---|---|---|
| `README.md` | 8 | `- **Manuscript:** _(link to be added on preprint posting)_` | the preprint URL |
| `README.md` | 283 | `the manuscript _(link pending; not yet posted)_` | the preprint URL |
| `CITATION.cff` | 63 | `preferred-citation.title` = `PREPRINT TITLE PENDING` | the preprint title |
| `CITATION.cff` | 69 | `preferred-citation.doi` = `10.0000/preprint.pending` | the preprint DOI |
| `CITATION.cff` | 70 | `preferred-citation.url` = `https://example.org/preprint-link-pending` | the preprint URL |

## Group B — release-gated (2)

These describe a GitHub **release**, which must not be created yet.

| File | Line | Field | Replace with |
|---|---|---|---|
| `CITATION.cff` | 53 | `version` = `0.0.0-pre-release` | the release tag, e.g. `1.0.0` |
| `CITATION.cff` | 54 | `date-released` = `2026-01-01` | the release date, `YYYY-MM-DD` |

**Ordering still matters.** Zenodo archives a GitHub release only if the
repository was linked to Zenodo *before* the release was created. The sequence
from here:

1. Link the repository to Zenodo in Zenodo's GitHub settings. *(The data
   deposit at 10.5281/zenodo.23146107 was uploaded manually and is unaffected
   by this; linking is only needed so that a future code release is archived.)*
2. Post the preprint and fill Group A.
3. Only then tag and publish the release, and fill Group B to match.

---

## Find them again

```bash
git grep -nE "DOI pending|DOI to be added|link pending|link to be added|zenodo\.0000000|preprint\.pending|preprint-link-pending|PREPRINT TITLE PENDING|0\.0\.0-pre-release|2026-01-01" -- README.md CITATION.cff
```

Expected now: **7 hits**. When it returns nothing, delete this file.

> Two corrections to the earlier version of this file, both found while filling
> the DOIs. It said "fourteen locations" when there were **16** (8 in each
> file), and its `git grep` one-liner omitted `DOI to be added`, so it silently
> missed `README.md` line 9 — one of the very DOI markers it was meant to
> track. The pattern above adds that alternative and restricts to the two files
> that carry placeholders, so the command no longer matches this file's own
> contents and inflates the count.

## On which DOI is used

`10.5281/zenodo.23146107` is the **version** DOI: it always resolves to this
exact deposit, whose checksums are the ones recorded in the deposit's own
`CHECKSUMS.sha256`. That is the right behaviour for a reproducibility record.

The **concept** DOI `10.5281/zenodo.23146106` resolves instead to whatever the
newest version is. If a v2 is ever deposited, the concept DOI would point at it
and the published checksums would no longer match. Swap it into
`CITATION.cff`'s `doi:` field only if you want "latest version" semantics there.
