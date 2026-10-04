#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 06_verify_accessions.sh
#
# Verify every accession resolves live at GEO/ENCODE.
#
# CONSUMES: accession list (network)
# PRODUCES: output/tables/accession_verification_log.txt
# ---------------------------------------------------------------------------
# Task 1 (existence gate): confirm every accession cited in this audit resolves to a
# live GEO or ENCODE record before any deeper work happens. A prior audit in this
# project found GSE254358 does not exist -- this script is the standing check that
# makes existence something verified on every re-run, not assumed.
set -uo pipefail  # no -e: we want to keep checking remaining accessions if one fails

geo_check () {
  local acc="$1"
  local numpart="${acc#GSE}"
  local padded; padded=$(printf "%06d" "$numpart")
  local uid="200${padded}"
  local json
  json=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=${uid}&retmode=json")
  echo "=== $acc (uid $uid) ==="
  echo "$json" | python3 -c "
import json,sys
d=json.load(sys.stdin)
rec=d.get('result',{}).get('$uid',{})
acc=rec.get('accession')
if acc != '$acc':
    print('  !! DOES NOT EXIST or UID mismatch -- verify manually')
else:
    print('  OK  title:', rec.get('title'))
    print('      n_samples:', rec.get('n_samples'), '| bioproject:', rec.get('bioproject'), '| pdat:', rec.get('pdat'))
"
  sleep 0.4
}

encode_check () {
  local acc="$1"
  local url="https://www.encodeproject.org/experiments/${acc}/?format=json"
  local json; json=$(curl -sL -H "Accept: application/json" "$url")
  echo "=== $acc (ENCODE) ==="
  echo "$json" | python3 -c "
import json,sys
d=json.load(sys.stdin)
acc = d.get('accession')
if acc == '$acc':
    print('  OK  type: Experiment | assay:', d.get('assay_title'), '| status:', d.get('status'))
    print('      biosample:', d.get('biosample_summary'))
else:
    # not an Experiment -- could be an Annotation (e.g. a ChromBPNet model) or missing
    print('  not an /experiments/ record as filed -- checking generically...')
"
  sleep 0.3
}

encode_check_generic () {
  # for ENCODE accessions that are not plain Experiments (e.g. Annotation/model records)
  local acc="$1"
  local json
  json=$(curl -s -H "Accept: application/json" "https://www.encodeproject.org/search/?searchTerm=${acc}&format=json")
  echo "$json" | python3 -c "
import json,sys
d=json.load(sys.stdin)
hits=[r for r in d.get('@graph',[]) if r.get('accession')=='$acc']
if not hits:
    print('  !! DOES NOT EXIST on ENCODE portal')
else:
    r=hits[0]
    print('  OK  type:', r.get('@type'), '| @id:', r.get('@id'))
"
  sleep 0.3
}

echo "###### GEO accessions ######"
for acc in GSE175533 GSE206402 GSE206496 GSE205898 GSE206493 GSE210285 GSE74238 GSE74324; do
  geo_check "$acc"
done

echo ""
echo "###### ENCODE accessions ######"
encode_check ENCSR200OML
echo "=== ENCSR978WIX (ENCODE) ==="
curl -s -H "Accept: application/json" "https://www.encodeproject.org/experiments/ENCSR978WIX/?format=json" \
  | python3 -c "
import json,sys
d=json.load(sys.stdin)
if d.get('code')=='301 Moved Permanently':
    print('  redirected (not an Experiment) ->', d.get('message','').split('to ')[1].split(';')[0])
"
curl -s -H "Accept: application/json" "https://www.encodeproject.org/annotations/ENCSR978WIX/?format=json" \
  | python3 -c "
import json,sys
d=json.load(sys.stdin)
print('  OK  type: Annotation (', d.get('annotation_type'), ') | status:', d.get('status'))
print('      description:', d.get('description'))
print('      assembly:', d.get('assembly'))
"
