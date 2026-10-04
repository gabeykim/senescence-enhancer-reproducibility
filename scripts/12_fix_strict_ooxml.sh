#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 12_fix_strict_ooxml.sh
#
# Repair strict-OOXML xlsx that pandas/openpyxl cannot open.
#
# CONSUMES: data/GSE175533/*.xlsx
# PRODUCES: data/GSE175533/_xlsx_fixed*/
# ---------------------------------------------------------------------------
# Both GSE175533_atac_peaks_sig.xlsx AND GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx
# were saved by Excel/Numbers in "Strict OOXML" conformance (workbook.xml
# declares conformance="strict" and uses the http://purl.oclc.org/ooxml/...
# namespace URIs instead of the usual "Transitional"
# http://schemas.openxmlformats.org/... ones). openpyxl (and most other
# libraries) only understand Transitional OOXML, so pd.ExcelFile(...) silently
# returns zero sheets with no error -- easy to mistake for "the file is empty"
# rather than "the file uses a namespace variant this library doesn't parse."
#
# Fix: rewrite the Strict namespace URIs to their Transitional equivalents in
# every XML/rels part, then re-zip. This is a standard, lossless fix (the
# underlying cell/row/sharedStrings structure is identical between the two
# conformance classes -- only the XML namespace URIs and relationship-type
# URIs differ) -- verified by manually diffing Strict vs expected Transitional
# namespace strings in _rels/.rels and xl/_rels/workbook.xml.rels.
#
# Usage: 12_fix_strict_ooxml.sh <input.xlsx> <output.xlsx>
set -euo pipefail
SRC="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
OUT_DIR="$(cd "$(dirname "$2")" && pwd)"
OUT="$OUT_DIR/$(basename "$2")"
WORKDIR="$OUT_DIR/_xlsx_fixed_$$"

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
cd "$WORKDIR"
unzip -o -q "$SRC"

find . -name "*.xml" -o -name "*.rels" | while read -r f; do
  sed -i '' \
    -e 's#http://purl.oclc.org/ooxml/spreadsheetml/main#http://schemas.openxmlformats.org/spreadsheetml/2006/main#g' \
    -e 's#http://purl.oclc.org/ooxml/officeDocument/relationships#http://schemas.openxmlformats.org/officeDocument/2006/relationships#g' \
    -e 's#http://purl.oclc.org/ooxml/package/relationships#http://schemas.openxmlformats.org/package/2006/relationships#g' \
    -e 's# conformance="strict"##g' \
    "$f"
done

rm -f "$OUT"
zip -q -r -X "$OUT" '[Content_Types].xml' _rels docProps xl
cd - > /dev/null
rm -rf "$WORKDIR"
echo "wrote $OUT"
python3 -c "
import pandas as pd
xl = pd.ExcelFile('$OUT', engine='openpyxl')
print('sheets:', xl.sheet_names)
"
