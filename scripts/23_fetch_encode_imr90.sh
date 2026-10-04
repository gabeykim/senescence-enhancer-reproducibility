#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 23_fetch_encode_imr90.sh
#
# Fetch the ENCODE IMR-90 ATAC experiment and annotation records.
#
# CONSUMES: encodeproject.org (network)
# PRODUCES: data/ENCSR200OML/, data/ENCSR978WIX/
# ---------------------------------------------------------------------------
# Fetch ENCODE portal JSON + small processed files for the IMR-90 ATAC-seq
# baseline (ENCSR200OML) and its ChromBPNet model annotation (ENCSR978WIX).
#
# Scope, per audit constraints: REST/JSON API only, no raw FASTQ/BAM.
# Processed peak/region files are downloaded only if file_size (per the
# experiment/annotation JSON) is comfortably under ~200MB; anything larger
# is left un-downloaded and just noted (size/format) in the report.
#
# Requires: curl, python3, tar, gunzip on PATH.

set -euo pipefail
cd "$(dirname "$0")/.."   # -> project root

DATA_ENC1="data/ENCSR200OML"
DATA_ENC2="data/ENCSR978WIX"
mkdir -p "$DATA_ENC1" "$DATA_ENC2"

echo "== Fetching experiment/annotation JSON =="
curl -sL -H "Accept: application/json" \
  "https://www.encodeproject.org/experiments/ENCSR200OML/?format=json" \
  -o "$DATA_ENC1/experiment.json"
curl -sL -H "Accept: application/json" \
  "https://www.encodeproject.org/annotations/ENCSR978WIX/?format=json" \
  -o "$DATA_ENC2/annotation.json"

python3 -c "import json; json.load(open('$DATA_ENC1/experiment.json')); print('experiment.json OK')"
python3 -c "import json; json.load(open('$DATA_ENC2/annotation.json')); print('annotation.json OK')"

echo "== Downloading small processed peak/region files (ENCSR200OML) =="
# ENCFF243NTP: preferred_default=True, output_type "pseudoreplicated peaks"
#   (pooled-pseudoreplicate overlap peak set; the ENCODE4-designated primary/
#   "optimal"-equivalent peak call for this experiment). ~6.5MB.
curl -sL "https://www.encodeproject.org/files/ENCFF243NTP/@@download/ENCFF243NTP.bed.gz" \
  -o "$DATA_ENC1/ENCFF243NTP_pseudoreplicated_peaks.bed.gz"
# ENCFF982UNH: "conservative IDR thresholded peaks" (true-replicate IDR, strict). ~3.7MB.
curl -sL "https://www.encodeproject.org/files/ENCFF982UNH/@@download/ENCFF982UNH.bed.gz" \
  -o "$DATA_ENC1/ENCFF982UNH_conservative_idr_peaks.bed.gz"
# ENCFF114GDS: "IDR thresholded peaks" (pooled-pseudoreplicate IDR). ~3.8MB.
curl -sL "https://www.encodeproject.org/files/ENCFF114GDS/@@download/ENCFF114GDS.bed.gz" \
  -o "$DATA_ENC1/ENCFF114GDS_idr_thresholded_peaks.bed.gz"

echo "== Downloading small processed files (ENCSR978WIX / ChromBPNet) =="
# ENCFF815MUM: "selected regions for predicted signal and sequence contribution
#   scores" (bed3, ~1.5MB) -- the region set ChromBPNet reports full per-base
#   contribution scores for (a subset, not genome-wide).
curl -sL "https://www.encodeproject.org/files/ENCFF815MUM/@@download/ENCFF815MUM.bed.gz" \
  -o "$DATA_ENC2/ENCFF815MUM_selected_regions.bed.gz"
# ENCFF802KJC: "training and test regions" tar (~72MB) -- per-fold peak/nonpeak
#   train/valid/test BED splits used to train the 5 ChromBPNet folds.
curl -sL "https://www.encodeproject.org/files/ENCFF802KJC/@@download/ENCFF802KJC.tar.gz" \
  -o "$DATA_ENC2/ENCFF802KJC_training_test_regions.tar.gz"

echo "== Verifying downloads =="
for f in "$DATA_ENC1"/*.bed.gz "$DATA_ENC2"/*.bed.gz; do
  gunzip -t "$f" && echo "OK (gzip valid): $f"
done
tar tzf "$DATA_ENC2/ENCFF802KJC_training_test_regions.tar.gz" > /dev/null \
  && echo "OK (tar valid): $DATA_ENC2/ENCFF802KJC_training_test_regions.tar.gz"

ls -la "$DATA_ENC1" "$DATA_ENC2"
echo "Done."
