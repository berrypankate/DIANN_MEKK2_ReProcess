#!/usr/bin/env bash

set -euo pipefail


# ==============================================================================
# LOAD CONFIG
# ==============================================================================

CONFIG="${1:?Usage: $0 config/<dataset>.env}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config not found: $CONFIG"
    exit 1
fi

source "$CONFIG"


# ==============================================================================
# PATHS
# ==============================================================================

QUANT="${PROJECT_DIR}/quant/${DATASET_ID}"

FINAL_DIR="${PROJECT_DIR}/results/${DATASET_ID}/final"

LOG_DIR="${PROJECT_DIR}/logs/${DATASET_ID}"

STATE_DIR="${PROJECT_DIR}/work/state/${DATASET_ID}"

COMPLETED="${STATE_DIR}/completed.txt"

mkdir -p "$FINAL_DIR" "$LOG_DIR"


# ==============================================================================
# VERIFY ALL MANIFEST FILES WERE PROCESSED
# ==============================================================================

TOTAL=$(grep -cv '^[[:space:]]*$' "$MANIFEST" || true)

DONE=$(grep -cv '^[[:space:]]*$' "$COMPLETED" || true)


echo
echo "=========================================================="
echo "FINALIZING $DATASET_ID"
echo "=========================================================="
echo "Manifest files:  $TOTAL"
echo "Completed files: $DONE"
echo


if [[ "$TOTAL" -ne "$DONE" ]]; then

    echo "ERROR:"
    echo "Not every manifest entry has completed batch processing."
    echo
    echo "Run:"
    echo "./bin/run_diann_batches.sh $CONFIG"
    echo

    exit 1
fi


QUANT_COUNT=$(find "$QUANT" \
    -maxdepth 1 \
    -type f \
    -name "*.quant" \
    | wc -l)


echo "Quant files:     $QUANT_COUNT"


if [[ "$QUANT_COUNT" -lt "$TOTAL" ]]; then

    echo
    echo "ERROR:"
    echo "There are fewer .quant files than manifest entries."
    echo
    echo "Expected >= $TOTAL"
    echo "Found       $QUANT_COUNT"
    echo

    exit 1
fi


# ==============================================================================
# BUILD THE ORIGINAL RAW FILE REFERENCES
#
# RAW files no longer need to exist if the matching .quant files are found.
# ==============================================================================

RAW_ARGS=()

while IFS= read -r ENTRY; do

    [[ -z "$ENTRY" ]] && continue

    FILE=$(basename "$ENTRY")

    RAW_ARGS+=(
        --f "/data/staging/${DATASET_ID}/${FILE}"
    )

done < "$MANIFEST"


# ==============================================================================
# OPTIONAL SETTINGS
# ==============================================================================


# ==============================================================================
# OUTPUT
# ==============================================================================

OUTPUT="/data/results/${DATASET_ID}/final/${DATASET_ID}_final.parquet"

LOG="${LOG_DIR}/${DATASET_ID}_final_$(date +%Y%m%d_%H%M%S).log"


# ==============================================================================
# FINAL DIA-NN PASS
# ==============================================================================

echo
echo "Running experiment-wide DIA-NN cross-run analysis..."
echo


sudo docker run --rm \
    -v "${LICENSE_FILE}:${LICENSE_TARGET}:ro" \
    -v "${PROJECT_DIR}:/data" \
    "$DIANN_IMAGE" \
    "$DIANN_BIN" \
        "${RAW_ARGS[@]}" \
        --lib "/data/${LIBRARY_FILE}" \
        --fasta "/data/${FASTA_FILE}" \
        --temp "/data/quant/${DATASET_ID}" \
        --use-quant \
        --out "$OUTPUT" \
        --threads "$THREADS" \
        --mass-acc-ms1 "$MS1_PPM" \
        --mass-acc "$MS2_PPM" \
        --window "$SCAN_WINDOW" \
        --cut "$CUT" \
        --missed-cleavages "$MISSED_CLEAVAGES" \
        --min-pep-len "$MIN_PEPTIDE_LENGTH" \
        --max-pep-len "$MAX_PEPTIDE_LENGTH" \
        --fixed-mod "$CARBAMIDOMETHYL_MOD" \
        --min-pr-charge "$MIN_CHARGE" \
        --max-pr-charge "$MAX_CHARGE" \
        --quant-ori-names \
        --qvalue "$QVALUE" \
        --verbose "$VERBOSE" \
        "${EXTRA_ARGS[@]}" \
    2>&1 | tee "$LOG"


STATUS=${PIPESTATUS[0]}


if [[ "$STATUS" -ne 0 ]]; then

    echo
    echo "ERROR: final DIA-NN analysis failed."
    echo
    echo "No .quant files have been deleted."
    echo
    echo "Log:"
    echo "$LOG"

    exit "$STATUS"
fi


echo
echo "=========================================================="
echo "FINAL ANALYSIS COMPLETE"
echo "=========================================================="
echo
echo "Main report:"
echo "${PROJECT_DIR}/results/${DATASET_ID}/final/${DATASET_ID}_final.parquet"
echo
echo "Keep the .quant files until QC is complete."
