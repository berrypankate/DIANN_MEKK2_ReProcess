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
# DERIVED PATHS
# ==============================================================================

STAGING="${PROJECT_DIR}/staging/${DATASET_ID}"
QUANT="${PROJECT_DIR}/quant/${DATASET_ID}"

RESULTS="${PROJECT_DIR}/results/${DATASET_ID}"
BATCH_RESULTS="${RESULTS}/batches"

LOG_DIR="${PROJECT_DIR}/logs/${DATASET_ID}"

STATE_DIR="${PROJECT_DIR}/work/state/${DATASET_ID}"
COMPLETED="${STATE_DIR}/completed.txt"

mkdir -p \
    "$STAGING" \
    "$QUANT" \
    "$BATCH_RESULTS" \
    "$LOG_DIR" \
    "$STATE_DIR"

touch "$COMPLETED"


# ==============================================================================
# VALIDATION
# ==============================================================================

[[ -f "$MANIFEST" ]] || {
    echo "ERROR: manifest not found:"
    echo "$MANIFEST"
    exit 1
}

[[ -f "$LICENSE_FILE" ]] || {
    echo "ERROR: license not found:"
    echo "$LICENSE_FILE"
    exit 1
}

[[ -f "${PROJECT_DIR}/${FASTA_FILE}" ]] || {
    echo "ERROR: FASTA not found:"
    echo "${PROJECT_DIR}/${FASTA_FILE}"
    exit 1
}

[[ -f "${PROJECT_DIR}/${LIBRARY_FILE}" ]] || {
    echo "ERROR: spectral library not found:"
    echo "${PROJECT_DIR}/${LIBRARY_FILE}"
    exit 1
}

if [[ "$MS1_PPM" == "0" || "$MS2_PPM" == "0" || "$SCAN_WINDOW" == "0" ]]; then
    echo
    echo "ERROR:"
    echo "MS1_PPM, MS2_PPM, and SCAN_WINDOW must be non-zero."
    echo
    echo "Current values:"
    echo "MS1_PPM=$MS1_PPM"
    echo "MS2_PPM=$MS2_PPM"
    echo "SCAN_WINDOW=$SCAN_WINDOW"
    exit 1
fi


# ==============================================================================
# BUILD MODIFICATION ARGUMENTS
# ==============================================================================

MOD_ARGS=()

if [[ "${USE_CARBAMIDOMETHYL:-false}" == "true" ]]; then
    MOD_ARGS+=(
        --fixed-mod "$CARBAMIDOMETHYL_MOD"
    )
fi

if [[ "${USE_PHOSPHO:-false}" == "true" ]]; then
    MOD_ARGS+=(
        --var-mod "$PHOSPHO_MOD"
    )
fi

if [[ "${USE_OXIDATION:-false}" == "true" ]]; then
    MOD_ARGS+=(
        --var-mod "$OXIDATION_MOD"
    )
fi

if [[ "${USE_NTERM_ACETYL:-false}" == "true" ]]; then
    MOD_ARGS+=(
        --var-mod "$NTERM_ACETYL_MOD"
    )
fi

if [[ "${USE_PHOSPHO:-false}" == "true" || \
      "${USE_OXIDATION:-false}" == "true" || \
      "${USE_NTERM_ACETYL:-false}" == "true" ]]; then

    MOD_ARGS+=(
        --var-mods "$MAX_VAR_MODS"
    )
fi


# ==============================================================================
# BUILD OPTIONAL FEATURE ARGUMENTS
# ==============================================================================

FEATURE_ARGS=()
# need to make use proteoforms to contain "true" so line 131 is true 
# and gives it back to 129 so that 132 

if [[ "${USE_PROTEOFORMS:-false}" == "true" ]]; then
    FEATURE_ARGS+=(--proteoforms)
fi

if [[ "${RELAXED_PROT_INF:-false}" == "true" ]]; then
    FEATURE_ARGS+=(--relaxed-prot-inf)
fi

if [[ "${NO_PROT_INF:-false}" == "true" ]]; then
    FEATURE_ARGS+=(--no-prot-inf)
fi
GENERATE_MATRICES="true"
# ==============================================================================
# CLEAN STALE RAW FILES
# ==============================================================================

echo "Cleaning staging directory..."

find "$STAGING" \
    -maxdepth 1 \
    -type f \
    -iname "*.raw" \
    -delete


# ==============================================================================
# BUILD REMAINING WORK LIST
# ==============================================================================

TODO_FILE=$(mktemp)

grep -Fvx -f "$COMPLETED" "$MANIFEST" > "$TODO_FILE" || true

TOTAL=$(grep -cv '^[[:space:]]*$' "$MANIFEST" || true)
DONE=$(grep -cv '^[[:space:]]*$' "$COMPLETED" || true)
TODO=$(grep -cv '^[[:space:]]*$' "$TODO_FILE" || true)


echo
echo "=========================================================="
echo "DIA-NN BATCH PROCESSING"
echo "=========================================================="
echo "Dataset:       $DATASET_ID"
echo "Total RAWs:    $TOTAL"
echo "Completed:     $DONE"
echo "Remaining:     $TODO"
echo "Batch size:    $BATCH_SIZE"
echo
echo "Library:"
echo "${PROJECT_DIR}/${LIBRARY_FILE}"
echo
echo "Mass accuracy:"
echo "  MS1:         ${MS1_PPM} ppm"
echo "  MS2:         ${MS2_PPM} ppm"
echo "  Scan window: ${SCAN_WINDOW}"
echo
echo "Mods:"
echo "  Carbamidomethyl: ${USE_CARBAMIDOMETHYL:-false}"
echo "  Phospho:         ${USE_PHOSPHO:-false}"
echo "  Oxidation:       ${USE_OXIDATION:-false}"
echo "  N-term acetyl:   ${USE_NTERM_ACETYL:-false}"
echo "=========================================================="

if [[ "$TODO" -eq 0 ]]; then
    echo "Nothing left to process."
    rm -f "$TODO_FILE"
    exit 0
fi


# ==============================================================================
# HELPER: MANIFEST ENTRY -> S3 URI
# ==============================================================================

make_s3_uri() {
    local ENTRY="$1"

    if [[ "$ENTRY" == s3://* ]]; then
        echo "$ENTRY"
    else
        echo "${S3_BUCKET}/${ENTRY}"
    fi
}


# ==============================================================================
# BATCH LOOP
# ==============================================================================

BATCH_NUMBER=0

while [[ -s "$TODO_FILE" ]]; do

    BATCH_NUMBER=$((BATCH_NUMBER + 1))

    BATCH_FILE=$(mktemp)

    head -n "$BATCH_SIZE" "$TODO_FILE" > "$BATCH_FILE"

    echo
    echo "=========================================================="
    echo "BATCH $BATCH_NUMBER"
    echo "=========================================================="

    cat "$BATCH_FILE"


    # ==========================================================================
    # DOWNLOAD CURRENT BATCH
    # ==========================================================================

    while IFS= read -r ENTRY; do

        [[ -z "$ENTRY" ]] && continue

        FILE=$(basename "$ENTRY")
        S3_URI=$(make_s3_uri "$ENTRY")

        echo
        echo "Downloading:"
        echo "$S3_URI"

        aws s3 cp \
            "$S3_URI" \
            "${STAGING}/${FILE}" \
            --only-show-errors

    done < "$BATCH_FILE"


    echo
    echo "Downloaded batch:"
    ls -lh "$STAGING"


    # ==========================================================================
    # BUILD --f ARGUMENTS
    # ==========================================================================

    RAW_ARGS=()

    while IFS= read -r ENTRY; do

        [[ -z "$ENTRY" ]] && continue

        FILE=$(basename "$ENTRY")

        RAW_ARGS+=(
            --f "/data/staging/${DATASET_ID}/${FILE}"
        )

    done < "$BATCH_FILE"


    # ==========================================================================
    # OUTPUT NAMES
    # ==========================================================================

    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

    OUTPUT="/data/results/${DATASET_ID}/batches/${DATASET_ID}_batch_${BATCH_NUMBER}_${TIMESTAMP}.parquet"

    LOG="${LOG_DIR}/${DATASET_ID}_batch_${BATCH_NUMBER}_${TIMESTAMP}.log"


    # ==========================================================================
    # RUN DIA-NN
    # ==========================================================================

    echo
    echo "Running DIA-NN..."
    echo

    set +e

    sudo docker run --rm \
        -v "${LICENSE_FILE}:${LICENSE_TARGET}:ro" \
        -v "${PROJECT_DIR}:/data" \
        "$DIANN_IMAGE" \
        "$DIANN_BIN" \
            "${RAW_ARGS[@]}" \
            --lib "/data/${LIBRARY_FILE}" \
            --fasta "/data/${FASTA_FILE}" \
            --temp "/data/quant/${DATASET_ID}" \
            --out "$OUTPUT" \
            --threads "$THREADS" \
            --mass-acc-ms1 "$MS1_PPM" \
            --mass-acc "$MS2_PPM" \
            --window "$SCAN_WINDOW" \
            --cut "$CUT" \
            --missed-cleavages "$MISSED_CLEAVAGES" \
            --min-pep-len "$MIN_PEPTIDE_LENGTH" \
            --max-pep-len "$MAX_PEPTIDE_LENGTH" \
            --min-pr-charge "$MIN_CHARGE" \
            --max-pr-charge "$MAX_CHARGE" \
            --quant-ori-names \
            --qvalue "$QVALUE" \
            --verbose "$VERBOSE" \
            "${MOD_ARGS[@]}" \
            "${FEATURE_ARGS[@]}" \
        2>&1 | tee "$LOG"

    STATUS=${PIPESTATUS[0]}

    set -e


    # ==========================================================================
    # STOP ON DIA-NN FAILURE
    # ==========================================================================

    if [[ "$STATUS" -ne 0 ]]; then

        echo
        echo "ERROR: DIA-NN failed."
        echo
        echo "RAW files are being retained."
        echo
        echo "Log:"
        echo "$LOG"

        rm -f "$BATCH_FILE" "$TODO_FILE"

        exit "$STATUS"
    fi


    # ==========================================================================
    # COLLECT / VERIFY .quant FILES
    # ==========================================================================

    echo
    echo "Collecting .quant files..."

    BATCH_OK=1

    while IFS= read -r ENTRY; do

        [[ -z "$ENTRY" ]] && continue

        FILE=$(basename "$ENTRY")

	# Remove .raw extension because DIA-NN writes:
	# sample.raw -> sample.quant
	STEM="${FILE%.raw}"

	DEST_QUANT="${QUANT}/${STEM}.quant"
	STAGING_QUANT="${STAGING}/${STEM}.quant"        

        # Already in persistent quant dir
        if [[ -f "$DEST_QUANT" ]]; then
            echo "OK: $DEST_QUANT"
            continue
        fi

        # DIA-NN may write next to RAW
        if [[ -f "$STAGING_QUANT" ]]; then

            echo "Moving .quant:"
            echo "  $STAGING_QUANT"
            echo "  -> $DEST_QUANT"

            mv "$STAGING_QUANT" "$DEST_QUANT"
        fi

        # Final verification
        if [[ ! -f "$DEST_QUANT" ]]; then

            echo
            echo "ERROR: missing .quant for:"
            echo "$FILE"

            BATCH_OK=0

        else

            echo "Verified: $DEST_QUANT"

        fi

    done < "$BATCH_FILE"


    # ==========================================================================
    # SAFETY STOP
    # ==========================================================================

    if [[ "$BATCH_OK" -ne 1 ]]; then

        echo
        echo "ERROR:"
        echo "One or more .quant files are missing."
        echo "RAW files will NOT be deleted."

        rm -f "$BATCH_FILE" "$TODO_FILE"

        exit 1
    fi


    # ==========================================================================
    # MARK SUCCESS
    # ==========================================================================

    cat "$BATCH_FILE" >> "$COMPLETED"

    sort -u "$COMPLETED" -o "$COMPLETED"


    # ==========================================================================
    # DELETE RAW FILES
    # ==========================================================================

    echo
    echo "All .quant files verified."
    echo "Deleting local RAW files..."

    while IFS= read -r ENTRY; do

        [[ -z "$ENTRY" ]] && continue

        FILE=$(basename "$ENTRY")

        rm -f "${STAGING}/${FILE}"

    done < "$BATCH_FILE"


    # ==========================================================================
    # REFRESH TODO LIST
    # ==========================================================================

    grep -Fvx -f "$COMPLETED" "$MANIFEST" \
        > "${TODO_FILE}.new" || true

    mv "${TODO_FILE}.new" "$TODO_FILE"

    rm -f "$BATCH_FILE"


    # ==========================================================================
    # STATUS
    # ==========================================================================

    echo
    echo "Disk usage:"
    df -h "$PROJECT_DIR"

    echo
    echo "Persistent .quant files:"
    find "$QUANT" \
        -maxdepth 1 \
        -type f \
        -name "*.quant" \
        | wc -l

done


rm -f "$TODO_FILE"

echo
echo "=========================================================="
echo "ALL ${DATASET_ID} BATCHES COMPLETE"
echo "=========================================================="
echo
echo "Quant directory:"
echo "$QUANT"
