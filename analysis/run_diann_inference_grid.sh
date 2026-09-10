#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="/home/ubuntu/external_data/DIANN_MEKK2_ReProcess"

PF_CONFIG="${PROJECT_DIR}/config/grid/MEKK2_PF.env"
NPI_CONFIG="${PROJECT_DIR}/config/grid/MEKK2_NPI.env"

cd "$PROJECT_DIR"

echo "============================================================"
echo "MEKK2 DIA-NN PROTEIN-INFERENCE BENCHMARK"
echo "============================================================"

# ============================================================
# 1. PROTEOFORMS
# ============================================================

echo
echo "============================================================"
echo "WORKFLOW:      PROTEOFORMS"
echo "EXPERIMENT ID: MEKK_PF"
echo "FLAG:          --proteoforms"
echo "============================================================"

bash bin/run_diann_batches.sh "$PF_CONFIG"
# contents of pf_config gonna go as args to script im running, syntax always command then arguments
# bash command, first arg name of script u run in bash, everything after foll

echo
echo "Finalizing MEKK_PF..."

bash bin/finalize_diann.sh "$PF_CONFIG" 

echo
echo "MEKK_PF COMPLETE"


# ============================================================
# 2. NO PROTEIN INFERENCE
# ============================================================

echo
echo "============================================================"
echo "WORKFLOW:      NO PROTEIN INFERENCE"
echo "EXPERIMENT ID: MEKK2_NPI"
echo "FLAG:          --no-prot-inf"
echo "============================================================"

bash bin/run_diann_batches.sh "$NPI_CONFIG"

echo
echo "Finalizing MEKK2_NPI..."

bash bin/finalize_diann.sh "$NPI_CONFIG"

echo
echo "MEKK2_NPI COMPLETE"


# ============================================================
# OUTPUT CHECK
# ============================================================

echo
echo "============================================================"
echo "FINAL OUTPUTS"
echo "============================================================"

ls -lh \
    results/MEKK_PF/final/MEKK_PF_final.parquet \
    results/MEKK2_NPI/final/MEKK2_NPI_final.parquet

echo
echo "ALL INFERENCE WORKFLOWS COMPLETE"