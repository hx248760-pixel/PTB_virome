#!/bin/bash
# ==============================================================================
# PTB Virome Project
# script/run_rf.sh
# Run random forest: feature importance, K-fold CV, model building, prediction.
# Called by 04c_random_forest.R, but can also be run standalone.
#
# Usage:  bash script/run_rf.sh <project_dir>
# ==============================================================================
set -euo pipefail

PROJECT_DIR="${1:-.}"
PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"

# ------------------------------------------------------------------------------
# Paths (mirror config.R)
# ------------------------------------------------------------------------------
RF_DIR="${PROJECT_DIR}/results/random_forest"
RF_PY="${PROJECT_DIR}/script/rf.py"
PYTHON="${RF_PYTHON:-python3}"
THREADS="${RF_THREADS:-10}"
N_EST="${RF_N_EST:-1000}"
K="${RF_K:-5}"

# Input files prepared by 04c_random_forest.R
VIR_PROFILE="${RF_DIR}/sig_vir_spe.profile"
BAC_PROFILE="${RF_DIR}/sig_bac_spe.profile"
MERGE_PROFILE="${RF_DIR}/merged.profile"
SAMPLE_INFO="${RF_DIR}/sample_info.txt"

# ------------------------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------------------------
for f in "${VIR_PROFILE}" "${BAC_PROFILE}" "${MERGE_PROFILE}" "${SAMPLE_INFO}"; do
    if [ ! -f "${f}" ]; then
        echo "[RF] Missing required file: ${f}" >&2
        echo "     Please run 04c_random_forest.R first to generate inputs." >&2
        exit 1
    fi
done

if [ ! -f "${RF_PY}" ]; then
    echo "[RF] rf.py not found at ${RF_PY}" >&2
    exit 1
fi

mkdir -p "${RF_DIR}"
cd "${RF_DIR}"

# ------------------------------------------------------------------------------
# 0. Seed list
# ------------------------------------------------------------------------------
if [ ! -f seed.list ]; then
    seq 1 10 > seed.list
fi

# ------------------------------------------------------------------------------
# 1. Feature importance (--fast)
# ------------------------------------------------------------------------------
echo "[RF] Computing feature importance..."
"${PYTHON}" "${RF_PY}" importance -i "${VIR_PROFILE}"   -g "${SAMPLE_INFO}" \
    -o vir   --fast --threads "${THREADS}" --nest "${N_EST}"
"${PYTHON}" "${RF_PY}" importance -i "${BAC_PROFILE}"   -g "${SAMPLE_INFO}" \
    -o bac   --fast --threads "${THREADS}" --nest "${N_EST}"
"${PYTHON}" "${RF_PY}" importance -i "${MERGE_PROFILE}" -g "${SAMPLE_INFO}" \
    -o merge --fast --threads "${THREADS}" --nest "${N_EST}"

# ------------------------------------------------------------------------------
# 2. K-Fold cross-validation predictions
# ------------------------------------------------------------------------------
mkdir -p vir_kf bac_kf merge_kf

TOP_N_VIR=$(seq 2 5 100; seq 101 50 806)
TOP_N_BAC=$(seq 2 5 81)
TOP_N_MERGE=$(seq 2 5 100; seq 101 50 887)

echo "[RF] Running K-Fold for virus..."
parallel -j "${THREADS}" "${PYTHON}" "${RF_PY}" KF \
    -i vir.sorted -g "${SAMPLE_INFO}" -K "${K}" --stratified \
    -o vir_kf/predict_top{1}_seed{2}.tsv --topN {1} --seed {2} \
    ::: ${TOP_N_VIR} :::: seed.list

echo "[RF] Running K-Fold for bacteria..."
parallel -j "${THREADS}" "${PYTHON}" "${RF_PY}" KF \
    -i bac.sorted -g "${SAMPLE_INFO}" -K "${K}" --stratified \
    -o bac_kf/predict_top{1}_seed{2}.tsv --topN {1} --seed {2} \
    ::: ${TOP_N_BAC} :::: seed.list

echo "[RF] Running K-Fold for merged..."
parallel -j "${THREADS}" "${PYTHON}" "${RF_PY}" KF \
    -i merge.sorted -g "${SAMPLE_INFO}" -K "${K}" --stratified \
    -o merge_kf/predict_top{1}_seed{2}.tsv --topN {1} --seed {2} \
    ::: ${TOP_N_MERGE} :::: seed.list

# ------------------------------------------------------------------------------
# 3. Merge K-Fold outputs into single auc files
# ------------------------------------------------------------------------------
merge_kf_output() {
    local indir="$1"
    local outfile="$2"

    local first_file
    first_file=$(ls "${indir}"/predict_top*_seed*.tsv 2>/dev/null \
                 | sort -V | head -n 1)

    if [ -z "${first_file}" ]; then
        echo "[RF] No K-Fold files found in ${indir}" >&2
        return 1
    fi

    head -1 "${first_file}" > "${outfile}"
    for f in "${indir}"/predict_top*_seed*.tsv; do
        tail -n +2 "${f}"
    done >> "${outfile}"
}

echo "[RF] Merging K-Fold results..."
merge_kf_output vir_kf   vir_auc_kf.txt
merge_kf_output bac_kf   bac_auc_kf.txt
merge_kf_output merge_kf merge_auc_kf.txt

echo "[RF] Done. Outputs in ${RF_DIR}"
