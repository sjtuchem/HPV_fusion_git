#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${1:-${SCRIPT_DIR}/config.sh}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "[ERROR] Config file not found: ${CONFIG_FILE}" >&2
    exit 1
fi

# shellcheck disable=SC1090
source "${CONFIG_FILE}"

mkdir -p "${RESULT_DIR}/logs"
RUN_TS="$(date '+%Y%m%d_%H%M%S')"
LOG_FILE="${RESULT_DIR}/logs/run_all_${RUN_TS}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================================"
echo "On-target rates generic pipeline"
echo "Start time   : $(date '+%F %T')"
echo "Config       : ${CONFIG_FILE}"
echo "RAW_DIR      : ${RAW_DIR}"
echo "CODE_DIR     : ${CODE_DIR}"
echo "RESULT_DIR   : ${RESULT_DIR}"
echo "PRIMER_EXCEL : ${PRIMER_EXCEL}"
echo "SAMPLE_TABLE : ${SAMPLE_TABLE}"
echo "============================================================"

for f in "${PRIMER_EXCEL}" "${SAMPLE_TABLE}"; do
    if [[ ! -s "${f}" ]]; then
        echo "[ERROR] Required file missing or empty: ${f}" >&2
        exit 1
    fi
done

if [[ ! -d "${RAW_DIR}" ]]; then
    echo "[ERROR] RAW_DIR not found: ${RAW_DIR}" >&2
    exit 1
fi

# Step 1: PEAR -> A/B merge -> fastp
bash "${SCRIPT_DIR}/01_preprocess_expression.sh" "${CONFIG_FILE}"

# Step 2: iFP / iFP_spacer / X7(RC) count + on-target rates
# 避免在 set -u 下展开空数组导致“未绑定变量”。
if [[ "${FORCE_RERUN_COUNT:-0}" == "1" ]]; then
    python3 "${SCRIPT_DIR}/02_calculate_on_target_rates.py" \
        --fastp-dir "${RESULT_DIR}/03_fastp" \
        --out-dir "${RESULT_DIR}/04_on_target_rates" \
        --primer-excel "${PRIMER_EXCEL}" \
        --primer-sheet "${PRIMER_SHEET:-Sheet1}" \
        --samples-table "${SAMPLE_TABLE}" \
        --anchor "${ANCHOR_RC:-CTTGCAC}" \
        --force
else
    python3 "${SCRIPT_DIR}/02_calculate_on_target_rates.py" \
        --fastp-dir "${RESULT_DIR}/03_fastp" \
        --out-dir "${RESULT_DIR}/04_on_target_rates" \
        --primer-excel "${PRIMER_EXCEL}" \
        --primer-sheet "${PRIMER_SHEET:-Sheet1}" \
        --samples-table "${SAMPLE_TABLE}" \
        --anchor "${ANCHOR_RC:-CTTGCAC}"
fi

echo "============================================================"
echo "All done."
echo "End time : $(date '+%F %T')"
echo "Log      : ${LOG_FILE}"
echo "Results  : ${RESULT_DIR}"
echo "============================================================"
