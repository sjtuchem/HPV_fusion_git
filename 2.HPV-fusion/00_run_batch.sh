#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${SCRIPT_DIR}/config.sh}"
if [[ ! -s "${CONFIG_FILE}" ]]; then
    echo "[ERROR] Config missing or empty: ${CONFIG_FILE}" >&2
    exit 1
fi
# shellcheck disable=SC1090
source "${CONFIG_FILE}"

INPUT_SAMPLE_LIST="${1:-${SAMPLE_LIST}}"
for command_name in python3 parallel; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "[ERROR] Required command not found: ${command_name}" >&2
        exit 1
    }
done
if ! [[ "${PARALLEL_JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] PARALLEL_JOBS must be a positive integer: ${PARALLEL_JOBS}" >&2
    exit 1
fi
if ! [[ "${THREADS_PER_SAMPLE}" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] THREADS_PER_SAMPLE must be a positive integer: ${THREADS_PER_SAMPLE}" >&2
    exit 1
fi

BATCH_ID="$(date '+%Y%m%d_%H%M%S')_$$"
BATCH_DIR="${RESULT_DIR}/_batch_logs/${BATCH_ID}"
mkdir -p "${BATCH_DIR}"
JOBS_FILE="${BATCH_DIR}/validated_jobs.tsv"
SAMPLE_NAMES_FILE="${BATCH_DIR}/sample_names.txt"
VALIDATION_REPORT="${BATCH_DIR}/sample_list_validation.json"
cp "${INPUT_SAMPLE_LIST}" "${BATCH_DIR}/sample.list.snapshot"

python3 "${CODE_DIR}/00_validate_sample_list.py" \
    --sample-list "${INPUT_SAMPLE_LIST}" \
    --output "${JOBS_FILE}" \
    --report "${VALIDATION_REPORT}"
SAMPLE_COUNT="$(wc -l < "${JOBS_FILE}" | tr -d ' ')"
cut -f 1 "${JOBS_FILE}" > "${SAMPLE_NAMES_FILE}"

echo "============================================================"
echo "[PIPELINE]           HPV-human fusion V4.2.0 targeted UMI A+B batch"
echo "[SAMPLE LIST]        ${INPUT_SAMPLE_LIST}"
echo "[SAMPLES]            ${SAMPLE_COUNT}"
echo "[PARALLEL JOBS]      ${PARALLEL_JOBS}"
echo "[THREADS / SAMPLE]   ${THREADS_PER_SAMPLE}"
echo "[MAX ALIGN THREADS]  $((PARALLEL_JOBS * THREADS_PER_SAMPLE))"
echo "[RESULT]             ${RESULT_DIR}/<sample>/"
echo "[BATCH LOG]          ${BATCH_DIR}"
echo "============================================================"

set +e
export CONFIG_FILE FORCE
export VALIDATED_JOBS_FILE="${JOBS_FILE}"
(
    cd "${CODE_DIR}"
    parallel \
        --jobs "${PARALLEL_JOBS}" \
        --line-buffer \
        --tagstring '[{}]' \
        --joblog "${BATCH_DIR}/parallel.joblog" \
        bash ./00_run_one_logged.sh {} \
        < "${SAMPLE_NAMES_FILE}"
) 2>&1 | tee "${BATCH_DIR}/batch.log"
PARALLEL_STATUS="${PIPESTATUS[0]}"
set -e

set +e
python3 "${CODE_DIR}/00_batch_summary.py" \
    --jobs "${JOBS_FILE}" \
    --result-dir "${RESULT_DIR}" \
    --joblog "${BATCH_DIR}/parallel.joblog" \
    --output-tsv "${BATCH_DIR}/batch_summary.tsv" \
    --output-json "${BATCH_DIR}/batch_summary.json"
SUMMARY_STATUS="$?"
set -e

echo "[SUMMARY] ${BATCH_DIR}/batch_summary.tsv"
echo "[JOBLOG]  ${BATCH_DIR}/parallel.joblog"
if [[ "${PARALLEL_STATUS}" -ne 0 || "${SUMMARY_STATUS}" -ne 0 ]]; then
    echo "[ERROR] One or more samples failed or remain incomplete; inspect the batch summary and per-sample logs." >&2
    exit 1
fi
echo "[DONE] All listed samples are COMPLETE_PASS"
