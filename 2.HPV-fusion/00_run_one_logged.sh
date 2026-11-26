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
if [[ "$#" -lt 1 ]]; then
    echo "[ERROR] Missing sample argument" >&2
    exit 2
fi
SAMPLE_NAME="$1"
if [[ -z "${SAMPLE_NAME}" || "${SAMPLE_NAME}" == */* ]]; then
    echo "[ERROR] Unsafe sample name: ${SAMPLE_NAME}" >&2
    exit 2
fi

# Batch mode passes only the ASCII sample name through GNU parallel. This keeps
# Chinese filesystem paths out of GNU parallel's command template, which is
# required by older Perl/parallel builds that replace non-ASCII bytes with �.
if [[ "$#" -eq 1 ]]; then
    if [[ -z "${VALIDATED_JOBS_FILE:-}" || ! -s "${VALIDATED_JOBS_FILE}" ]]; then
        echo "[ERROR] VALIDATED_JOBS_FILE is required when only a sample name is supplied" >&2
        exit 2
    fi
    JOB_LINE="$(awk -F '\t' -v wanted="${SAMPLE_NAME}" '
        $1 == wanted { print; found=1; exit }
        END { if (!found) exit 1 }
    ' "${VALIDATED_JOBS_FILE}")" || {
        echo "[ERROR] Sample ${SAMPLE_NAME} is absent from ${VALIDATED_JOBS_FILE}" >&2
        exit 2
    }
    IFS=$'\t' read -r SAMPLE_NAME MERGED_FASTQ RAW_R1_A RAW_R2_A RAW_R1_B RAW_R2_B <<< "${JOB_LINE}"
    set -- "${SAMPLE_NAME}" "${MERGED_FASTQ}" "${RAW_R1_A}" "${RAW_R2_A}" "${RAW_R1_B}" "${RAW_R2_B}"
fi

SAMPLE_LOG_DIR="${RESULT_DIR}/${SAMPLE_NAME}/logs"
mkdir -p "${SAMPLE_LOG_DIR}"
bash "${CODE_DIR}/01_run_one_sample.sh" "$@" \
    2>&1 | tee -a "${SAMPLE_LOG_DIR}/${SAMPLE_NAME}.pipeline.log"
