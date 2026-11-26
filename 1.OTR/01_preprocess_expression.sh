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

PEAR_OUT_DIR="${RESULT_DIR}/01_pear"
MERGE_OUT_DIR="${RESULT_DIR}/02_merge_AB"
FASTP_OUT_DIR="${RESULT_DIR}/03_fastp"
LOG_DIR="${RESULT_DIR}/logs"
mkdir -p "${PEAR_OUT_DIR}" "${MERGE_OUT_DIR}" "${FASTP_OUT_DIR}" "${LOG_DIR}"

for cmd in pear fastp parallel gzip; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] ${cmd} not found in PATH" >&2
        exit 1
    }
done

if [[ ! -s "${SAMPLE_TABLE}" ]]; then
    echo "[ERROR] SAMPLE_TABLE missing or empty: ${SAMPLE_TABLE}" >&2
    exit 1
fi

# 输出启用样本（跳过表头、空行和#注释）
get_enabled_samples() {
    awk -F '\t' '
        NR==1 && $1=="sample_id" {next}
        /^[[:space:]]*#/ {next}
        NF==0 {next}
        $4=="1" || tolower($4)=="true" || tolower($4)=="yes" {print $1"\t"$2"\t"$3}
    ' "${SAMPLE_TABLE}"
}

find_one_fastq() {
    local sample_dir="$1"
    local mate="$2"
    local -a files=()

    mapfile -t files < <(
        find "${sample_dir}" -maxdepth 1 -type f \
            -name "*_R${mate}.fastq.gz" \
            ! -name "*.md5" \
            -print | sort
    )

    if [[ "${#files[@]}" -ne 1 ]]; then
        echo "[ERROR] Expected exactly 1 R${mate} FASTQ in: ${sample_dir}" >&2
        echo "[ERROR] Found ${#files[@]} file(s):" >&2
        printf '  %s\n' "${files[@]:-}" >&2
        return 1
    fi

    printf '%s\n' "${files[0]}"
}

run_pear_tube() {
    local sample_id="$1"
    local tube="$2"
    local rel_dir="$3"
    local sample_dir="${RAW_DIR}/${rel_dir}"

    if [[ ! -d "${sample_dir}" ]]; then
        echo "[ERROR] Sample directory not found: ${sample_dir}" >&2
        return 1
    fi

    local r1 r2
    r1="$(find_one_fastq "${sample_dir}" 1)"
    r2="$(find_one_fastq "${sample_dir}" 2)"

    local out_dir="${PEAR_OUT_DIR}/${sample_id}/${tube}"
    mkdir -p "${out_dir}"

    local prefix="${out_dir}/${sample_id}_${tube}"
    local final="${out_dir}/${sample_id}_${tube}_assembled_plus_unassembled_R1.fastq.gz"
    local assembled_gz="${prefix}.assembled.fastq.gz"
    local unassembled_fwd_gz="${prefix}.unassembled.forward.fastq.gz"

    echo "------------------------------------------------------------"
    echo "[PEAR] Sample=${sample_id} Tube=${tube}"
    echo "[PEAR] R1=${r1}"
    echo "[PEAR] R2=${r2}"
    echo "------------------------------------------------------------"

    # 通用跳过逻辑：最终 PEAR 单端文件已经存在且非空，就认为本管完成。
    # 不再因为中间 assembled/unassembled 文件缺失而强制重跑。
    if [[ "${FORCE_RERUN_PEAR:-0}" == "0" && -s "${final}" ]]; then
        echo "[SKIP] PEAR final file already exists: ${sample_id} ${tube}"
        printf '%s\n' "${final}"
        return 0
    fi

    rm -f \
        "${prefix}.assembled.fastq" "${prefix}.assembled.fastq.gz" \
        "${prefix}.unassembled.forward.fastq" "${prefix}.unassembled.forward.fastq.gz" \
        "${prefix}.unassembled.reverse.fastq" "${prefix}.unassembled.reverse.fastq.gz" \
        "${prefix}.discarded.fastq" "${prefix}.discarded.fastq.gz" \
        "${final}"

    pear \
        -f "${r1}" \
        -r "${r2}" \
        -o "${prefix}" \
        -v "${PEAR_MIN_OVERLAP:-30}" \
        -m "${PEAR_MAX_ASSEMBLY_LENGTH:-604}" \
        -n "${PEAR_MIN_ASSEMBLY_LENGTH:-30}" \
        -q "${PEAR_QUALITY_THRESHOLD:-20}" \
        -t "${PEAR_TRIM_LENGTH:-30}" \
        -u "${PEAR_P_VALUE:-0.01}" \
        -j "${THREADS:-20}" \
        2>&1 | tee "${out_dir}/${sample_id}_${tube}_pear.log"

    if [[ ! -e "${prefix}.assembled.fastq" ]]; then
        echo "[ERROR] PEAR assembled output missing: ${prefix}.assembled.fastq" >&2
        return 1
    fi
    if [[ ! -e "${prefix}.unassembled.forward.fastq" ]]; then
        echo "[ERROR] PEAR unassembled.forward output missing: ${prefix}.unassembled.forward.fastq" >&2
        return 1
    fi

    cat "${prefix}.assembled.fastq" "${prefix}.unassembled.forward.fastq" | gzip > "${final}"

    local f
    for f in \
        "${prefix}.assembled.fastq" \
        "${prefix}.unassembled.forward.fastq" \
        "${prefix}.unassembled.reverse.fastq" \
        "${prefix}.discarded.fastq"; do
        if [[ -e "${f}" ]]; then
            gzip -f "${f}"
        fi
    done

    if [[ ! -s "${final}" ]]; then
        echo "[ERROR] Empty PEAR final file: ${final}" >&2
        return 1
    fi

    echo "[DONE] PEAR ${sample_id} ${tube}: ${final}"
    printf '%s\n' "${final}"
}

process_one_sample() {
    local sample_id="$1"
    local a_rel_dir="$2"
    local b_rel_dir="$3"

    echo "============================================================"
    echo "[SAMPLE] ${sample_id}"
    echo "[A_DIR]  ${a_rel_dir}"
    echo "[B_DIR]  ${b_rel_dir}"
    echo "============================================================"

    # 固定输出路径，因此不依赖测序公司文件名前缀
    run_pear_tube "${sample_id}" "A" "${a_rel_dir}" >/dev/null
    run_pear_tube "${sample_id}" "B" "${b_rel_dir}" >/dev/null

    local a_file="${PEAR_OUT_DIR}/${sample_id}/A/${sample_id}_A_assembled_plus_unassembled_R1.fastq.gz"
    local b_file="${PEAR_OUT_DIR}/${sample_id}/B/${sample_id}_B_assembled_plus_unassembled_R1.fastq.gz"
    local merged_file="${MERGE_OUT_DIR}/${sample_id}_A_B_merged.fastq.gz"

    if [[ ! -s "${a_file}" || ! -s "${b_file}" ]]; then
        echo "[ERROR] PEAR A/B result missing for ${sample_id}" >&2
        return 1
    fi

    if [[ "${FORCE_RERUN_MERGE:-0}" == "0" && -s "${merged_file}" ]]; then
        echo "[SKIP] A+B merge already exists: ${merged_file}"
    else
        rm -f "${merged_file}"
        gzip -cd "${a_file}" "${b_file}" | gzip > "${merged_file}"
        if [[ ! -s "${merged_file}" ]]; then
            echo "[ERROR] A+B merged file is empty: ${merged_file}" >&2
            return 1
        fi
        echo "[DONE] A+B merged: ${merged_file}"
    fi

    local fastp_out="${FASTP_OUT_DIR}/${sample_id}_qc1.fastq.gz"
    local fastp_html="${FASTP_OUT_DIR}/${sample_id}_fastp_report.html"
    local fastp_json="${FASTP_OUT_DIR}/${sample_id}_fastp_report.json"
    local fastp_log="${FASTP_OUT_DIR}/${sample_id}_fastp.log"

    if [[ "${FORCE_RERUN_FASTP:-0}" == "0" && -s "${fastp_out}" ]]; then
        echo "[SKIP] fastp already exists: ${fastp_out}"
    else
        rm -f "${fastp_out}" "${fastp_html}" "${fastp_json}" "${fastp_log}"
        fastp \
            -i "${merged_file}" \
            -o "${fastp_out}" \
            -w "${THREADS:-20}" \
            --qualified_quality_phred "${FASTP_QUALIFIED_QUALITY_PHRED:-20}" \
            --unqualified_percent_limit "${FASTP_UNQUALIFIED_PERCENT_LIMIT:-20}" \
            --length_required "${FASTP_LENGTH_REQUIRED:-30}" \
            --length_limit "${FASTP_LENGTH_LIMIT:-604}" \
            --disable_adapter_trimming \
            -h "${fastp_html}" \
            -j "${fastp_json}" \
            2>&1 | tee "${fastp_log}"

        if [[ ! -s "${fastp_out}" ]]; then
            echo "[ERROR] fastp output missing or empty: ${fastp_out}" >&2
            return 1
        fi
        echo "[DONE] fastp: ${fastp_out}"
    fi

    echo "[DONE] Sample complete: ${sample_id}"
}

export -f find_one_fastq run_pear_tube process_one_sample
export RAW_DIR RESULT_DIR PEAR_OUT_DIR MERGE_OUT_DIR FASTP_OUT_DIR
export THREADS FORCE_RERUN_PEAR FORCE_RERUN_MERGE FORCE_RERUN_FASTP
export PEAR_MIN_OVERLAP PEAR_MAX_ASSEMBLY_LENGTH PEAR_MIN_ASSEMBLY_LENGTH
export PEAR_QUALITY_THRESHOLD PEAR_TRIM_LENGTH PEAR_P_VALUE
export FASTP_QUALIFIED_QUALITY_PHRED FASTP_UNQUALIFIED_PERCENT_LIMIT
export FASTP_LENGTH_REQUIRED FASTP_LENGTH_LIMIT
export SHELL=/bin/bash

ENABLED_COUNT="$(get_enabled_samples | awk 'END{print NR+0}')"
if [[ "${ENABLED_COUNT}" -eq 0 ]]; then
    echo "[ERROR] No enabled samples found in ${SAMPLE_TABLE}" >&2
    exit 1
fi

echo "============================================================"
echo "Preprocessing enabled samples: ${ENABLED_COUNT}"
echo "GNU parallel jobs=${PARALLEL_JOBS:-5}, PEAR/fastp threads=${THREADS:-20}"
echo "============================================================"

get_enabled_samples | \
    parallel --colsep '\t' --halt now,fail=1 \
        -j "${PARALLEL_JOBS:-5}" \
        --joblog "${LOG_DIR}/preprocess_parallel.joblog" \
        process_one_sample {1} {2} {3}

echo "============================================================"
echo "Preprocessing done."
echo "PEAR   : ${PEAR_OUT_DIR}"
echo "A+B    : ${MERGE_OUT_DIR}"
echo "fastp  : ${FASTP_OUT_DIR}"
echo "============================================================"
