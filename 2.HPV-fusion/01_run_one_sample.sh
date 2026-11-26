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

if [[ "${VERIFY_RAW_MD5}" != "auto" && "${VERIFY_RAW_MD5}" != "0" && "${VERIFY_RAW_MD5}" != "1" ]]; then
    echo "[ERROR] VERIFY_RAW_MD5 must be auto, 0 or 1: ${VERIFY_RAW_MD5}" >&2
    exit 2
fi
if [[ "${FORCE}" != "0" && "${FORCE}" != "1" ]]; then
    echo "[ERROR] FORCE must be 0 or 1: ${FORCE}" >&2
    exit 2
fi
if [[ "${REQUIRE_UMI}" != "0" && "${REQUIRE_UMI}" != "1" ]]; then
    echo "[ERROR] REQUIRE_UMI must be 0 or 1: ${REQUIRE_UMI}" >&2
    exit 2
fi

if [[ "$#" -ne 4 && "$#" -ne 6 ]]; then
    echo "Usage: bash 01_run_one_sample.sh SAMPLE MERGED_FASTQ RAW_R1_A RAW_R2_A [RAW_R1_B RAW_R2_B]" >&2
    exit 2
fi
SAMPLE="$1"
MERGED_FASTQ="$2"
R1_A="$3"
R2_A="$4"
R1_B="${5:-}"
R2_B="${6:-}"
[[ "${R1_B}" == "." ]] && R1_B=""
[[ "${R2_B}" == "." ]] && R2_B=""
if [[ -z "${SAMPLE}" || "${SAMPLE}" == "." || "${SAMPLE}" == ".." || "${SAMPLE}" == */* ]]; then
    echo "[ERROR] Unsafe sample name: ${SAMPLE}" >&2
    exit 2
fi
if { [[ -z "${R1_B}" ]] && [[ -n "${R2_B}" ]]; } || \
   { [[ -n "${R1_B}" ]] && [[ -z "${R2_B}" ]]; }; then
    echo "[ERROR] RAW_R1_B and RAW_R2_B must both be supplied or both be omitted" >&2
    exit 2
fi
if [[ "${REQUIRE_UMI}" == "1" && -z "${R1_B}" ]]; then
    echo "[ERROR] This targeted UMI release requires both raw tube A and raw tube B" >&2
    exit 2
fi

PIPELINE_VERSION="4.3.0-targeted-UMI-AB"
THREADS="${THREADS_PER_SAMPLE}"
SAMPLE_ROOT="${RESULT_DIR}/${SAMPLE}"
WORK_DIR="${SAMPLE_ROOT}/work"
FINAL_DIR="${SAMPLE_ROOT}/final"
AUDIT_DIR="${SAMPLE_ROOT}/audit"
LOG_DIR="${SAMPLE_ROOT}/logs"
STATE_DIR="${LOG_DIR}/state"
STRICT_CALL_DIR="${AUDIT_DIR}/strict"
RESCUE_DIR="${AUDIT_DIR}/technical_rescue_v4_1"
RESCUE_CALL_DIR="${RESCUE_DIR}/call"
mkdir -p \
    "${WORK_DIR}" "${FINAL_DIR}" "${AUDIT_DIR}" "${LOG_DIR}" "${STATE_DIR}" \
    "${WORK_DIR}/star" "${STRICT_CALL_DIR}" "${RESCUE_DIR}/star" "${RESCUE_CALL_DIR}"

if command -v flock >/dev/null 2>&1; then
    exec 9> "${SAMPLE_ROOT}/.pipeline.lock"
    if ! flock -n 9; then
        echo "[ERROR] Sample ${SAMPLE} is already running (lock: ${SAMPLE_ROOT}/.pipeline.lock)" >&2
        exit 1
    fi
fi

TIMING_FILE="${LOG_DIR}/${SAMPLE}.runtime_metrics.tsv"
RUN_CONFIG_FILE="${LOG_DIR}/${SAMPLE}.run_config.txt"
RUN_MANIFEST="${LOG_DIR}/${SAMPLE}.run_manifest.json"
CURRENT_MANIFEST="${LOG_DIR}/.${SAMPLE}.run_manifest.current.$$.json"
COMPLETION_FILE="${LOG_DIR}/${SAMPLE}.complete"
FINAL_EVENTS="${FINAL_DIR}/${SAMPLE}.fusion_events.final.tsv"
FINAL_SUMMARY="${FINAL_DIR}/${SAMPLE}.fusion_summary.final.json"

record_timing() {
    local stage_name="$1"
    local stage_status="$2"
    local stage_start="$3"
    local stage_end
    stage_end="$(date +%s)"
    printf '%s\t%s\t%d\n' "${stage_name}" "${stage_status}" "$((stage_end - stage_start))" >> "${TIMING_FILE}"
}

stage_ready() {
    local marker="$1"
    shift
    [[ "${FORCE}" != "1" && -f "${marker}" ]] || return 1
    local value
    for value in "$@"; do
        [[ -s "${value}" ]] || return 1
    done
}

mark_complete() {
    local marker="$1"
    printf 'completed_at=%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" > "${marker}"
}

final_ready() {
    local required
    local final_outputs=(
        "${FINAL_DIR}/${SAMPLE}.fusion_events.final.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_events.main.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_events.candidate.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_events.review.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_loci.final.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_loci.main.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_loci.candidate.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_loci.review.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_evidence.final.tsv.gz"
        "${FINAL_DIR}/${SAMPLE}.fusion_unique_sequences.final.tsv.gz"
        "${FINAL_DIR}/${SAMPLE}.fusion_sequence_members.final.tsv.gz"
        "${FINAL_DIR}/${SAMPLE}.fusion_umi_families.exact.tsv.gz"
        "${FINAL_DIR}/${SAMPLE}.fusion_umi_families.final.tsv.gz"
        "${FINAL_DIR}/${SAMPLE}.fusion_umi_families.corrected.tsv.gz"
        "${FINAL_DIR}/${SAMPLE}.fusion_filter_summary.tsv"
        "${FINAL_DIR}/${SAMPLE}.fusion_summary.final.json"
        "${LOG_DIR}/11_validate_results.log"
    )
    for required in "${final_outputs[@]}"; do
        [[ -s "${required}" ]] || return 1
    done
    grep -q '"status": "PASS"' "${LOG_DIR}/11_validate_results.log"
}

for command_name in python3; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "[ERROR] Required command not found: ${command_name}" >&2
        exit 1
    }
done

required_files=(
    "${CODE_DIR}/00_build_run_manifest.py"
    "${CODE_DIR}/00_verify_md5.py"
    "${CODE_DIR}/00_validate_umi_summary.py"
    "${CODE_DIR}/fusion_common.py"
    "${CODE_DIR}/01_prepare_unique_sequences.py"
    "${CODE_DIR}/02_rotate_hpv_reference.py"
    "${CODE_DIR}/02_select_hpv_candidates.py"
    "${CODE_DIR}/03_call_and_aggregate_fusions.py"
    "${CODE_DIR}/04_prepare_technical_rescue.py"
    "${CODE_DIR}/05_finalize_fusion_results.py"
    "${CODE_DIR}/98_validate_results.py"
    "${MERGED_FASTQ}" "${R1_A}" "${R2_A}"
    "${HPV_REF}" "${GENE_BED}" "${HPV_TYPE_MAP}"
)
if [[ -n "${R1_B}" ]]; then
    required_files+=("${R1_B}" "${R2_B}")
fi
for required_file in "${required_files[@]}"; do
    if [[ ! -s "${required_file}" ]]; then
        echo "[ERROR] Required file missing or empty: ${required_file}" >&2
        exit 1
    fi
done

manifest_args=(
    --sample "${SAMPLE}" --pipeline-version "${PIPELINE_VERSION}"
    --input "merged_fastq=${MERGED_FASTQ}"
    --input "raw_r1_a=${R1_A}" --input "raw_r2_a=${R2_A}"
    --reference "hpv_ref=${HPV_REF}"
    --reference "gene_bed=${GENE_BED}"
    --reference "hpv_type_map=${HPV_TYPE_MAP}"
    --parameter "prefilter_min_hpv=${PREFILTER_MIN_HPV}"
    --parameter "prefilter_min_flank=${PREFILTER_MIN_FLANK}"
    --parameter "thresholds=${THRESHOLDS}"
    --parameter "primary_threshold=${PRIMARY_THRESHOLD}"
    --parameter "min_human_mapq=${MIN_HUMAN_MAPQ}"
    --parameter "min_hpv_mapq=${MIN_HPV_MAPQ}"
    --parameter "max_junction_gap=${MAX_JUNCTION_GAP}"
    --parameter "full_human_fraction=${FULL_HUMAN_FRACTION}"
    --parameter "max_missing_raw_id_rate=${MAX_MISSING_RAW_ID_RATE}"
    --parameter "main_min_raw_reads=${MAIN_MIN_RAW_READS}"
    --parameter "main_min_umi_families=${MAIN_MIN_UMI_FAMILIES}"
    --parameter "candidate_min_raw_reads=${CANDIDATE_MIN_RAW_READS}"
    --parameter "candidate_min_umi_families=${CANDIDATE_MIN_UMI_FAMILIES}"
    --parameter "enable_technical_rescue=${ENABLE_TECHNICAL_RESCUE}"
    --parameter "rescue_terminal_slack=${RESCUE_TERMINAL_SLACK}"
    --parameter "rescue_min_trimmed_length=${RESCUE_MIN_TRIMMED_LENGTH}"
    --parameter "locus_human_window=${LOCUS_HUMAN_WINDOW}"
    --parameter "locus_hpv_window=${LOCUS_HPV_WINDOW}"
    --parameter "umi_correction_distance=${UMI_CORRECTION_DISTANCE}"
    --parameter "require_umi=${REQUIRE_UMI}"
)
if [[ -n "${R1_B}" ]]; then
    manifest_args+=(--input "raw_r1_b=${R1_B}" --input "raw_r2_b=${R2_B}")
fi
for reference_file in \
    "${HPV_REF}.amb" "${HPV_REF}.ann" "${HPV_REF}.bwt" "${HPV_REF}.pac" "${HPV_REF}.sa" \
    "${STAR_INDEX}/Genome" "${STAR_INDEX}/SA" "${STAR_INDEX}/SAindex"; do
    manifest_args+=(--reference "reference_$(basename "${reference_file}")=${reference_file}")
done
for script_file in \
    "${CODE_DIR}/00_build_run_manifest.py" "${CODE_DIR}/00_verify_md5.py" \
    "${CODE_DIR}/00_validate_umi_summary.py" "${CODE_DIR}/fusion_common.py" \
    "${CODE_DIR}/01_prepare_unique_sequences.py" "${CODE_DIR}/02_rotate_hpv_reference.py" \
    "${CODE_DIR}/02_select_hpv_candidates.py" \
    "${CODE_DIR}/03_call_and_aggregate_fusions.py" "${CODE_DIR}/04_prepare_technical_rescue.py" \
    "${CODE_DIR}/05_finalize_fusion_results.py" "${CODE_DIR}/98_validate_results.py" \
    "${CODE_DIR}/01_run_one_sample.sh"; do
    manifest_args+=(--script "$(basename "${script_file}")=${script_file}")
done
CURRENT_SIGNATURE="$(python3 "${CODE_DIR}/00_build_run_manifest.py" \
    "${manifest_args[@]}" --output "${CURRENT_MANIFEST}")"

PREVIOUS_SIGNATURE=""
if [[ -s "${RUN_MANIFEST}" ]]; then
    PREVIOUS_SIGNATURE="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["signature"])' "${RUN_MANIFEST}")"
fi
COMPLETE_SIGNATURE=""
if [[ -s "${COMPLETION_FILE}" ]]; then
    COMPLETE_SIGNATURE="$(sed -n 's/^signature=//p' "${COMPLETION_FILE}" | head -n 1)"
fi
if [[ "${FORCE}" != "1" && -n "${COMPLETE_SIGNATURE}" && \
      "${COMPLETE_SIGNATURE}" == "${CURRENT_SIGNATURE}" ]] && final_ready; then
    rm -f "${CURRENT_MANIFEST}"
    echo "[SKIP COMPLETE] ${SAMPLE}: validated results already exist and the run signature matches"
    echo "[FINAL] ${FINAL_DIR}"
    exit 0
fi
if [[ "${FORCE}" != "1" && -n "${PREVIOUS_SIGNATURE}" && \
      "${PREVIOUS_SIGNATURE}" != "${CURRENT_SIGNATURE}" ]]; then
    rm -f "${CURRENT_MANIFEST}"
    echo "[ERROR] Existing ${SAMPLE} intermediates were created from different inputs, code, references, or thresholds." >&2
    echo "[ERROR] Use a new sample name/result directory, or rerun intentionally with FORCE=1." >&2
    exit 1
fi
mv "${CURRENT_MANIFEST}" "${RUN_MANIFEST}"
rm -f "${COMPLETION_FILE}"

for command_name in bwa STAR samtools gzip; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "[ERROR] Required command not found: ${command_name}" >&2
        exit 1
    }
done
for bwa_suffix in amb ann bwt pac sa; do
    if [[ ! -s "${HPV_REF}.${bwa_suffix}" ]]; then
        echo "[ERROR] BWA index file missing: ${HPV_REF}.${bwa_suffix}" >&2
        exit 1
    fi
done
for star_file in Genome SA SAindex; do
    if [[ ! -s "${STAR_INDEX}/${star_file}" ]]; then
        echo "[ERROR] STAR index file missing: ${STAR_INDEX}/${star_file}" >&2
        exit 1
    fi
done

PIPELINE_START_EPOCH="$(date +%s)"
printf 'stage\tstatus\tduration_seconds\n' > "${TIMING_FILE}"
cat > "${RUN_CONFIG_FILE}" <<EOF
pipeline=HPV-human fusion batch supplied-merged run
pipeline_version=${PIPELINE_VERSION}
started_at=$(date '+%Y-%m-%dT%H:%M:%S%z')
config_file=${CONFIG_FILE}
sample=${SAMPLE}
merged_fastq=${MERGED_FASTQ}
r1_a=${R1_A}
r2_a=${R2_A}
r1_b=${R1_B}
r2_b=${R2_B}
hpv_ref=${HPV_REF}
star_index=${STAR_INDEX}
gene_bed=${GENE_BED}
hpv_type_map=${HPV_TYPE_MAP}
threads=${THREADS}
thresholds=${THRESHOLDS}
primary_threshold=${PRIMARY_THRESHOLD}
force=${FORCE}
require_umi=${REQUIRE_UMI}
EOF

echo "============================================================"
echo "[PIPELINE]     HPV-human fusion supplied-merged run ${PIPELINE_VERSION}"
echo "[SAMPLE]       ${SAMPLE}"
echo "[MERGED FASTQ] ${MERGED_FASTQ} (supplied input; no merge is performed)"
echo "[RAW A]        ${R1_A} | ${R2_A}"
echo "[RAW B]        ${R1_B:-.} | ${R2_B:-.}"
echo "[FINAL]        ${FINAL_DIR}"
echo "[REUSE POLICY] current formal result directory only"
echo "============================================================"

STEP_START="$(date +%s)"
INPUT_MARKER="${STATE_DIR}/01_input_validation.done"
raw_fastqs=("${R1_A}" "${R2_A}")
if [[ -n "${R1_B}" ]]; then
    raw_fastqs+=("${R1_B}" "${R2_B}")
fi
if stage_ready "${INPUT_MARKER}" "${MERGED_FASTQ}" "${raw_fastqs[@]}"; then
    STEP_STATUS="REUSE"
    echo "[SKIP 1/11] Input FASTQ validation already passed"
else
    STEP_STATUS="RUN"
    rm -f "${INPUT_MARKER}"
    echo "[STEP 1/11] Validate supplied merged FASTQ and corresponding raw FASTQ"
    gzip -t "${MERGED_FASTQ}"
    gzip -t "${raw_fastqs[@]}"
    if [[ "${VERIFY_RAW_MD5}" != "0" ]]; then
        md5_args=()
        missing_md5=0
        for raw_fastq in "${raw_fastqs[@]}"; do
            md5_file=""
            for candidate in "${raw_fastq}.md5" "${raw_fastq}.md5sum" "$(dirname "${raw_fastq}")/MD5.txt"; do
                if [[ -s "${candidate}" ]]; then
                    md5_file="${candidate}"
                    break
                fi
            done
            if [[ -n "${md5_file}" ]]; then
                md5_args+=(--pair "${raw_fastq}" "${md5_file}")
            else
                missing_md5=$((missing_md5 + 1))
                echo "[WARN] No MD5 sidecar/MD5.txt found for ${raw_fastq}; gzip integrity was checked"
            fi
        done
        if [[ "${VERIFY_RAW_MD5}" == "1" && "${missing_md5}" -gt 0 ]]; then
            echo "[ERROR] VERIFY_RAW_MD5=1 but ${missing_md5} raw FASTQ files have no MD5 entry" >&2
            exit 1
        fi
        if [[ "${#md5_args[@]}" -gt 0 ]]; then
            python3 "${CODE_DIR}/00_verify_md5.py" "${md5_args[@]}" \
                | tee "${LOG_DIR}/01_raw_md5.json"
        else
            printf '{"status":"UNAVAILABLE","reason":"no MD5 sidecars; gzip integrity PASS"}\n' \
                > "${LOG_DIR}/01_raw_md5.json"
        fi
    else
        printf '{"status":"DISABLED","reason":"VERIFY_RAW_MD5=0"}\n' \
            > "${LOG_DIR}/01_raw_md5.json"
    fi
    mark_complete "${INPUT_MARKER}"
fi
record_timing "01_input_validation" "${STEP_STATUS}" "${STEP_START}"

UNIQUE_FASTA="${WORK_DIR}/${SAMPLE}.unique_sequences.fa"
UNIQUE_METADATA="${WORK_DIR}/${SAMPLE}.unique_sequences.tsv.gz"
PREPARE_SUMMARY="${WORK_DIR}/${SAMPLE}.prepare_summary.json"
STEP_START="$(date +%s)"
PREPARE_MARKER="${STATE_DIR}/02_prepare_unique_sequences.done"
if stage_ready "${PREPARE_MARKER}" "${UNIQUE_FASTA}" "${UNIQUE_METADATA}" "${PREPARE_SUMMARY}"; then
    STEP_STATUS="REUSE"
    echo "[SKIP 2/11] Reuse exact unique sequences and UMI metadata from this formal run"
else
    STEP_STATUS="RUN"
    rm -f "${PREPARE_MARKER}"
    echo "[STEP 2/11] Extract paired-read UMI and collapse identical merged sequences"
    prepare_args=(
        --sample "${SAMPLE}" \
        --merged-fastq "${MERGED_FASTQ}" \
        --r1-a "${R1_A}" --r2-a "${R2_A}" \
        --out-dir "${WORK_DIR}" \
        --max-missing-raw-id-rate "${MAX_MISSING_RAW_ID_RATE}"
    )
    if [[ -n "${R1_B}" ]]; then
        prepare_args+=(--r1-b "${R1_B}" --r2-b "${R2_B}")
    fi
    python3 "${CODE_DIR}/01_prepare_unique_sequences.py" "${prepare_args[@]}" \
        2>&1 | tee "${LOG_DIR}/02_prepare_unique_sequences.log"
    gzip -t "${UNIQUE_METADATA}"
    mark_complete "${PREPARE_MARKER}"
fi
record_timing "02_prepare_unique_sequences" "${STEP_STATUS}" "${STEP_START}"

if [[ "${REQUIRE_UMI}" == "1" ]]; then
    echo "[CHECK 2/11] Require UMI-positive raw tube A, raw tube B and merged reads"
    python3 "${CODE_DIR}/00_validate_umi_summary.py" \
        --summary "${PREPARE_SUMMARY}" \
        --required-tube A --required-tube B \
        | tee "${LOG_DIR}/02_umi_requirement.json"
fi

HPV_ALIGNMENTS="${WORK_DIR}/${SAMPLE}.hpv_candidates.tsv"
STAR_QUERIES="${WORK_DIR}/${SAMPLE}.star_queries.fa"
HPV_SUMMARY="${WORK_DIR}/${SAMPLE}.hpv_candidate_summary.json"
ROTATED_HPV_REF=""

prepare_rotated_hpv_index() {
    [[ -z "${ROTATED_HPV_REF}" ]] || return 0
    ROTATED_HPV_DIR="$(mktemp -d "${WORK_DIR}/.hpv-circular.XXXXXXXX")"
    trap 'rm -rf -- "${ROTATED_HPV_DIR}"' EXIT
    ROTATED_HPV_REF="${ROTATED_HPV_DIR}/rotated.fa"
    python3 "${CODE_DIR}/02_rotate_hpv_reference.py" "${HPV_REF}" > "${ROTATED_HPV_REF}"
    bwa index "${ROTATED_HPV_REF}" > "${LOG_DIR}/hpv_circular_index.log" 2>&1
}

align_hpv() {
    local query_fasta="$1"
    bwa mem -Y -t "${THREADS}" "${ROTATED_HPV_REF}" "${query_fasta}" &&
        printf '@CO\tHPV_ORIGINAL_REFERENCE\n' &&
        bwa mem -Y -t "${THREADS}" "${HPV_REF}" "${query_fasta}"
}

STEP_START="$(date +%s)"
HPV_MARKER="${STATE_DIR}/03_hpv_candidate_alignment.done"
if stage_ready "${HPV_MARKER}" "${HPV_ALIGNMENTS}" "${HPV_SUMMARY}"; then
    STEP_STATUS="REUSE"
    echo "[SKIP 3/11] Reuse HPV candidate alignment from this formal run"
else
    STEP_STATUS="RUN"
    rm -f "${HPV_MARKER}"
    echo "[STEP 3/11] Align exact unique sequences to original and rotated HPV references"
    prepare_rotated_hpv_index
    align_hpv "${UNIQUE_FASTA}" \
        2> "${LOG_DIR}/03_bwa_hpv.stderr.log" \
        | python3 "${CODE_DIR}/02_select_hpv_candidates.py" \
            --sam - --rotated-first \
            --unique-fasta "${UNIQUE_FASTA}" \
            --out-alignments "${HPV_ALIGNMENTS}" \
            --out-star-fasta "${STAR_QUERIES}" \
            --summary "${HPV_SUMMARY}" \
            --prefilter-min-hpv "${PREFILTER_MIN_HPV}" \
            --prefilter-min-flank "${PREFILTER_MIN_FLANK}" \
        2>&1 | tee "${LOG_DIR}/03_select_hpv_candidates.log"
    mark_complete "${HPV_MARKER}"
fi
record_timing "03_hpv_candidate_alignment" "${STEP_STATUS}" "${STEP_START}"
NO_HPV_CANDIDATES=0
if [[ ! -s "${STAR_QUERIES}" ]]; then
    NO_HPV_CANDIDATES=1
    echo "[INFO] No HPV partial-alignment candidates; a valid zero-event result will be produced"
fi

STAR_PREFIX="${WORK_DIR}/star/${SAMPLE}."
STAR_BAM="${STAR_PREFIX}Aligned.out.bam"
STEP_START="$(date +%s)"
STAR_MARKER="${STATE_DIR}/04_human_star_alignment.done"
if [[ "${NO_HPV_CANDIDATES}" == "1" ]]; then
    STEP_STATUS="NO_INPUT"
    echo "[SKIP 4/11] No HPV candidates; human STAR alignment is unnecessary"
    mark_complete "${STAR_MARKER}"
elif stage_ready "${STAR_MARKER}" "${STAR_BAM}"; then
    STEP_STATUS="REUSE"
    echo "[SKIP 4/11] Reuse human STAR BAM from this formal run"
else
    STEP_STATUS="RUN"
    rm -f "${STAR_MARKER}"
    echo "[STEP 4/11] Align complete HPV candidates and non-HPV flanks to human"
    STAR \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${STAR_QUERIES}" \
        --runThreadN "${THREADS}" \
        --outTmpDir "${WORK_DIR}/star/${SAMPLE}.STARtmp.${PIPELINE_START_EPOCH}.$$" \
        --alignEndsType Local \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM HardClip \
        --chimMultimapNmax 10 \
        --outFilterMultimapNmax 10 \
        --outSAMmultNmax 10 \
        --outSAMmapqUnique 60 \
        --outSAMprimaryFlag OneBestScore \
        --outSAMattributes NH HI AS nM NM \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outFileNamePrefix "${STAR_PREFIX}" \
        > "${LOG_DIR}/04_star.stdout.log" \
        2> "${LOG_DIR}/04_star.stderr.log"
    samtools quickcheck -v "${STAR_BAM}"
    mark_complete "${STAR_MARKER}"
fi
record_timing "04_human_star_alignment" "${STEP_STATUS}" "${STEP_START}"

STRICT_EVIDENCE="${STRICT_CALL_DIR}/${SAMPLE}.fusion_evidence.tsv.gz"
STRICT_SUMMARY="${STRICT_CALL_DIR}/${SAMPLE}.fusion_summary.json"
STEP_START="$(date +%s)"
STRICT_MARKER="${STATE_DIR}/05_strict_evidence.done"
if stage_ready "${STRICT_MARKER}" "${STRICT_EVIDENCE}" "${STRICT_SUMMARY}"; then
    STEP_STATUS="REUSE"
    echo "[SKIP 5/11] Reuse complete strict evidence from this formal run"
else
    STEP_STATUS="RUN"
    rm -f "${STRICT_MARKER}"
    echo "[STEP 5/11] Build complete strict evidence audit"
    strict_args=(
            --sample "${SAMPLE}" \
            --hpv-alignments "${HPV_ALIGNMENTS}" \
            --sequence-metadata "${UNIQUE_METADATA}" \
            --hpv-type-map "${HPV_TYPE_MAP}" \
            --gene-bed "${GENE_BED}" \
            --out-dir "${STRICT_CALL_DIR}" \
            --thresholds "${THRESHOLDS}" \
            --min-human-mapq "${MIN_HUMAN_MAPQ}" \
            --min-hpv-mapq "${MIN_HPV_MAPQ}" \
            --max-junction-gap "${MAX_JUNCTION_GAP}" \
            --full-human-fraction "${FULL_HUMAN_FRACTION}" \
            --audit-only
    )
    if [[ "${NO_HPV_CANDIDATES}" == "1" ]]; then
        python3 "${CODE_DIR}/03_call_and_aggregate_fusions.py" \
            --human-sam /dev/null "${strict_args[@]}" \
            2>&1 | tee "${LOG_DIR}/05_strict_evidence.log"
    else
        samtools view -h "${STAR_BAM}" \
            | python3 "${CODE_DIR}/03_call_and_aggregate_fusions.py" \
                --human-sam - "${strict_args[@]}" \
            2>&1 | tee "${LOG_DIR}/05_strict_evidence.log"
    fi
    gzip -t "${STRICT_EVIDENCE}"
    mark_complete "${STRICT_MARKER}"
fi
record_timing "05_strict_evidence" "${STEP_STATUS}" "${STEP_START}"

RESCUE_FASTA="${RESCUE_DIR}/${SAMPLE}.technical_rescue.fa"
RESCUE_METADATA="${RESCUE_DIR}/${SAMPLE}.technical_rescue.tsv.gz"
RESCUE_MAP="${RESCUE_DIR}/${SAMPLE}.technical_rescue_map.tsv"
RESCUE_PREP_SUMMARY="${RESCUE_DIR}/${SAMPLE}.technical_rescue_summary.json"
RESCUE_HPV_ALIGNMENTS="${RESCUE_DIR}/${SAMPLE}.hpv_candidates.tsv"
RESCUE_STAR_QUERIES="${RESCUE_DIR}/${SAMPLE}.star_queries.fa"
RESCUE_STAR_PREFIX="${RESCUE_DIR}/star/${SAMPLE}."
RESCUE_STAR_BAM="${RESCUE_STAR_PREFIX}Aligned.out.bam"
RESCUE_EVIDENCE="${RESCUE_CALL_DIR}/${SAMPLE}.fusion_evidence.tsv.gz"

STEP_START="$(date +%s)"
RESCUE_PREP_MARKER="${STATE_DIR}/06_prepare_rescue.done"
if [[ "${ENABLE_TECHNICAL_RESCUE}" != "1" ]]; then
    STEP_STATUS="DISABLED"
    echo "[SKIP 6/11] Technical rescue disabled by configuration"
elif stage_ready "${RESCUE_PREP_MARKER}" "${RESCUE_METADATA}" "${RESCUE_MAP}" "${RESCUE_PREP_SUMMARY}"; then
    STEP_STATUS="REUSE"
    echo "[SKIP 6/11] Reuse normalized rescue set from this formal run"
else
    STEP_STATUS="RUN"
    rm -f "${RESCUE_PREP_MARKER}"
    echo "[STEP 6/11] Trim terminal technical cassettes and collapse identical biological sequences"
    python3 "${CODE_DIR}/04_prepare_technical_rescue.py" \
        --strict-evidence "${STRICT_EVIDENCE}" \
        --sequence-metadata "${UNIQUE_METADATA}" \
        --out-fasta "${RESCUE_FASTA}" \
        --out-metadata "${RESCUE_METADATA}" \
        --out-map "${RESCUE_MAP}" \
        --summary "${RESCUE_PREP_SUMMARY}" \
        --threshold "${PRIMARY_THRESHOLD}" \
        --min-human-mapq "${MIN_HUMAN_MAPQ}" \
        --min-hpv-mapq "${MIN_HPV_MAPQ}" \
        --terminal-slack "${RESCUE_TERMINAL_SLACK}" \
        --min-trimmed-length "${RESCUE_MIN_TRIMMED_LENGTH}" \
        2>&1 | tee "${LOG_DIR}/06_prepare_technical_rescue.log"
    gzip -t "${RESCUE_METADATA}"
    mark_complete "${RESCUE_PREP_MARKER}"
fi
record_timing "06_prepare_normalized_rescue" "${STEP_STATUS}" "${STEP_START}"

RESCUE_COUNT=0
if [[ "${ENABLE_TECHNICAL_RESCUE}" == "1" && -s "${RESCUE_PREP_SUMMARY}" ]]; then
    RESCUE_COUNT="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["rescue_sequences_written"])' "${RESCUE_PREP_SUMMARY}")"
fi

if [[ "${RESCUE_COUNT}" -gt 0 ]]; then
    STEP_START="$(date +%s)"
    RESCUE_HPV_MARKER="${STATE_DIR}/07_rescue_hpv_alignment.done"
    if stage_ready "${RESCUE_HPV_MARKER}" "${RESCUE_HPV_ALIGNMENTS}" "${RESCUE_DIR}/${SAMPLE}.hpv_candidate_summary.json"; then
        STEP_STATUS="REUSE"
        echo "[SKIP 7/11] Reuse rescue HPV alignment from this formal run"
    else
        STEP_STATUS="RUN"
        rm -f "${RESCUE_HPV_MARKER}"
        echo "[STEP 7/11] Re-align ${RESCUE_COUNT} unique trimmed rescue sequences to HPV"
        prepare_rotated_hpv_index
        align_hpv "${RESCUE_FASTA}" \
            2> "${LOG_DIR}/07_rescue_bwa_hpv.stderr.log" \
            | python3 "${CODE_DIR}/02_select_hpv_candidates.py" \
                --sam - --rotated-first \
                --unique-fasta "${RESCUE_FASTA}" \
                --out-alignments "${RESCUE_HPV_ALIGNMENTS}" \
                --out-star-fasta "${RESCUE_STAR_QUERIES}" \
                --summary "${RESCUE_DIR}/${SAMPLE}.hpv_candidate_summary.json" \
                --prefilter-min-hpv "${PREFILTER_MIN_HPV}" \
                --prefilter-min-flank "${PREFILTER_MIN_FLANK}" \
            2>&1 | tee "${LOG_DIR}/07_select_rescue_hpv_candidates.log"
        mark_complete "${RESCUE_HPV_MARKER}"
    fi
    record_timing "07_rescue_hpv_alignment" "${STEP_STATUS}" "${STEP_START}"

    if [[ ! -s "${RESCUE_STAR_QUERIES}" ]]; then
        echo "[SKIP 8/11] Trimmed rescue sequences produced no HPV candidates"
        echo "[SKIP 9/11] Trimmed rescue sequences produced no HPV candidates"
        record_timing "08_rescue_human_alignment" "NO_INPUT" "$(date +%s)"
        record_timing "09_rescue_breakpoint_call" "NO_INPUT" "$(date +%s)"
        RESCUE_COUNT=0
    else
        STEP_START="$(date +%s)"
        RESCUE_STAR_MARKER="${STATE_DIR}/08_rescue_human_alignment.done"
        if stage_ready "${RESCUE_STAR_MARKER}" "${RESCUE_STAR_BAM}"; then
            STEP_STATUS="REUSE"
            echo "[SKIP 8/11] Reuse rescue human STAR BAM from this formal run"
        else
            STEP_STATUS="RUN"
            rm -f "${RESCUE_STAR_MARKER}"
            echo "[STEP 8/11] Align trimmed rescue candidates and flanks to human"
            STAR \
                --genomeDir "${STAR_INDEX}" \
                --readFilesIn "${RESCUE_STAR_QUERIES}" \
                --runThreadN "${THREADS}" \
                --outTmpDir "${RESCUE_DIR}/star/${SAMPLE}.STARtmp.${PIPELINE_START_EPOCH}.$$" \
                --alignEndsType Local \
                --chimSegmentMin 10 \
                --chimOutType WithinBAM HardClip \
                --chimMultimapNmax 10 \
                --outFilterMultimapNmax 10 \
                --outSAMmultNmax 10 \
                --outSAMmapqUnique 60 \
                --outSAMprimaryFlag OneBestScore \
                --outSAMattributes NH HI AS nM NM \
                --outSAMtype BAM Unsorted \
                --outSAMunmapped Within \
                --outFileNamePrefix "${RESCUE_STAR_PREFIX}" \
                > "${LOG_DIR}/08_rescue_star.stdout.log" \
                2> "${LOG_DIR}/08_rescue_star.stderr.log"
            samtools quickcheck -v "${RESCUE_STAR_BAM}"
            mark_complete "${RESCUE_STAR_MARKER}"
        fi
        record_timing "08_rescue_human_alignment" "${STEP_STATUS}" "${STEP_START}"

        STEP_START="$(date +%s)"
        RESCUE_CALL_MARKER="${STATE_DIR}/09_rescue_breakpoint_call.done"
        if stage_ready "${RESCUE_CALL_MARKER}" "${RESCUE_EVIDENCE}" "${RESCUE_CALL_DIR}/${SAMPLE}.fusion_summary.json"; then
            STEP_STATUS="REUSE"
            echo "[SKIP 9/11] Reuse rescue breakpoint evidence from this formal run"
        else
            STEP_STATUS="RUN"
            rm -f "${RESCUE_CALL_MARKER}"
            echo "[STEP 9/11] Recalculate breakpoints from physically trimmed rescue sequences"
            samtools view -h "${RESCUE_STAR_BAM}" \
                | python3 "${CODE_DIR}/03_call_and_aggregate_fusions.py" \
                    --sample "${SAMPLE}" \
                    --human-sam - \
                    --hpv-alignments "${RESCUE_HPV_ALIGNMENTS}" \
                    --sequence-metadata "${RESCUE_METADATA}" \
                    --hpv-type-map "${HPV_TYPE_MAP}" \
                    --gene-bed "${GENE_BED}" \
                    --out-dir "${RESCUE_CALL_DIR}" \
                    --thresholds "${THRESHOLDS}" \
                    --min-human-mapq "${MIN_HUMAN_MAPQ}" \
                    --min-hpv-mapq "${MIN_HPV_MAPQ}" \
                    --max-junction-gap "${MAX_JUNCTION_GAP}" \
                    --full-human-fraction "${FULL_HUMAN_FRACTION}" \
                    --audit-only \
                2>&1 | tee "${LOG_DIR}/09_call_rescue_evidence.log"
            gzip -t "${RESCUE_EVIDENCE}"
            mark_complete "${RESCUE_CALL_MARKER}"
        fi
        record_timing "09_rescue_breakpoint_call" "${STEP_STATUS}" "${STEP_START}"
    fi
else
    echo "[SKIP 7/11] No rescue sequences"
    echo "[SKIP 8/11] No rescue sequences"
    echo "[SKIP 9/11] No rescue sequences"
    record_timing "07_rescue_hpv_alignment" "NO_INPUT" "$(date +%s)"
    record_timing "08_rescue_human_alignment" "NO_INPUT" "$(date +%s)"
    record_timing "09_rescue_breakpoint_call" "NO_INPUT" "$(date +%s)"
fi

STEP_START="$(date +%s)"
FINAL_MARKER="${STATE_DIR}/10_finalize_results.done"
if stage_ready "${FINAL_MARKER}" \
    "${FINAL_EVENTS}" "${FINAL_SUMMARY}" \
    "${FINAL_DIR}/${SAMPLE}.fusion_events.main.tsv" \
    "${FINAL_DIR}/${SAMPLE}.fusion_events.candidate.tsv" \
    "${FINAL_DIR}/${SAMPLE}.fusion_events.review.tsv" \
    "${FINAL_DIR}/${SAMPLE}.fusion_loci.final.tsv" \
    "${FINAL_DIR}/${SAMPLE}.fusion_loci.main.tsv" \
    "${FINAL_DIR}/${SAMPLE}.fusion_loci.candidate.tsv" \
    "${FINAL_DIR}/${SAMPLE}.fusion_loci.review.tsv" \
    "${FINAL_DIR}/${SAMPLE}.fusion_evidence.final.tsv.gz" \
    "${FINAL_DIR}/${SAMPLE}.fusion_unique_sequences.final.tsv.gz" \
    "${FINAL_DIR}/${SAMPLE}.fusion_sequence_members.final.tsv.gz" \
    "${FINAL_DIR}/${SAMPLE}.fusion_umi_families.exact.tsv.gz" \
    "${FINAL_DIR}/${SAMPLE}.fusion_umi_families.final.tsv.gz" \
    "${FINAL_DIR}/${SAMPLE}.fusion_umi_families.corrected.tsv.gz" \
    "${FINAL_DIR}/${SAMPLE}.fusion_filter_summary.tsv"; then
    STEP_STATUS="REUSE"
    echo "[SKIP 10/11] Reuse compact final results from this formal run"
else
    STEP_STATUS="RUN"
    rm -f "${FINAL_MARKER}"
    echo "[STEP 10/11] Build compact event, locus, sequence, read and UMI tables"
    finalize_args=(
        --sample "${SAMPLE}"
        --strict-evidence "${STRICT_EVIDENCE}"
        --sequence-metadata "${UNIQUE_METADATA}"
        --out-dir "${FINAL_DIR}"
        --thresholds "${THRESHOLDS}"
        --primary-threshold "${PRIMARY_THRESHOLD}"
        --min-human-mapq "${MIN_HUMAN_MAPQ}"
        --min-hpv-mapq "${MIN_HPV_MAPQ}"
        --main-min-raw-reads "${MAIN_MIN_RAW_READS}"
        --main-min-umi-families "${MAIN_MIN_UMI_FAMILIES}"
        --candidate-min-raw-reads "${CANDIDATE_MIN_RAW_READS}"
        --candidate-min-umi-families "${CANDIDATE_MIN_UMI_FAMILIES}"
        --locus-human-window "${LOCUS_HUMAN_WINDOW}"
        --locus-hpv-window "${LOCUS_HPV_WINDOW}"
        --umi-correction-distance "${UMI_CORRECTION_DISTANCE}"
        --pipeline-version "${PIPELINE_VERSION}"
    )
    if [[ "${RESCUE_COUNT}" -gt 0 && -s "${RESCUE_EVIDENCE}" ]]; then
        finalize_args+=(--rescue-evidence "${RESCUE_EVIDENCE}" --rescue-map "${RESCUE_MAP}")
    fi
    python3 "${CODE_DIR}/05_finalize_fusion_results.py" "${finalize_args[@]}" \
        2>&1 | tee "${LOG_DIR}/10_finalize.log"
    mark_complete "${FINAL_MARKER}"
fi
record_timing "10_finalize_results" "${STEP_STATUS}" "${STEP_START}"

echo "[STEP 11/11] Independently validate IDs, coordinates, read/UMI aggregation and tiers"
STEP_START="$(date +%s)"
python3 "${CODE_DIR}/98_validate_results.py" \
    --sample "${SAMPLE}" \
    --final-dir "${FINAL_DIR}" \
    --primary-threshold "${PRIMARY_THRESHOLD}" \
    --main-min-raw-reads "${MAIN_MIN_RAW_READS}" \
    --main-min-umi-families "${MAIN_MIN_UMI_FAMILIES}" \
    --candidate-min-raw-reads "${CANDIDATE_MIN_RAW_READS}" \
    --candidate-min-umi-families "${CANDIDATE_MIN_UMI_FAMILIES}" \
    2>&1 | tee "${LOG_DIR}/11_validate_results.log"
record_timing "11_validate_results" "RUN" "${STEP_START}"
printf 'TOTAL\tRUN\t%d\n' "$(( $(date +%s) - PIPELINE_START_EPOCH ))" >> "${TIMING_FILE}"
COMPLETION_TMP="${COMPLETION_FILE}.tmp.$$"
{
    printf 'status=PASS\n'
    printf 'signature=%s\n' "${CURRENT_SIGNATURE}"
    printf 'completed_at=%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
} > "${COMPLETION_TMP}"
mv "${COMPLETION_TMP}" "${COMPLETION_FILE}"

echo "[DONE] ${SAMPLE} supplied-merged run completed and validated"
echo "[MAIN]      ${FINAL_DIR}/${SAMPLE}.fusion_events.main.tsv"
echo "[CANDIDATE] ${FINAL_DIR}/${SAMPLE}.fusion_events.candidate.tsv"
echo "[REVIEW]    ${FINAL_DIR}/${SAMPLE}.fusion_events.review.tsv"
echo "[ALL 30BP]  ${FINAL_EVENTS}"
echo "[LOCI]      ${FINAL_DIR}/${SAMPLE}.fusion_loci.final.tsv"
echo "[SEQUENCES] ${FINAL_DIR}/${SAMPLE}.fusion_unique_sequences.final.tsv.gz"
echo "[UMI EXACT] ${FINAL_DIR}/${SAMPLE}.fusion_umi_families.exact.tsv.gz"
echo "[UMI CORR.] ${FINAL_DIR}/${SAMPLE}.fusion_umi_families.corrected.tsv.gz"
echo "[SUMMARY]   ${FINAL_SUMMARY}"
echo "[TIMING]    ${TIMING_FILE}"
