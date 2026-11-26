#!/usr/bin/env bash
# HPV-human fusion targeted-library batch configuration (V4.2.0-targeted-UMI-AB).
# Analysis inputs are read from sample.list; values here may be overridden at launch.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-/data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/HPV_fusion/20260917_HPV_fusion_code_analysis}"
CODE_DIR="${CODE_DIR:-${PIPELINE_ROOT}/code}"
RESULT_DIR="${RESULT_DIR:-${PIPELINE_ROOT}/result}"
SAMPLE_LIST="${SAMPLE_LIST:-${CODE_DIR}/sample.list}"

# Eight samples run simultaneously. Two alignment threads per sample gives a
# default ceiling of about 16 alignment threads. Adjust only for the server's
# actual CPU and RAM capacity.
PARALLEL_JOBS="${PARALLEL_JOBS:-8}"
THREADS_PER_SAMPLE="${THREADS_PER_SAMPLE:-2}"

HPV_REF="${HPV_REF:-/data/lindf/4_ref_database/HPV_ref/hpv_allTypes_maskedByTRF.fasta}"
STAR_INDEX="${STAR_INDEX:-/data/lindf/4_ref_database/star_ref/genome_dna_fa_index}"
GENE_BED="${GENE_BED:-/data/public/SequencingData/xiuxh/ENSEMBL_GTFandFASTA/genes_position_chr.tab}"
HPV_TYPE_MAP="${HPV_TYPE_MAP:-${CODE_DIR}/hpv_GenBankID_fromIHPRC.csv}"

# auto: validate an MD5 file when found; always run gzip integrity checks.
# 1: require an MD5 entry for every raw FASTQ. 0: gzip integrity checks only.
VERIFY_RAW_MD5="${VERIFY_RAW_MD5:-1}"

# This release is dedicated to targeted libraries carrying H15-X7 UMIs.
# Both raw tube A and tube B must yield UMI-eligible reads before alignment.
REQUIRE_UMI="${REQUIRE_UMI:-1}"

PREFILTER_MIN_HPV="${PREFILTER_MIN_HPV:-15}"
PREFILTER_MIN_FLANK="${PREFILTER_MIN_FLANK:-15}"
THRESHOLDS="${THRESHOLDS:-20,30,40,50}"
PRIMARY_THRESHOLD="${PRIMARY_THRESHOLD:-30}"
MIN_HUMAN_MAPQ="${MIN_HUMAN_MAPQ:-20}"
MIN_HPV_MAPQ="${MIN_HPV_MAPQ:-0}"
MAX_JUNCTION_GAP="${MAX_JUNCTION_GAP:-20}"
FULL_HUMAN_FRACTION="${FULL_HUMAN_FRACTION:-0.98}"
MAX_MISSING_RAW_ID_RATE="${MAX_MISSING_RAW_ID_RATE:-0.001}"

MAIN_MIN_RAW_READS="${MAIN_MIN_RAW_READS:-3}"
MAIN_MIN_UMI_FAMILIES="${MAIN_MIN_UMI_FAMILIES:-2}"
CANDIDATE_MIN_RAW_READS="${CANDIDATE_MIN_RAW_READS:-2}"
CANDIDATE_MIN_UMI_FAMILIES="${CANDIDATE_MIN_UMI_FAMILIES:-1}"

ENABLE_TECHNICAL_RESCUE="${ENABLE_TECHNICAL_RESCUE:-1}"
RESCUE_TERMINAL_SLACK="${RESCUE_TERMINAL_SLACK:-10}"
RESCUE_MIN_TRIMMED_LENGTH="${RESCUE_MIN_TRIMMED_LENGTH:-60}"
LOCUS_HUMAN_WINDOW="${LOCUS_HUMAN_WINDOW:-5}"
LOCUS_HPV_WINDOW="${LOCUS_HPV_WINDOW:-5}"
UMI_CORRECTION_DISTANCE="${UMI_CORRECTION_DISTANCE:-1}"

# 0: resume matching incomplete samples and skip matching PASS samples.
# 1: rerun every stage for listed samples and replace their outputs in place.
FORCE="${FORCE:-0}"
