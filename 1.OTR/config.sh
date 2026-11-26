#!/usr/bin/env bash
# ============================================================
# On-target rates 通用流程配置文件
# 以后有新批次时，原则上只需要修改本文件 + samples.tsv
# ============================================================

# 1) 原始 FASTQ 根目录
RAW_DIR="/data2/lindf/1_rawdata_fastq/HPV_fusion_lib/Targeted_Library_Rawdata/Targeted_Library_Tissue/Lib_Expression"

# 2) 本批代码目录（config、samples.tsv、primer Excel、脚本均放这里）
CODE_DIR="/data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/On_target_rates/code"

# 3) 本批结果目录
RESULT_DIR="/data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/On_target_rates/result"

# 4) Primer 表与样本表
PRIMER_EXCEL="${CODE_DIR}/20250520_126_plex.xlsx"
SAMPLE_TABLE="${CODE_DIR}/samples.tsv"
PRIMER_SHEET="Sheet1"

# 5) ARRP/X7 参数
ANCHOR_RC="CTTGCAC"

# 6) 并行参数
# 总 CPU 需求大约 = PARALLEL_JOBS × THREADS（PEAR/fastp阶段）
PARALLEL_JOBS=5
THREADS=20

# 7) PEAR 参数（沿用原脚本）
PEAR_MIN_OVERLAP=30
PEAR_MAX_ASSEMBLY_LENGTH=604
PEAR_MIN_ASSEMBLY_LENGTH=30
PEAR_QUALITY_THRESHOLD=20
PEAR_TRIM_LENGTH=30
PEAR_P_VALUE=0.01

# 8) fastp 参数（沿用原脚本；明确禁用 adapter trimming）
FASTP_QUALIFIED_QUALITY_PHRED=20
FASTP_UNQUALIFIED_PERCENT_LIMIT=20
FASTP_LENGTH_REQUIRED=30
FASTP_LENGTH_LIMIT=604

# 9) 是否强制重跑：0=已有且有效则跳过；1=强制重跑
FORCE_RERUN_PEAR=0
FORCE_RERUN_MERGE=0
FORCE_RERUN_FASTP=0
FORCE_RERUN_COUNT=0
