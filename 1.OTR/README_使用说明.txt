On_target_rates 通用流程（2026-09-01）
=====================================

一、这次默认分析的数据
----------------------
原始数据根目录：
/data2/lindf/1_rawdata_fastq/HPV_fusion_lib/Targeted_Library_Rawdata/Targeted_Library_Tissue/Lib_Expression

只处理 samples.tsv 中 enabled=1 的8个生物学样本：
Tissue1, Tissue15, Tissue26, Tissue41, Tissue50, Tissue80, Tissue81, TissueC2@2

每个样本均明确指定 A 管和 B 管，因此不会自动把以下历史目录带入：
20250814_99plex_sample1_tisuue
20250926_133plex_sample18_tissue
20250926_133plex_sample19_tissue

Primer：20250520_126_plex.xlsx
（已经去掉原先重复的 HPV16-3(+)；代码运行前还会再次自动检查重复。）


二、服务器上的放置位置
----------------------
把压缩包内所有文件解压到：
/data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/On_target_rates/code

结果会写到：
/data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/On_target_rates/result

推荐结构：
On_target_rates/
├── code/
│   ├── 00_run_all.sh
│   ├── 01_preprocess_expression.sh
│   ├── 02_calculate_on_target_rates.py
│   ├── config.sh
│   ├── samples.tsv
│   └── 20250520_126_plex.xlsx
└── result/


三、运行方法
------------
cd /data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/On_target_rates/code

chmod +x 00_run_all.sh 01_preprocess_expression.sh 02_calculate_on_target_rates.py

nohup bash 00_run_all.sh config.sh > run_on_target_rates.nohup.log 2>&1 &

echo $!

tail -f run_on_target_rates.nohup.log

注意：00_run_all.sh 自己也会在 result/logs/ 中保存带时间戳的总日志。


四、完整流程
------------
每个生物学样本：

A管 R1/R2 --PEAR--> assembled + unassembled.forward
B管 R1/R2 --PEAR--> assembled + unassembled.forward
                       ↓
                    A+B合并
                       ↓
                    fastp过滤
                       ↓
          iFP / iFP_spacer / CTTGCAC统计
                       ↓
                 on-target rates

PEAR和fastp参数沿用旧脚本。
fastp明确使用 --disable_adapter_trimming，不剪掉 CTTGCAC。


五、主要结果目录
----------------
result/01_pear/
  每个样本A/B管分别PEAR后的结果和日志。

result/02_merge_AB/
  Tissue1_A_B_merged.fastq.gz 等。

result/03_fastp/
  Tissue1_qc1.fastq.gz 等，以及fastp html/json/log。

result/04_on_target_rates/
  primer_search_list.csv
      126个target拆成126个iFP + 126个iFP_spacer，共252条搜索序列。

  primer_audit.txt
      primer数量及重复检查结果。

  sum_onTargetRate_ALL_samples_long.csv
      最推荐用于后续R分析的长表。
      包含：sample_id, gene, number, target_type, seq,
            count, count_anchor, total_reads,
            on_target_rate_pct,
            on_target_anchor_rate_pct,
            anchor_within_target_pct。

  sum_onTargetRate_ALL_samples_matrix.csv
      与旧代码思路接近的矩阵表。
      每个样本单元格为 count|count_anchor；最后一行是 TOTAL_READS。

  overall_onTargetRate_ALL_samples.csv
      每个样本整体的“任意iFP/任意iFP_spacer”union统计，避免把各target直接相加造成重复计数。

  per_sample/<sample_id>/
      每个样本的详细结果、legacy风格lib文件和analysis_manifest.json。


六、count / count_anchor 的含义
-------------------------------
count：
某条read只要包含当前 iFP 或 iFP_spacer 序列，就对该target计1次。
同一read里即使该序列出现多次，也只计1次。

count_anchor：
上述read中，同时包含 CTTGCAC 的read数。

因此：
on_target_rate_pct = count / total_reads × 100
on_target_anchor_rate_pct = count_anchor / total_reads × 100
anchor_within_target_pct = count_anchor / count × 100

注意：iFP_spacer通常包含iFP，因此同一个target的iFP count和iFP_spacer count不能直接相加。


七、以后有新数据时怎么改
------------------------
核心程序原则上不需要修改。

1. 修改 config.sh：
   RAW_DIR      新原始数据根目录
   CODE_DIR     新项目code目录
   RESULT_DIR   新项目result目录
   PRIMER_EXCEL 新primer表路径

2. 修改 samples.tsv：
   每行写：
   sample_id    A_dir    B_dir    enabled

   例如：
   Tissue100    Sample_xxxTissue100A    Sample_xxxTissue100B    1

3. 运行：
   bash 00_run_all.sh config.sh

FASTQ文件名不要求必须带 _combined；程序只要求每个A/B目录中各自“正好有一个 *_R1.fastq.gz 和一个 *_R2.fastq.gz”。
因此既支持：
xxx_combined_R1.fastq.gz
也支持：
lib1_R1.fastq.gz


八、旧结果安全机制
------------------
旧Python脚本只要发现结果文件存在就skip，因此从127 plex改成126 plex后可能错误复用旧结果。
新版为每个样本保存 analysis_manifest.json，记录：
- primer Excel SHA256
- anchor
- FASTQ路径/大小/mtime
- 算法版本

只有这些全部一致才会skip。
只要更换primer Excel（例如127 plex -> 126 plex），会自动重新计算，不会静默复用旧结果。

如明确想强制重跑，可在config.sh设置：
FORCE_RERUN_PEAR=1
FORCE_RERUN_MERGE=1
FORCE_RERUN_FASTP=1
FORCE_RERUN_COUNT=1

通常不需要全部设为1；例如只想重新统计126 plex，只设 FORCE_RERUN_COUNT=1 即可。


============================================================
V2 更新：存在结果即跳过 + Bash bug 修复
============================================================
1. 修复 00_run_all.sh 中 COUNT_FORCE_ARG 空数组在 set -u 下导致的：
   COUNT_FORCE_ARG[@]: 为绑定变量

2. 默认 FORCE_RERUN_*=0 时：
   - PEAR：最终 *_assembled_plus_unassembled_R1.fastq.gz 存在且非空 -> 跳过
   - A+B：*_A_B_merged.fastq.gz 存在且非空 -> 跳过
   - fastp：*_qc1.fastq.gz 存在且非空 -> 跳过
   - target count：每个样本 sum_onTargetRate / lib / overall 三个结果均存在且非空 -> 跳过

3. 因此流程可以安全重复执行。已经完成的步骤不会重新计算，缺失的步骤会继续补跑。

4. 如果你修改了 FASTQ、primer Excel 或 anchor，并且明确希望覆盖旧结果，请在 config.sh 中把对应项设为1：
   FORCE_RERUN_PEAR=1
   FORCE_RERUN_MERGE=1
   FORCE_RERUN_FASTP=1
   FORCE_RERUN_COUNT=1
   正常情况下全部保持0。
