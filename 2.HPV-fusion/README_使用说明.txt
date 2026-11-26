HPV-human融合批量分析 V4.2.0-targeted-UMI-AB（2026-09-17）
============================================================

一、适用数据和本版输入
----------------------
本包专用于带H15-X7 UMI的组织靶向建库数据，并按原Tissue1正式V4的方式同时使用A、B两管：

1. A管原始R1/R2：提取A管UMI并建立read ID索引。
2. B管原始R1/R2：提取B管UMI并建立read ID索引。
3. 同一样本已经merge好的A+B单端FASTQ：按read ID关联回A/B管UMI。
4. 不执行merge；相同merged序列只比对一次，但完整保留raw read数、UMI家族数及A/B来源。
5. UMI家族定义含tube字段，所以A、B两管即使出现相同UMI也不会错误合并。

本包预填8个样本：Tissue1、Tissue15、Tissue26、Tissue41、Tissue50、Tissue80、Tissue81、TissueC2@2。

二、安装目录
------------
将压缩包内容解压到：

/data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/HPV_fusion/20260917_HPV_fusion_code_analysis

解压后应为：

20260917_HPV_fusion_code_analysis/
├── code/
└── result/

三、sample.list
---------------
文件：code/sample.list。必须用Tab分隔，六列为：

sample  merged_fastq  raw_r1_a  raw_r2_a  raw_r1_b  raw_r2_b

本版要求A、B两管四个原始FASTQ全部填写，且merged_fastq必须对应同一样本A+B的数据。所有路径均使用绝对路径；样本名不能重复；以#开头的行不运行。

四、运行
--------
cd /data2/lindf/project/HPV_project/20260831_国妇婴组织样本_靶向建库_数据分析/HPV_fusion/20260917_HPV_fusion_code_analysis/code

先做本地逻辑自检：

python3 99_self_test.py
python3 99_batch_self_test.py

后台批量运行：

nohup bash 00_run_batch.sh > ../result/run_batch.nohup.log 2>&1 &

查看日志：

tail -f ../result/run_batch.nohup.log

默认PARALLEL_JOBS=8、THREADS_PER_SAMPLE=2，即8个样本并行、每个样本2个比对线程。若RNA-seq的8任务批次仍在运行，建议等它结束后再启动本批，避免两个批次同时大量读取FASTQ和写中间文件造成磁盘I/O竞争。

五、运行前保护
--------------
1. merged与A/B四个raw FASTQ均执行gzip完整性检查。
2. VERIFY_RAW_MD5=1：四个raw FASTQ必须存在并通过.md5校验。
3. REQUIRE_UMI=1：步骤2结束后，A管、B管及merged关联结果都必须检出可用UMI，否则样本在BWA/STAR比对前失败。
4. merged read ID关联raw A+B的缺失率不得超过MAX_MISSING_RAW_ID_RATE（默认0.001）。
5. 运行签名记录输入、参考、代码和阈值；同名旧结果不匹配时拒绝混用。

六、断点续跑与存在即跳过
--------------------------
重复执行00_run_batch.sh时：

1. 最终文件齐全、独立校验PASS且运行签名一致：显示[SKIP COMPLETE]。
2. 同签名未完成样本：按logs/state中的步骤标记续跑。
3. 输入、代码、参考或阈值变化：报签名不一致，不会误用旧中间文件。

确需覆盖重跑时可用FORCE=1，但它会重跑传入样本表中的全部样本：

FORCE=1 bash 00_run_batch.sh /absolute/path/one_sample.list

七、结果位置
------------
每个样本写入result/<sample>/：

work/    唯一序列、UMI元数据、HPV候选和STAR中间文件
audit/   严格证据及技术序列救援审计文件
final/   正式结果
logs/    日志、耗时、状态、MD5及运行签名

常用正式文件：

<sample>.fusion_events.main.tsv               主结果事件
<sample>.fusion_events.candidate.tsv          候选事件
<sample>.fusion_events.review.tsv             低支持/需复核事件
<sample>.fusion_events.final.tsv              30 bp门槛下全部事件
<sample>.fusion_loci.final.tsv                邻近断点聚合后的位点
<sample>.fusion_unique_sequences.final.tsv.gz 去冗余序列及read/UMI计数
<sample>.fusion_umi_families.exact.tsv.gz     精确UMI家族
<sample>.fusion_umi_families.corrected.tsv.gz 纠错后UMI家族（辅助）
<sample>.fusion_summary.final.json            结果汇总

批次状态位于result/_batch_logs/<日期时间_PID>/batch_summary.tsv。单样本详细日志位于result/<sample>/logs/<sample>.pipeline.log；耗时位于<sample>.runtime_metrics.tsv。

八、主要分析顺序
----------------
1. 检查FASTQ、gzip与raw MD5。
2. 从raw A/B提取UMI，按read ID关联merged，并合并完全相同序列。
3. 强制检查A、B及merged确实具有可用UMI。
4. 唯一序列比对HPV，筛选部分比对候选。
5. 候选序列比对人基因组，计算HPV/人断点及30/40/50 bp证据。
6. 对疑似技术序列干扰的候选做裁剪救援。
7. 按事件、位点、唯一序列、raw reads和UMI家族聚合并分级。
8. 独立校验所有正式输出，PASS后才写完成标志。

九、失败排查
------------
1. 看result/_batch_logs/最近一次/batch_summary.tsv。
2. 看同目录parallel.joblog中Exitval非0的样本。
3. 看result/<sample>/logs/<sample>.pipeline.log。
4. UMI问题重点看work/<sample>.prepare_summary.json和logs/02_umi_requirement.json。
5. 修复后直接重跑；同签名样本会断点续跑。
