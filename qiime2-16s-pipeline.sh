#!/bin/bash
##########################################################
# QIIME2 16S 扩增子分析流程 (通用模板)
# 适用: 16S 双端测序数据
# 环境: qc_preprocess (质控/剪接) + qiime2-2025.7 (核心分析)
# 数据库: SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza
##########################################################
# 首次使用: 先运行 qiime2-16s-pipeline_install.sh 安装环境,
#           再运行本脚本分析数据
##########################################################

## 0. 参数设置 (用户根据实际修改)
wd=~/amplicon_analysis
metadata=metadata.txt
# classifiers/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza 需预先下载至 classifiers/ 目录
# 下载: https://pan.baidu.com/s/1yFOvSJXofc6H7AahHSJq2A?pwd=5gdz
#     或从 SILVA 官网获取: https://www.arb-silva.de/
# 注: SILVA138.2 分类器仅兼容 qiime2 ≥ 2023 版本 (本流程以 qiime2-2025.7 为例)
#     若使用 qiime2-2023.2, 请改用 silva-138-99-nb-classifier.qza

## 1. 目录初始化与文件放置
# 运行本脚本前, 请将以下文件放到指定位置:
#   - 原始测序数据: *.fq.gz 放入 seq/ 目录 (命名格式: 样本名_1.fq.gz / 样本名_2.fq.gz)
#   - 元数据文件: metadata.txt 放在工作目录下
#     格式: TSV (制表符分隔), 第一行为表头, 第一列为样本ID, 必须含 Group 分组列
#     示例 (参考: examples/metadata.txt):
#       SampleID  Group  ...
#       WT1       WT     ...
#       KO1       KO     ...
mkdir -p ${wd}/{seq,trimmed,qiime2,logs,results/{fastqc_raw,fastqc_trimmed,cutadapt_logs,export}}
cd ${wd}

# ln /path/to/raw/seq/*.fq.gz seq/
# ln /path/to/metadata.txt ./

## 2. 环境验证
conda activate qc_preprocess
fastqc --version
cutadapt --version
multiqc --version

conda activate qiime2-2025.7
qiime --version

## 3. 原始数据质控
conda activate qc_preprocess
cd ${wd}

fastqc seq/*_1.fq.gz seq/*_2.fq.gz -t 8 -o results/fastqc_raw
multiqc results/fastqc_raw/ -o results/fastqc_raw/ -n multiqc_report_raw.html
# 检查: results/fastqc_raw/multiqc_report_raw.html

## 4. 引物切除 
conda activate qc_preprocess
cd ${wd}

for f in seq/*_1.fq.gz; do
    base=$(basename "${f}" _1.fq.gz)
    r="seq/${base}_2.fq.gz"
    [ ! -f "${r}" ] && echo "警告: 反向文件缺失 ${r}" && continue

    cutadapt \
        -g GTGCCAGCMGCCGCGG \
        -G CCGTCAATTCMTTTRAGTTT \
        --pair-filter=any \
        --minimum-length 150 \
        --quality-cutoff 20,20 \
        --max-n 0 \
        -j 8 \
        -o trimmed/"${base}_1.fq.gz" \
        -p trimmed/"${base}_2.fq.gz" \
        "${f}" "${r}" \
        > results/cutadapt_logs/"${base}.log" 2>&1
done
# 检查: results/cutadapt_logs/ 下各样本日志, 关注 Pairs written 占比

## 4a. 引物切除结果汇总
conda activate qc_preprocess
cd ${wd}

echo "=== 序列通过率统计 ===" > results/cutadapt_logs/summary_report.txt

for f in seq/*_1.fq.gz; do
    base=$(basename "${f}" _1.fq.gz)
    log="results/cutadapt_logs/${base}.log"

    [ ! -f "${log}" ] && echo "警告: 日志文件缺失 ${log}" | tee -a results/cutadapt_logs/summary_report.txt && continue

    echo "样本: ${base}" >> results/cutadapt_logs/summary_report.txt
    grep "Total read pairs processed" "${log}" >> results/cutadapt_logs/summary_report.txt
    grep "Pairs written" "${log}" >> results/cutadapt_logs/summary_report.txt
    echo "" >> results/cutadapt_logs/summary_report.txt
done
# 检查: results/cutadapt_logs/summary_report.txt 各样本的 Pairs written 占比

## 5. 修剪后质控
conda activate qc_preprocess
cd ${wd}

fastqc trimmed/*_1.fq.gz trimmed/*_2.fq.gz -t 8 -o results/fastqc_trimmed
multiqc results/fastqc_trimmed/ -o results/fastqc_trimmed/ -n multiqc_report_trimmed.html
# 检查: results/fastqc_trimmed/multiqc_report_trimmed.html, 平均质量值需 ≥Q25

## 6. 生成 QIIME2 manifest 文件
conda activate qiime2-2025.7
cd ${wd}

awk -v pwd="$PWD" 'BEGIN {FS=OFS="\t"}
    NR==1 {print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"}
    NR>1 {print $1, pwd "/trimmed/" $1 "_1.fq.gz", pwd "/trimmed/" $1 "_2.fq.gz"}' \
    ${metadata} > manifest
# 检查: head manifest, 确认路径和样本 ID 正确

## 7. 导入 QIIME2
conda activate qiime2-2025.7
cd ${wd}

time qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest \
    --output-path qiime2/demux.qza \
    --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
    --i-data qiime2/demux.qza \
    --o-visualization qiime2/demux.qzv
# 检查: qiime2/demux.qzv, 查看各样本测序量分布

## 8. DADA2 去噪 (参数根据数据质量调整)
conda activate qiime2-2025.7
cd ${wd}

nohup qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qiime2/demux.qza \
    --p-trim-left-f 0 --p-trim-left-r 0 \
    --p-trunc-len-f 230 --p-trunc-len-r 217 \
    --p-max-ee-f 5.0 --p-max-ee-r 5.0 \
    --p-trunc-q 10 \
    --p-n-threads 8 \
    --p-chimera-method consensus \
    --o-table qiime2/table.qza \
    --o-representative-sequences qiime2/rep-seqs.qza \
    --o-denoising-stats qiime2/denoising-stats.qza \
    > logs/dada2.log 2>&1 &
echo "DADA2 已后台运行, 查看进度: tail -f logs/dada2.log"
echo "完成后继续运行后续统计命令"
# 注: 'pooled' 选项在 QIIME2 2025.4 中移除, 2025.7 仅支持 consensus 和 none

qiime metadata tabulate \
    --m-input-file qiime2/denoising-stats.qza \
    --o-visualization qiime2/denoising-stats.qzv

qiime feature-table summarize \
    --i-table qiime2/table.qza \
    --o-visualization qiime2/table.qzv \
    --m-sample-metadata-file ${metadata}

qiime feature-table tabulate-seqs \
    --i-data qiime2/rep-seqs.qza \
    --o-visualization qiime2/rep-seqs.qzv
# 检查: qiime2/table.qzv 样本测序深度, qiime2/denoising-stats.qzv 去噪效率

## 9. 物种注释
conda activate qiime2-2025.7
cd ${wd}

nohup qiime feature-classifier classify-sklearn \
    --i-classifier classifiers/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza \
    --i-reads qiime2/rep-seqs.qza \
    --o-classification qiime2/taxonomy.qza \
    --p-n-jobs 8 \
    > logs/classify.log 2>&1 &
echo "物种注释已后台运行, 查看进度: tail -f logs/classify.log"

qiime metadata tabulate \
    --m-input-file qiime2/taxonomy.qza \
    --o-visualization qiime2/taxonomy.qzv

## 10. 系统发育树构建
conda activate qiime2-2025.7
cd ${wd}

nohup qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences qiime2/rep-seqs.qza \
    --o-alignment qiime2/aligned-rep-seqs.qza \
    --o-masked-alignment qiime2/masked-aligned-rep-seqs.qza \
    --o-tree qiime2/unrooted-tree.qza \
    --o-rooted-tree qiime2/rooted-tree.qza \
    --p-n-threads 8 \
    > logs/tree.log 2>&1 &
echo "系统发育树已后台运行, 查看进度: tail -f logs/tree.log"

## 11. Alpha 稀疏曲线 (选择抽平深度)
conda activate qiime2-2025.7
cd ${wd}

# --p-max-depth 参考 table.qzv 最大样本序列数, 请根据实际修改
qiime diversity alpha-rarefaction \
    --i-table qiime2/table.qza \
    --i-phylogeny qiime2/rooted-tree.qza \
    --p-max-depth 11576 \
    --m-metadata-file ${metadata} \
    --o-visualization results/alpha-rarefaction.qzv
# 检查: results/alpha-rarefaction.qzv, 观察曲线平台期
# 平台期对应深度即为合适的抽平深度, 同时参考 table.qzv 的样本序列数分布
# 【将下方 sampling_depth 替换为选择的值】

sampling_depth=11576

## 12. 核心多样性分析 (Alpha + Beta)
conda activate qiime2-2025.7
cd ${wd}

rm -rf results/core-metrics-results
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny qiime2/rooted-tree.qza \
    --i-table qiime2/table.qza \
    --p-sampling-depth ${sampling_depth} \
    --m-metadata-file ${metadata} \
    --output-dir results/core-metrics-results
# 输出: faith_pd / shannon / observed_features / evenness (alpha)
#       unweighted_unifrac / weighted_unifrac / bray_curtis / jaccard (beta)

## 13. 物种组成柱状图
conda activate qiime2-2025.7
cd ${wd}

qiime taxa barplot \
    --i-table qiime2/table.qza \
    --i-taxonomy qiime2/taxonomy.qza \
    --m-metadata-file ${metadata} \
    --o-visualization results/taxa-bar-plots.qzv

## 14. 导出 QIIME2 结果为纯文本格式 (供 R/Python 等下游分析使用)
conda activate qiime2-2025.7
cd ${wd}

mkdir -p results/export

# 14a. 抽平 ASV 丰度表 (先导出→改名，避免被后续 feature-table.biom 覆盖)
qiime tools export \
    --input-path results/core-metrics-results/rarefied_table.qza \
    --output-path results/export
mv results/export/feature-table.biom results/export/rarefied_table.biom
biom convert -i results/export/rarefied_table.biom \
    -o results/export/rarefied_table.tsv --to-tsv
# 输出: results/export/rarefied_table.tsv (抽平后的 ASV 丰度表)

# 14b. 原始 ASV/OTU 特征表
qiime tools export --input-path qiime2/table.qza \
    --output-path results/export
biom convert -i results/export/feature-table.biom \
    -o results/export/feature-table.tsv --to-tsv
# 查看总序列数（total_frequency）
biom summarize-table -i results/export/feature-table.biom
# 输出: results/export/feature-table.tsv (行为 ASV, 列为样本, 值为序列计数)

# 14c. 物种注释
qiime tools export --input-path qiime2/taxonomy.qza \
    --output-path results/export
# 输出: results/export/taxonomy.tsv (ASV ID → 物种分类 + 置信度)

# 14d. 代表序列
qiime tools export --input-path qiime2/rep-seqs.qza \
    --output-path results/export
# 输出: results/export/dna-sequences.fasta

# 14e. 去噪统计
qiime tools export --input-path qiime2/denoising-stats.qza \
    --output-path results/export

# 14f. Alpha 多样性 (Faith PD / Shannon / Observed Features / Evenness)
qiime tools export \
    --input-path results/core-metrics-results/faith_pd_vector.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/shannon_vector.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/observed_features_vector.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/evenness_vector.qza \
    --output-path results/export

# 14g. Beta 多样性距离矩阵
qiime tools export \
    --input-path results/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/weighted_unifrac_distance_matrix.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/bray_curtis_distance_matrix.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/jaccard_distance_matrix.qza \
    --output-path results/export

# 14h. PCoA 坐标
qiime tools export \
    --input-path results/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/weighted_unifrac_pcoa_results.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/bray_curtis_pcoa_results.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/jaccard_pcoa_results.qza \
    --output-path results/export

# 14i. 系统发育树 (Newick 格式)
qiime tools export \
    --input-path qiime2/rooted-tree.qza \
    --output-path results/export

## 15. 一站式 R 可视化 (需先将 R 脚本放到项目根目录)
conda activate qiime2-2025.7
cd ${wd}

# 将 examples/QIIME2_16S_visualization.R 复制到 $(pwd)
# 修改 setwd() 后即可运行, 脚本将自动读取 results/export/ 中的数据
# 生成图表至 results/export/ 下各子目录, 并准备 FAPROTAX/PICRUSt2 输入文件
if [ -f "QIIME2_16S_visualization.R" ]; then
    Rscript QIIME2_16S_visualization.R
else
    echo "请先将 examples/QIIME2_16S_visualization.R 复制到 $(pwd)"
fi

## 16. FAPROTAX 功能预测分析
conda activate qiime2-2025.7
cd ${wd}

cd ${wd}/results/export/faprotax

# FAPROTAX 脚本路径 (用户按实际修改)
# 下载: http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Download
sd=~/db/EasyMicrobiome/script/FAPROTAX_1.2.12

# 制作输入 OTU 表 (添加物种注释到 BIOM)
biom add-metadata \
    -i rarefied_table.biom \
    --observation-metadata-fp taxonomy.tsv \
    -o rarefied_tax.biom \
    --sc-separated taxonomy \
    --observation-header OTUID,taxonomy

# FAPROTAX 功能预测
python ${sd}/collapse_table.py \
    -i rarefied_tax.biom \
    -g ${sd}/FAPROTAX.txt \
    --collapse_by_metadata 'taxonomy' \
    -v --force \
    -o faprotax.txt \
    -r faprotax_report.txt

# 制作功能注释有无矩阵
grep '*' -B 1 faprotax_report.txt | grep -v -P '^--$' > faprotax_report.clean
perl ${sd}/../faprotax_report_sum.pl \
    -i faprotax_report.clean \
    -o faprotax_report

echo "FAPROTAX 完成: results/export/faprotax/faprotax_report.txt"
echo "  功能有无矩阵: results/export/faprotax/faprotax_report.mat"
echo "  功能丰度表:   results/export/faprotax/faprotax.txt"

## 17. PICRUSt2 功能预测分析
conda activate picrust2
cd ${wd}

cd ${wd}/results/export/picrust2

# 运行 PICRUSt2 管道 (后台运行，等待完成后继续)
nohup picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.tsv \
    -o ./out -p 16 > picrust2.log 2>&1 &
echo "PICRUSt2 已后台运行, 查看进度: tail -f results/export/picrust2/picrust2.log"

# 添加 EC/KO/Pathway 注释 (等待 PICRUSt2 管道完成后执行)
cd out
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o METACYC.tsv
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  -o EC.tsv
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  -o KO.tsv

# KEGG 按层级合并 (需 EasyMicrobiome 数据库)
db=~/db/EasyMicrobiome/
python3 ${db}/script/summarizeAbundance.py \
    -i KO.tsv \
    -m ${db}/kegg/KO1-4.txt \
    -c 2,3,4 -s ',+,+,' -n raw \
    -o KEGG
wc -l KEGG*

# 导出与功能预测完成: results/export/
# 关键文件:
#   rarefied_table.tsv     — 抽平 ASV 丰度表
#   feature-table.tsv      — ASV/OTU 丰度表
#   taxonomy.tsv           — 物种注释
#   dna-sequences.fasta    — 代表序列
#   *.tsv (alpha)          — Alpha 多样性 (Faith PD, Shannon 等)
#   *_distance_matrix.tsv  — Beta 多样性距离矩阵
#   *_pcoa_results.tsv     — PCoA 降维坐标
#   faprotax/              — FAPROTAX 功能预测结果
#   picrust2/              — PICRUSt2 功能预测结果
#
# 一站式可视化: 将 QIIME2_16S_visualization.R 放在项目根目录 (与 metadata.txt 同级),
# 修改 setwd() 后运行 Rscript QIIME2_16S_visualization.R 即可自动读取 results/export/
# 并生成图表至 results/export/ 下各子目录, 无需手动搬运文件

# 输出项目路径与文件结构
cd ${wd}
find . | sed -e 's/[^/]*\//|   /g' -e 's/| *[^|]*$/|-- &/' > Project_file_structure.log
