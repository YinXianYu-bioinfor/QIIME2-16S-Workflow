# QIIME2-16S-Workflow

A QIIME2-based 16S rRNA amplicon analysis pipeline for paired-end sequencing data.

## Overview

This pipeline covers the complete 16S amplicon data analysis workflow:

1. **QC & Primer Trimming** — FastQC quality inspection + cutadapt primer removal
2. **Denoising** — DADA2 denoising (error correction, ASV inference, chimera removal)
3. **Taxonomy Classification** — SILVA-138 classifier-based species annotation
4. **Phylogenetic Tree** — MAFFT alignment + FastTree phylogeny
5. **Diversity Analysis** — Alpha rarefaction, core metrics (Faith PD, Shannon, UniFrac), PCoA
6. **Export** — Export QIIME2 artifacts to text/TSV/FASTA for downstream analysis

## Prerequisites

- Miniconda (or Anaconda)
- Two conda environments (see install script):
  - `qc_preprocess`: FastQC, cutadapt, MultiQC
  - `qiime2-2025.7`: QIIME2 amplicon distribution (2023.2 available as fallback)
- SILVA-138.2 classifier (`SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza`)
  - **Compatibility**: SILVA138.2 requires qiime2 ≥ 2023.x (this pipeline uses 2025.7)
  - For qiime2-2023.2, use `silva-138-99-nb-classifier.qza` instead

## Quick Start

### 1. Environment Setup

Open the install script and **manually run the commands step by step**:
- `qiime2-16s-pipeline_install_en.sh` (English)
- `qiime2-16s-pipeline_install.sh` (中文版)

Some steps (e.g., `conda activate`) require interactive shell execution and cannot be run as a single `bash script.sh`.

### 2. Data Preparation

Place raw paired-end FASTQ files in `seq/` directory with naming format:
```
seq/<sample>_1.fq.gz
seq/<sample>_2.fq.gz
```

Place metadata file (`metadata.txt`) in the working directory (TSV format, first column = sample ID, must include Group column).

### 3. Run Pipeline

Edit working directory and parameters in `qiime2-16s-pipeline.sh` (or `qiime2-16s-pipeline_en.sh` for English), then execute step by step:

```bash
# Step 0: Edit parameters in the script (wd, metadata path, etc.)

# Step 1-2: Directory setup & environment validation
# Step 3-5: QC & primer trimming
# Step 6-7: QIIME2 import
# Step 8: DADA2 denoising
# Step 9: Taxonomy classification
# Step 10: Phylogenetic tree
# Step 11-13: Diversity analysis
# Step 14: Export results
```

**Note**: This pipeline is designed for **manual step-by-step execution**. Check output quality at each stage before proceeding to the next.

## Project Structure

### Runtime directory structure (after pipeline completes)

```
project_root/                         # your working directory (wd)
├── metadata.txt                      # sample metadata (TSV, required)
├── manifest                          # manifest file (auto-generated)
├── QIIME2_16S_visualization.R        # R visualization script (copy here)
├── seq/                              # raw paired-end FASTQ files
│   ├── sample1_1.fq.gz
│   ├── sample1_2.fq.gz
│   └── ...
├── trimmed/                          # primer-trimmed reads (cutadapt output)
│   ├── sample1_1.fq.gz
│   ├── sample1_2.fq.gz
│   └── ...
├── classifiers/                      # taxonomy classifier databases
│   ├── SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza
│   ├── silva-138-99-nb-classifier.qza
│   └── gg_2022_10_backbone_full_length.nb.qza
├── qiime2/                           # QIIME2 artifacts (.qza / .qzv)
│   ├── demux.qza / demux.qzv         # demultiplexing summary
│   ├── table.qza / table.qzv         # feature table
│   ├── denoising-stats.qza / .qzv    # DADA2 denoising statistics
│   ├── rep-seqs.qza / rep-seqs.qzv   # representative sequences
│   ├── aligned-rep-seqs.qza          # MAFFT alignment
│   ├── masked-aligned-rep-seqs.qza   # masked alignment
│   ├── unrooted-tree.qza             # FastTree unrooted tree
│   ├── rooted-tree.qza               # rooted tree (midpoint)
│   └── taxonomy.qza / taxonomy.qzv   # taxonomy classification
├── results/
│   ├── fastqc_raw/                   # FastQC reports (raw reads)
│   │   ├── *_fastqc.html / *.zip
│   │   ├── multiqc_report_raw.html
│   │   └── multiqc_report_raw_data/
│   ├── cutadapt_logs/                # per-sample primer trimming logs
│   │   ├── sample1.log
│   │   ├── summary_report.txt
│   │   └── ...
│   ├── fastqc_trimmed/               # FastQC reports (after trimming)
│   │   ├── *_fastqc.html / *.zip
│   │   ├── multiqc_report_trimmed.html
│   │   └── multiqc_report_trimmed_data/
│   ├── alpha-rarefaction.qzv         # alpha rarefaction curve
│   ├── taxa-bar-plots.qzv            # taxonomy bar plot
│   ├── core-metrics-results/         # core diversity metrics
│   │   ├── rarefied_table.qza
│   │   ├── shannon_vector.qza
│   │   ├── observed_features_vector.qza
│   │   ├── faith_pd_vector.qza
│   │   ├── evenness_vector.qza
│   │   ├── bray_curtis_distance_matrix.qza
│   │   ├── jaccard_distance_matrix.qza
│   │   ├── unweighted_unifrac_distance_matrix.qza
│   │   ├── weighted_unifrac_distance_matrix.qza
│   │   ├── bray_curtis_pcoa_results.qza
│   │   ├── jaccard_pcoa_results.qza
│   │   ├── unweighted_unifrac_pcoa_results.qza
│   │   ├── weighted_unifrac_pcoa_results.qza
│   │   ├── bray_curtis_emperor.qzv
│   │   ├── jaccard_emperor.qzv
│   │   ├── unweighted_unifrac_emperor.qzv
│   │   └── weighted_unifrac_emperor.qzv
│   └── export/                       # exported data (read by R script)
│       ├── feature-table.tsv         # raw ASV abundance table
│       ├── feature-table.biom        # raw ASV abundance table (BIOM)
│       ├── taxonomy.tsv              # species annotation
│       ├── dna-sequences.fasta       # representative sequences
│       ├── rarefied_table.tsv        # rarefied ASV abundance table
│       ├── rarefied_table.biom       # rarefied ASV abundance table (BIOM)
│       ├── alpha-diversity.tsv       # alpha diversity metrics
│       ├── distance-matrix.tsv       # beta diversity distance matrix
│       ├── ordination.txt            # PCoA ordination coordinates
│       ├── stats.tsv                 # DADA2 denoising stats
│       ├── tree.nwk                  # phylogenetic tree (Newick)
│       ├── alpha/                    # alpha diversity boxplots *
│       │   ├── alpha_diversity_boxplot.pdf
│       │   └── rarefaction_curves.pdf
│       ├── beta/                     # beta diversity PCoA plots *
│       │   ├── beta_diversity_pcoa_bray_curtis.pdf
│       │   └── beta_diversity_pcoa_jaccard.pdf
│       ├── taxa/                     # phylum composition plots *
│       │   ├── phylum_stacked_barplot.pdf
│       │   └── phylum_abundance_barchart.pdf
│       ├── heatmap/                  # genus-level heatmap *
│       │   └── genus_heatmap.pdf
│       ├── faprotax/                 # FAPROTAX functional prediction
│       │   ├── faprotax.txt
│       │   ├── faprotax_report.txt
│       │   ├── faprotax_report.clean
│       │   ├── faprotax_report.mat
│       │   ├── faprotax_report.func_otu
│       │   ├── faprotax_report.otu_func
│       │   ├── taxonomy.tsv
│       │   ├── rarefied_table.biom
│       │   └── rarefied_tax.biom
│       ├── picrust2/                 # PICRUSt2 functional prediction
│       │   ├── feature-table.tsv
│       │   ├── dna-sequences.fasta
│       │   └── out/
│       │       ├── EC.tsv / KO.tsv
│       │       ├── EC_predicted.tsv.gz / KO_predicted.tsv.gz
│       │       ├── EC_metagenome_out/
│       │       │   ├── pred_metagenome_unstrat.tsv.gz
│       │       │   ├── seqtab_norm.tsv.gz
│       │       │   └── weighted_nsti.tsv.gz
│       │       ├── KO_metagenome_out/
│       │       │   ├── pred_metagenome_unstrat.tsv.gz
│       │       │   ├── seqtab_norm.tsv.gz
│       │       │   └── weighted_nsti.tsv.gz
│       │       ├── pathways_out/
│       │       │   └── path_abun_unstrat.tsv.gz
│       │       ├── METACYC.tsv
│       │       ├── marker_predicted_and_nsti.tsv.gz
│       │       ├── KEGG.Pathway.raw.txt
│       │       ├── KEGG.PathwayL1.raw.txt
│       │       ├── KEGG.PathwayL2.raw.txt
│       │       └── out.tre
│       └── feature_tables/           # processed tables *
│           ├── taxonomy_processed.tsv
│           ├── feature_table_with_taxonomy.tsv
│           ├── genus_abundance.tsv
│           └── alpha_diversity_metrics.tsv
└── logs/                             # pipeline run logs
    ├── dada2.log
    ├── classify.log
    └── tree.log
```

> **Note:** Directories marked with `*` (alpha/, beta/, taxa/, heatmap/, feature_tables/) are generated by the R visualization script, not by the shell pipeline directly.

### Repository structure (this repository)

```
QIIME2-16S-Workflow/
├── qiime2-16s-pipeline.sh           # Main analysis pipeline (中文版)
├── qiime2-16s-pipeline_en.sh        # Main analysis pipeline (English)
├── qiime2-16s-pipeline_install.sh   # Environment installation (中文版)
├── qiime2-16s-pipeline_install_en.sh# Environment installation (English)
├── examples/                         # Example results for preview only
│   ├── metadata.txt                  # Sample metadata reference
│   ├── manifest                      # Sample manifest reference
│   ├── QIIME2_16S_visualization.R    # R visualization script (中文版)
│   ├── QIIME2_16S_visualization_EN.R # R visualization script (English)
│   ├── export/                       # Preview of exported results
│   ├── qiime2/                       # Preview of QIIME2 visualizations
│   └── logs/                         # Example pipeline run logs
├── LICENSE
├── .gitignore
└── README.md
```

## Example Data

Dataset: *Arabidopsis thaliana* rhizosphere microbiome (16S rRNA, V3-V4 region)
- 18 samples (WT/KO/OE groups, 6 replicates each)
- Platform: Illumina HiSeq 2500, PE250
- Published under CRA002352

View `.qzv` files online at: https://view.qiime2.org/

> **⚠️ Note on example data:** The `examples/export/` directory in this repository contains a **lightweight preview subset** of pipeline output — it is intended **only** to demonstrate visualization quality and script behavior. **Do not use these files for actual analysis or as reference data.** Real pipeline runs generate the complete output in `results/export/` (see [Project Structure](#project-structure) and [R Visualization](#r-visualization) sections). The shell pipeline script (`qiime2-16s-pipeline.sh`) produces the full set of output files automatically.

## R Visualization

After the pipeline completes Step 14 (Export), you can use the R scripts for publication-ready visualizations.

> **Path note:** The R scripts in `examples/` use `export/` as the data path — this matches the simplified demo layout in this repo (`examples/export/`). In real pipeline runs, QIIME2 exports data to **`results/export/`** (see directory structure below). When using the R script with actual pipeline output, update the config section paths from `export/` to `results/export/`.

### Real pipeline directory structure

After running the pipeline, your project directory structure is detailed in the [Project Structure](#project-structure) section above. The key layout for R visualization is:

```
project_root/                         # your working directory (wd)
├── metadata.txt
├── QIIME2_16S_visualization_EN.R     # copy R script here
├── results/
│   └── export/                       # pipeline output, read by R script
│       ├── feature-table.tsv         # input for R visualization
│       ├── taxonomy.tsv
│       ├── dna-sequences.fasta
│       ├── rarefied_table.tsv / .biom
│       ├── alpha/                    # R script outputs here
│       ├── beta/
│       ├── taxa/
│       ├── heatmap/
│       ├── faprotax/
│       ├── picrust2/
│       └── feature_tables/
└── logs/
```

> **⚠️ Important:** The `examples/export/` directory in this repository is a **lightweight preview** of selected pipeline outputs. It is intended only to demonstrate visualization quality. **Do not use `examples/export/` files for actual analysis or as reference data.** Always use the full `results/export/` generated by your own pipeline run.

### Workflow

> **Note:** The R script can be run **directly on the server** as part of the pipeline (section 15 of the shell script). No need to download files to your local machine.

```bash
# 1. Copy the R script to your project root (same level as metadata.txt):
#    cp examples/QIIME2_16S_visualization_EN.R /path/to/project/

# 2. Open the copied script and modify:
#    - setwd() → set to your project root
#    - Config paths: change "export/" to "results/export/"
#      e.g., feature_table_file <- "results/export/feature-table.tsv"

# 3. Run directly on the server:
#    Rscript QIIME2_16S_visualization_EN.R
```

### Output preview

The R scripts generate 6 types of plots and 4 processed data tables under `export/` (or `results/export/` in real runs):
- **alpha/** — Alpha diversity boxplots (Shannon, Observed features, Faith PD, Evenness)
- **beta/** — PCoA ordination plots (Bray-Curtis, Jaccard)
- **taxa/** — Phylum-level stacked barplot and abundance barchart
- **heatmap/** — Genus-level abundance heatmap (Top N genera)
- **faprotax/** — FAPROTAX functional prediction input tables
- **picrust2/** — PICRUSt2 input files (BIOM + FASTA)
- **feature_tables/** — Processed taxonomy table and genus abundance table

You can preview the generated PDFs in this repository under `examples/export/alpha/`, `examples/export/beta/`, etc. to evaluate visualization quality before running on your own data. However, **remember that `examples/export/` is a preview subset — your actual pipeline output will be in `results/export/` and will contain many more files.**

## Notes

- DADA2 parameters (`--p-trunc-len-f`, `--p-trunc-len-r`, `--p-max-ee`) should be adjusted based on your sequencing data quality
- `--p-sampling-depth` for core diversity metrics should be determined from alpha rarefaction curve and table.qzv
- DADA2 chimera method: qiime2-2025.7 supports `consensus` and `none` only (`pooled` was removed in QIIME2 2025.4)

## License

[MIT](LICENSE)

---

## 中文说明

### 项目简介

基于 QIIME2 的 16S rRNA 扩增子分析流程，适用于双端测序数据。

### 分析流程

1. **质控与引物切除** — FastQC 质量检查 + cutadapt 引物去除
2. **去噪** — DADA2 去噪（纠错、ASV 推断、嵌合体去除）
3. **物种注释** — 基于 SILVA-138 分类器的物种注释
4. **系统发育树** — MAFFT 比对 + FastTree 建树
5. **多样性分析** — Alpha 稀疏曲线、核心多样性指标（Faith PD、Shannon、UniFrac）、PCoA
6. **结果导出** — 将 QIIME2 产物导出为文本/TSV/FASTA 格式，供 R/Python 下游分析

### 环境要求

- Miniconda（或 Anaconda）
- 两个 conda 环境（详见安装脚本）：
  - `qc_preprocess`：FastQC、cutadapt、MultiQC
  - `qiime2-2025.7`：QIIME2 amplicon 发行版（备选 2023.2）
- SILVA-138.2 分类器（`SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza`）
  - **兼容性**：SILVA138.2 要求 qiime2 ≥ 2023 版本（本流程使用 2025.7）
  - 若使用 qiime2-2023.2，请改用 `silva-138-99-nb-classifier.qza`

### 快速开始

#### 1. 环境安装

打开安装脚本，**逐条手动执行**其中的命令：
- `qiime2-16s-pipeline_install.sh`（中文版）
- `qiime2-16s-pipeline_install_en.sh`（English）

部分步骤（如 `conda activate`）需在交互式 shell 中执行，无法通过 `bash script.sh` 一键运行。

#### 2. 数据准备

将双端 FASTQ 文件放入 `seq/` 目录，命名格式：
```
seq/<样本名>_1.fq.gz
seq/<样本名>_2.fq.gz
```

将元数据文件 `metadata.txt` 放在工作目录下（TSV 格式，第一列为样本 ID，必须包含 Group 分组列）。

#### 3. 运行流程

编辑 `qiime2-16s-pipeline.sh` 中的工作目录和参数，然后按步骤逐条执行：

```bash
# 第0步：修改脚本中的参数（wd、metadata 路径等）
# 第1-2步：目录初始化与环境验证
# 第3-5步：质控与引物切除
# 第6-7步：导入 QIIME2
# 第8步：DADA2 去噪
# 第9步：物种注释
# 第10步：系统发育树
# 第11-13步：多样性分析
# 第14步：结果导出
```

**注意**：本流程设计为**手动分步执行**，每完成一步请检查输出质量，确认无误后再继续下一步。

### 示例数据

数据集：拟南芥根际微生物组（16S rRNA，V3-V4 区）
- 18 个样本（WT/KO/OE 三组，每组 6 个重复）
- 测序平台：Illumina HiSeq 2500，PE250
- 数据编号：CRA002352

`.qzv` 可视化文件可在 https://view.qiime2.org/ 在线查看。

> **⚠️ 关于示例数据的说明：** 本仓库 `examples/export/` 目录仅包含 pipeline 产出的**轻量预览子集**，**仅供预览图表效果和脚本表现，不可用于实际分析或作为参考数据**。实际运行 pipeline 时，完整输出位于 `results/export/`，包含所有分析文件（参见[项目结构](#项目结构)和 [R 可视化](#r-可视化)章节）。shell 管道脚本（`qiime2-16s-pipeline.sh`）会自动生成全部输出文件。

### R 可视化

流程第 14 步（导出）完成后，可使用 R 脚本进行出版级可视化。

> **路径说明：** `examples/` 中的 R 脚本使用 `export/` 路径，与本仓库简化的示例目录结构（`examples/export/`）一致。实际运行 pipeline 时，QIIME2 导出数据位于 **`results/export/`**。将 R 脚本用于真实数据时，需将配置区的路径从 `export/` 改为 `results/export/`（参见下方目录结构）。

### 真实运行时的目录结构

pipeline 运行完成后，完整的目录结构请参见上方的[项目结构](#项目结构)章节。R 可视化涉及的关键路径如下：

```
项目根目录/                          # 你的工作目录 (wd)
├── metadata.txt
├── QIIME2_16S_visualization.R       # R 脚本放在这里
├── results/
│   └── export/                      # pipeline 导出数据，R 脚本读取此处
│       ├── feature-table.tsv         # R 可视化输入文件
│       ├── taxonomy.tsv
│       ├── dna-sequences.fasta
│       ├── rarefied_table.tsv / .biom
│       ├── alpha/                   # R 脚本输出图表至此
│       ├── beta/
│       ├── taxa/
│       ├── heatmap/
│       ├── faprotax/
│       ├── picrust2/
│       └── feature_tables/
└── logs/
```

> **⚠️ 重要提示：** 本仓库中的 `examples/export/` 目录仅为**精简预览**，展示部分 pipeline 输出效果，**不可用于实际分析或作为参考数据**。请始终使用你自己运行 pipeline 生成的完整 `results/export/` 目录。

#### 操作流程

> **注意：** R 脚本可以直接在服务器上作为 pipeline 的一部分运行（shell 脚本第 15 节），无需下载到本地。

```bash
# 1. 将 R 脚本复制到项目根目录（与 metadata.txt 同级）:
#    cp examples/QIIME2_16S_visualization.R /path/to/project/

# 2. 打开复制的脚本，修改以下内容:
#    - setwd() → 设置为项目根目录
#    - 配置区路径：将 "export/" 改为 "results/export/"
#      例如: feature_table_file <- "results/export/feature-table.tsv"

# 3. 直接在服务器上运行:
#    Rscript QIIME2_16S_visualization.R
```

#### 产出预览

R 脚本在 `export/`（真实运行时为 `results/export/`）下生成 6 类图表和 4 张处理后的数据表：
- **alpha/** — Alpha 多样性箱线图（Shannon、Observed features、Faith PD、Evenness）
- **beta/** — PCoA 降维图（Bray-Curtis、Jaccard）
- **taxa/** — 门水平堆叠柱状图和丰度条形图
- **heatmap/** — 属水平丰度热图（Top N 属）
- **faprotax/** — FAPROTAX 功能预测输入表
- **picrust2/** — PICRUSt2 输入文件（BIOM + FASTA）
- **feature_tables/** — 处理后的分类学表和属水平丰度表

本仓库 `examples/export/` 下包含示例输出 PDF，可预览可视化效果。**但请注意：`examples/export/` 仅为预览子集，你实际运行 pipeline 后，完整输出在 `results/export/` 中。**

### 注意事项

- DADA2 参数（`--p-trunc-len-f`、`--p-trunc-len-r`、`--p-max-ee`）需根据实际测序数据质量调整
- `--p-sampling-depth`（抽平深度）需根据 alpha 稀疏曲线和 table.qzv 确定
- DADA2 嵌合体方法：qiime2-2025.7 仅支持 `consensus` 和 `none`（`pooled` 已在 QIIME2 2025.4 移除）
