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

```
QIIME2-16S-Workflow/
├── qiime2-16s-pipeline.sh          # Main analysis pipeline (中文版)
├── qiime2-16s-pipeline_en.sh       # Main analysis pipeline (English)
├── qiime2-16s-pipeline_install.sh  # Environment installation (中文版)
├── qiime2-16s-pipeline_install_en.sh # Environment installation (English)
├── examples/                        # Example analysis results
│   ├── metadata.txt                 # Sample metadata reference
│   ├── manifest                     # Sample manifest reference
│   ├── export/                      # Final exported results + R visualization
│   │   ├── QIIME2_16S_visualization.R     # R visualization script (中文版)
│   │   ├── QIIME2_16S_visualization_EN.R  # R visualization script (English)
│   │   ├── data/                           # QIIME2 exported data files
│   │   │   ├── feature-table.tsv           # Raw ASV abundance table
│   │   │   ├── feature-table.biom          # Raw ASV abundance table (BIOM)
│   │   │   ├── rarefied_table.tsv          # Rarefied ASV abundance table
│   │   │   ├── rarefied_table.biom         # Rarefied ASV abundance table (BIOM)
│   │   │   ├── taxonomy.tsv                # Species annotation
│   │   │   ├── dna-sequences.fasta         # Representative sequences
│   │   │   ├── alpha-diversity.tsv         # Alpha diversity metrics
│   │   │   ├── distance-matrix.tsv         # Beta diversity distances
│   │   │   ├── ordination.txt              # PCoA coordinates
│   │   │   ├── stats.tsv                   # Denoising statistics
│   │   │   └── tree.nwk                    # Phylogenetic tree (Newick)
│   │   ├── alpha/                   # Alpha diversity boxplots
│   │   ├── beta/                    # Beta diversity PCoA plots
│   │   ├── taxa/                    # Phylum composition plots
│   │   ├── heatmap/                 # Genus-level heatmap
│   │   ├── faprotax/                # FAPROTAX input tables
│   │   ├── picrust2/                # PICRUSt2 input files
│   │   └── feature_tables/          # Processed taxonomy & abundance tables
│   ├── qiime2/                      # QIIME2 visualizations (.qzv)
│   │   ├── demux.qzv               # Demultiplexing summary
│   │   ├── denoising-stats.qzv     # DADA2 denoising statistics
│   │   ├── table.qzv               # Feature table summary
│   │   └── taxonomy.qzv            # Taxonomy classification
│   └── logs/                        # Pipeline run logs
│       ├── dada2.log               # DADA2 denoising log
│       ├── classify.log            # Taxonomy classification log
│       └── tree.log                # Phylogenetic tree log
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

## R Visualization

After the pipeline completes Step 14 (Export), you can use the R scripts for publication-ready visualizations.

### Workflow

```bash
# 1. Download exported files from the server
#    (results/export/ from the pipeline working directory)
#    Save them to a local directory (e.g., data/)

# 2. Place the R script at the same level as data/
#    export/
#    ├── QIIME2_16S_visualization.R
#    ├── data/
#    │   ├── feature-table.tsv
#    │   ├── taxonomy.tsv
#    │   └── ...

# 3. Edit the script configuration:
#    - setwd() → your local export/ path
#    - metadata_file → path to metadata.txt

# 4. Run the script:
#    Rscript QIIME2_16S_visualization.R
```

### Output preview

The R scripts generate 6 types of plots and 4 processed data tables under `export/`:
- **alpha/** — Alpha diversity boxplots (Shannon, Observed features, Faith PD, Evenness)
- **beta/** — PCoA ordination plots (Bray-Curtis, Jaccard)
- **taxa/** — Phylum-level stacked barplot and abundance barchart
- **heatmap/** — Genus-level abundance heatmap (Top N genera)
- **faprotax/** — FAPROTAX functional prediction input tables
- **picrust2/** — PICRUSt2 input files (BIOM + FASTA)
- **feature_tables/** — Processed taxonomy table and genus abundance table

You can preview the generated PDFs in this repository under `examples/export/` to evaluate visualization quality before running on your own data.

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

### R 可视化

流程第 14 步（导出）完成后，可使用 R 脚本进行出版级可视化。

#### 操作流程

```bash
# 1. 从服务器下载导出文件（工作目录下的 results/export/）
#    保存到本地目录（如 data/）

# 2. 将 R 脚本放在与 data/ 同级的目录下
#    export/
#    ├── QIIME2_16S_visualization.R
#    ├── data/
#    │   ├── feature-table.tsv
#    │   ├── taxonomy.tsv
#    │   └── ...

# 3. 修改脚本配置：
#    - setwd() → 本地 export/ 的实际路径
#    - metadata_file → metadata.txt 路径

# 4. 运行脚本：
#    Rscript QIIME2_16S_visualization.R
```

#### 产出预览

R 脚本在 `export/` 下生成 6 类图表和 4 张处理后的数据表：
- **alpha/** — Alpha 多样性箱线图（Shannon、Observed features、Faith PD、Evenness）
- **beta/** — PCoA 降维图（Bray-Curtis、Jaccard）
- **taxa/** — 门水平堆叠柱状图和丰度条形图
- **heatmap/** — 属水平丰度热图（Top N 属）
- **faprotax/** — FAPROTAX 功能预测输入表
- **picrust2/** — PICRUSt2 输入文件（BIOM + FASTA）
- **feature_tables/** — 处理后的分类学表和属水平丰度表

本仓库 `examples/export/` 下包含示例输出 PDF，可预览可视化效果。

### 注意事项

- DADA2 参数（`--p-trunc-len-f`、`--p-trunc-len-r`、`--p-max-ee`）需根据实际测序数据质量调整
- `--p-sampling-depth`（抽平深度）需根据 alpha 稀疏曲线和 table.qzv 确定
- DADA2 嵌合体方法：qiime2-2025.7 仅支持 `consensus` 和 `none`（`pooled` 已在 QIIME2 2025.4 移除）
