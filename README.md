# QIIME2-16S-Workflow

A QIIME2-based 16S rRNA amplicon analysis pipeline for paired-end sequencing data, with R visualization scripts for publication-ready figures.

## Overview

This pipeline covers the complete 16S amplicon data analysis workflow:

1. **QC & Primer Trimming** — FastQC quality inspection + cutadapt primer removal
2. **Denoising** — DADA2 (error correction, ASV inference, chimera removal)
3. **Taxonomy Classification** — SILVA-138 classifier-based species annotation
4. **Phylogenetic Tree** — MAFFT alignment + FastTree
5. **Diversity Analysis** — Alpha rarefaction, core metrics (Faith PD, Shannon, UniFrac), PCoA
6. **Export** — Export QIIME2 artifacts to TSV/FASTA/Biom for downstream analysis
7. **Visualization** — R scripts for publication-ready plots (alpha/beta diversity, taxonomy, heatmap, PICRUSt2, FAPROTAX)

## Repository Contents

```
QIIME2-16S-Workflow/
├── scripts/                          # All analysis scripts (see table below)
├── metadata.txt                      # Example sample metadata (TSV)
├── qiime2/                           # Example .qzv visualizations (view on view.qiime2.org)
├── results/export/                   # Example pipeline outputs (preview only)
│   ├── alpha/                        #   Alpha diversity boxplots
│   ├── beta/                         #   Beta diversity PCoA plots
│   ├── taxa/                         #   Phylum composition plots
│   ├── heatmap/                      #   Genus-level heatmap
│   ├── feature_tables/               #   Processed abundance tables
│   ├── faprotax/                     #   FAPROTAX functional prediction
│   └── picrust2/                     #   PICRUSt2 functional prediction
├── Project_file_structure.log        # Auto-generated directory tree
├── LICENSE
├── .gitignore
└── README.md
```

> **Note:** The `results/export/` directory in this repo serves as a **preview** of pipeline outputs. Full runtime outputs include additional files (see [Runtime Output Structure](#runtime-output-structure)). Real pipeline runs generate fresh results in `results/export/`.

## Scripts

All scripts are stored in `scripts/`. Each has a Chinese (`.R` / `.sh`) and English (`_EN.R` / `_en.sh`) version.

| Script | Purpose |
|--------|---------|
| `qiime2-16s-pipeline.sh` | Main analysis pipeline — QC, DADA2, taxonomy, tree, diversity, export |
| `qiime2-16s-pipeline_install.sh` | Environment setup — create conda environments for QIIME2 and QC tools |
| `QIIME2_16S_visualization.R` | Publication-ready R plots — alpha/beta diversity, phylum, genus heatmap |
| `picrust2_visualization.R` | PICRUSt2 functional prediction visualization — NSTI, KEGG, KO, EC, MetaCyc |
| `faprotax_visualization.R` | FAPROTAX ecological function visualization — cycles, PCoA, OTU contribution |

> **Usage:** Copy the needed script(s) to your project working directory (same level as `metadata.txt`). The scripts default to `results/export/` paths — no manual path adjustment needed.

## Quick Start

### 1. Environment Setup

Open `scripts/qiime2-16s-pipeline_install.sh` (or `_en.sh`) and **manually run commands step by step**. Some steps (`conda activate`) require interactive shell execution.

### 2. Data Preparation

Place raw paired-end FASTQ files in `seq/`:
```
seq/<sample>_1.fq.gz
seq/<sample>_2.fq.gz
```

Place `metadata.txt` in the working directory (TSV, first column = sample ID, must include `Group` column).

### 3. Run Pipeline

Edit `scripts/qiime2-16s-pipeline.sh` (or `_en.sh`) — set `wd`, metadata path, and DADA2 parameters — then execute **step by step**:

```bash
# Steps 1-2:  Directory setup & environment validation
# Steps 3-5:  QC & primer trimming (FastQC + cutadapt + MultiQC)
# Steps 6-7:  QIIME2 import & demultiplexing
# Step 8:     DADA2 denoising
# Step 9:     Taxonomy classification (SILVA)
# Step 10:    Phylogenetic tree (MAFFT + FastTree)
# Steps 11-13: Diversity analysis (rarefaction, core metrics)
# Step 14:    Export results to TSV/FASTA/Biom
# Step 15:    R visualization
# Step 16:    PICRUSt2 & FAPROTAX functional prediction + visualization
```

> The pipeline is designed for **manual step-by-step execution**. Check output quality at each stage before proceeding.

## Prerequisites

- **Miniconda** (or Anaconda)
- **Two conda environments:**
  - `qc_preprocess` — FastQC, cutadapt, MultiQC
  - `qiime2-2025.7` — QIIME2 amplicon distribution (2023.2 available as fallback)
- **SILVA-138.2 classifier** — `SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza` (compatible with qiime2 ≥ 2023.x)
  - For qiime2-2023.2, use `silva-138-99-nb-classifier.qza` instead

## Runtime Output Structure

After the pipeline completes, your working directory contains:

```
project_root/                        # Your working directory (wd)
├── metadata.txt                     # Sample metadata
├── manifest                         # Auto-generated manifest file
├── seq/                             # Raw FASTQ files (user-provided)
│   ├── sample1_1.fq.gz
│   ├── sample1_2.fq.gz
│   └── ...
├── trimmed/                         # Primer-trimmed reads (cutadapt)
│   ├── sample1_1.fq.gz
│   ├── sample1_2.fq.gz
│   └── ...
├── classifiers/                     # Taxonomy classifier databases
│   ├── SILVA138.2_*_classifier_full-length.qza
│   ├── silva-138-99-nb-classifier.qza
│   └── gg_2022_10_backbone_full_length.nb.qza
├── qiime2/                          # QIIME2 artifacts (.qza / .qzv)
│   ├── demux.qza / demux.qzv          # Demultiplexing summary
│   ├── table.qza / table.qzv          # Feature table
│   ├── denoising-stats.qza / .qzv     # DADA2 statistics
│   ├── rep-seqs.qza / rep-seqs.qzv    # Representative sequences
│   ├── aligned-rep-seqs.qza           # MAFFT alignment
│   ├── masked-aligned-rep-seqs.qza    # Masked alignment
│   ├── unrooted-tree.qza              # FastTree phylogeny
│   ├── rooted-tree.qza                # Rooted tree (midpoint)
│   └── taxonomy.qza / taxonomy.qzv    # Taxonomy classification
├── results/
│   ├── fastqc_raw/                   # FastQC reports (raw reads)
│   │   ├── *_fastqc.html / *.zip
│   │   ├── multiqc_report_raw.html
│   │   └── multiqc_report_raw_data/
│   ├── fastqc_trimmed/               # FastQC reports (trimmed)
│   │   ├── *_fastqc.html / *.zip
│   │   ├── multiqc_report_trimmed.html
│   │   └── multiqc_report_trimmed_data/
│   ├── cutadapt_logs/                # Per-sample primer trimming logs
│   │   ├── sample1.log / sample2.log / ...
│   │   └── summary_report.txt
│   ├── alpha-rarefaction.qzv         # Alpha rarefaction curve
│   ├── taxa-bar-plots.qzv            # Taxonomy bar plot
│   ├── core-metrics-results/         # Core diversity metrics (.qza / .qzv)
│   │   ├── rarefied_table.qza
│   │   ├── shannon_vector.qza / observed_features_vector.qza
│   │   ├── faith_pd_vector.qza / evenness_vector.qza
│   │   ├── bray_curtis_distance_matrix.qza / jaccard_*.qza
│   │   ├── unweighted_unifrac_distance_matrix.qza / weighted_*.qza
│   │   ├── bray_curtis_pcoa_results.qza / jaccard_*.qza
│   │   ├── unweighted_unifrac_pcoa_results.qza / weighted_*.qza
│   │   ├── bray_curtis_emperor.qzv / jaccard_*.qzv
│   │   └── unweighted_unifrac_emperor.qzv / weighted_*.qzv
│   └── export/                       # Exported data (read by R scripts)
│       ├── feature-table.tsv / .biom   # ASV abundance table
│       ├── taxonomy.tsv                # Species annotation
│       ├── dna-sequences.fasta         # Representative sequences
│       ├── rarefied_table.tsv / .biom  # Rarefied ASV table
│       ├── alpha-diversity.tsv         # Alpha diversity metrics
│       ├── stats.tsv                   # DADA2 denoising stats
│       ├── tree.nwk                    # Phylogenetic tree (Newick)
│       ├── bray_curtis_distance_matrix.tsv / jaccard_*.tsv
│       ├── unweighted_unifrac_distance_matrix.tsv / weighted_*.tsv
│       ├── bray_curtis_pcoa_results.txt / jaccard_*.txt
│       ├── unweighted_unifrac_pcoa_results.txt / weighted_*.txt
│       ├── alpha/                     # Alpha diversity boxplots †
│       │   ├── alpha_diversity_boxplot.pdf     # Combined 3×2 grid
│       │   ├── alpha_boxplot_shannon.pdf       # Individual plots
│       │   ├── alpha_boxplot_chao1.pdf
│       │   ├── alpha_boxplot_simpson.pdf
│       │   ├── alpha_boxplot_pielou_evenness.pdf
│       │   ├── alpha_boxplot_observed_features.pdf
│       │   └── rarefaction_curves.pdf
│       ├── beta/                      # Beta diversity PCoA plots †
│       │   ├── beta_diversity_pcoa_bray_curtis.pdf
│       │   ├── beta_diversity_pcoa_jaccard.pdf
│       │   ├── beta_diversity_pcoa_unweighted_unifrac.pdf
│       │   ├── beta_diversity_pcoa_weighted_unifrac.pdf
│       │   ├── *pcoa_coords.tsv / *pcoa_variance.tsv
│       │   └── *permanova.tsv
│       ├── taxa/                      # Phylum composition †
│       │   ├── phylum_stacked_barplot.pdf
│       │   └── phylum_abundance_barchart.pdf
│       ├── heatmap/                   # Genus-level abundance heatmap †
│       │   └── genus_heatmap.pdf
│       ├── feature_tables/            # Processed tables †
│       │   ├── taxonomy_processed.tsv
│       │   ├── feature_table_with_taxonomy.tsv
│       │   ├── genus_abundance.tsv
│       │   ├── phylum_relative_abundance.tsv
│       │   ├── alpha_diversity_metrics.tsv
│       │   └── alpha_diversity_statistics.tsv
│       ├── faprotax/                  # FAPROTAX functional prediction
│       │   ├── faprotax.txt / faprotax_report.txt
│       │   ├── faprotax_report.{clean,mat,func_otu,otu_func}
│       │   ├── taxonomy.tsv
│       │   ├── rarefied_table.biom / rarefied_tax.biom
│       │   └── faprotax_visualization/   # FAPROTAX plots † (see below)
│       └── picrust2/                  # PICRUSt2 functional prediction
│           ├── feature-table.tsv
│           ├── dna-sequences.fasta
│           ├── picrust2_visualization/   # PICRUSt2 plots † (see below)
│           └── out/                      # Full PICRUSt2 output
│               ├── EC.tsv / KO.tsv / METACYC.tsv
│               ├── EC_predicted.tsv.gz / KO_predicted.tsv.gz
│               ├── EC_metagenome_out/ / KO_metagenome_out/
│               ├── pathways_out/ / KEGG.Pathway*.raw.txt
│               ├── marker_predicted_and_nsti.tsv.gz
│               └── intermediate/          # Per-sample intermediate files
└── logs/                            # Pipeline run logs
    ├── dada2.log / classify.log / tree.log
    └── *visualization.log
```

> **†** Directories marked with `†` are generated by R visualization scripts (Steps 15–16 of the pipeline), not by the shell pipeline directly.

## R Visualization Outputs

### QIIME2 16S Visualization (`QIIME2_16S_visualization.R`)

Generates alpha/beta diversity, phylum composition, genus heatmap, and processed tables under `results/export/`. Run via:

```bash
Rscript QIIME2_16S_visualization_EN.R    # English version
Rscript QIIME2_16S_visualization.R       # Chinese version
```

The R script is automatically called in pipeline Step 15 — no manual intervention needed.

### PICRUSt2 Visualization (`picrust2_visualization.R`)

Functional prediction analysis across **5 levels**:

| Module | Output | Key Files |
|--------|--------|-----------|
| NSTI | Quality assessment | Boxplot + values + Dunn posthoc |
| KEGG L1/L2 | Hierarchy overview | Stacked barplot (L1) + heatmap (L2) |
| KEGG Pathway | Differential analysis | Composition, α/β diversity, DESeq2 LRT, volcano plots |
| KO | PCA + differential | PCA ordination, DESeq2 LRT, volcano plots |
| EC | Enzyme profiling | Class stacked barplot, DESeq2 LRT, volcano plots |

```bash
# Run after pipeline Step 16 completes:
Rscript picrust2_visualization_EN.R    # English
Rscript picrust2_visualization.R       # Chinese
```

### FAPROTAX Visualization (`faprotax_visualization.R`)

Ecological function profiling across **8 modules**:

| Module | Content |
|--------|---------|
| 01 | Global ecological cycle composition (stacked barplot) |
| 02 | Top function barplot |
| 03 | Cycle-level differential analysis (K-W + Dunn posthoc) |
| 04 | (Reserved) |
| 05 | Bray-Curtis PCoA + PERMANOVA + pairwise PERMANOVA |
| 06 | High-CV function heatmap |
| 07 | OTU contribution tracing per function |
| 08 | Summary report |

```bash
Rscript faprotax_visualization_EN.R    # English
Rscript faprotax_visualization.R       # Chinese
```

> **R package dependencies:** `ggplot2`, `tidyr`, `dplyr`, `readr`, `vegan`, `ape`, `pheatmap`, `RColorBrewer`, `DESeq2`, `ggrepel`, `ggsci`, `viridis`, `rstatix`, `ggpubr`, `FSA`, `reshape2`

## Example Data

- **Dataset:** *Arabidopsis thaliana* rhizosphere microbiome (16S rRNA, V3-V4 region)
- **18 samples** (WT/KO/OE groups, 6 replicates each)
- **Platform:** Illumina HiSeq 2500, PE250
- **Accession:** CRA002352

View `.qzv` files at [view.qiime2.org](https://view.qiime2.org).

## Notes

- DADA2 parameters (`--p-trunc-len-f`, `--p-trunc-len-r`, `--p-max-ee`) should be adjusted based on sequencing data quality
- `--p-sampling-depth` for core metrics should be determined from alpha rarefaction curve and `table.qzv`
- QIIME2 2025.7 DADA2 chimera method: supports `consensus` and `none` only (`pooled` removed in QIIME2 2025.4)

## License

[MIT](LICENSE)

---

## 中文说明

### 项目简介

基于 QIIME2 的 16S rRNA 扩增子分析流程，适用于双端测序数据。包含完整的质控、去噪、物种注释、多样性分析、结果导出和 R 语言出版级可视化。

### 分析流程

1. **质控与引物切除** — FastQC + cutadapt
2. **去噪** — DADA2（纠错、ASV 推断、嵌合体去除）
3. **物种注释** — SILVA-138 分类器
4. **系统发育树** — MAFFT 比对 + FastTree
5. **多样性分析** — Alpha 稀疏曲线、核心多样性指标（Faith PD、Shannon、UniFrac）、PCoA
6. **结果导出** — 导出为 TSV/FASTA/Biom 格式
7. **可视化** — R 脚本生成出版级图表

### 仓库结构

```
QIIME2-16S-Workflow/
├── scripts/                          # 全部分析脚本（见下方表格）
├── metadata.txt                      # 示例元数据文件
├── qiime2/                           # 示例 .qzv 可视化文件
├── results/export/                   # 示例 pipeline 输出（预览用）
├── Project_file_structure.log        # 自动生成的目录树
├── LICENSE
├── .gitignore
└── README.md
```

### 脚本一览

| 脚本 | 功能 |
|------|------|
| `qiime2-16s-pipeline.sh` | 主分析流程 — 质控、DADA2、物种注释、建树、多样性、导出 |
| `qiime2-16s-pipeline_install.sh` | 环境安装 — 创建 QIIME2 和 QC 工具的 conda 环境 |
| `QIIME2_16S_visualization.R` | 出版级可视化 — α/β 多样性、门水平、属水平热图 |
| `picrust2_visualization.R` | PICRUSt2 功能预测可视化 — NSTI、KEGG、KO、EC |
| `faprotax_visualization.R` | FAPROTAX 生态功能可视化 — 循环、PCoA、OTU 贡献 |

> **使用方法：** 将所需脚本复制到项目工作目录（与 `metadata.txt` 同级）。脚本默认使用 `results/export/` 路径，无需手动调整。

### 快速开始

1. **环境安装：** 打开 `scripts/qiime2-16s-pipeline_install.sh`，逐条手动执行命令（`conda activate` 需要交互式 shell）
2. **数据准备：** 将双端 FASTQ 文件放入 `seq/` 目录，命名格式 `seq/<样本名>_1.fq.gz` 和 `seq/<样本名>_2.fq.gz`；将 `metadata.txt` 放在工作目录下
3. **运行流程：** 编辑 `scripts/qiime2-16s-pipeline.sh` 中的参数，按步骤逐条执行

> 本流程设计为**手动分步执行**，每完成一步请检查输出质量后再继续。

### 运行环境

- Miniconda + 两个 conda 环境：`qc_preprocess`（FastQC、cutadapt、MultiQC）和 `qiime2-2025.7`
- SILVA-138.2 分类器（要求 qiime2 ≥ 2023 版本）

### 运行时目录结构

pipeline 运行完成后，工作目录的结构参见上方英文部分的 [Runtime Output Structure](#runtime-output-structure)。关键词说明：

- `.qza` — QIIME2 二进制数据文件
- `.qzv` — QIIME2 可视化文件（可在 view.qiime2.org 在线查看）
- `results/export/` — 导出数据，供 R 脚本读取
- `alpha/`、`beta/`、`taxa/`、`heatmap/`、`feature_tables/` — R 脚本生成的图表
- `faprotax/`、`picrust2/` — 功能预测输入/输出

### R 可视化产出概览

**QIIME2 基础可视化：** α 多样性箱线图（Shannon、Chao1、Simpson、Pielou Evenness、Observed features）+ β 多样性 PCoA（Bray-Curtis、Jaccard、Unweighted/Weighted UniFrac）+ 门水平堆叠图 + 属水平热图

**PICRUSt2 功能预测：** NSTI 质量评估、KEGG 通路层级（L1/L2）、Pathway/KO/EC 差异分析（DESeq2 LRT + 火山图）、KO PCA 降维

**FAPROTAX 生态功能：** 全局生态循环组成、循环级差异分析（K-W）、PCoA + PERMANOVA、高变异功能热图、OTU 贡献追溯

### 示例数据

- **数据集：** 拟南芥根际微生物组（16S rRNA，V3-V4 区）
- **18 个样本**（WT/KO/OE 三组，每组 6 个重复）
- **测序平台：** Illumina HiSeq 2500，PE250
- **数据编号：** CRA002352

### 注意事项

- DADA2 参数需根据实际测序质量调整
- 抽平深度（`--p-sampling-depth`）需根据 α 稀疏曲线确定
- QIIME2 2025.7 仅支持 `consensus` 和 `none` 嵌合体检测方法

### 开源协议

[MIT](LICENSE)
