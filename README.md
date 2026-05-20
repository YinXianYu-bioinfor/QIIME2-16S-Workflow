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
  - `qiime2-2023.2`: QIIME2 amplicon distribution
- SILVA-138 classifier (`silva-138-99-nb-classifier.qza`)

## Quick Start

### 1. Environment Setup

```bash
bash qiime2-16s-pipeline_install.sh
```

### 2. Data Preparation

Place raw paired-end FASTQ files in `seq/` directory with naming format:
```
seq/<sample>_1.fq.gz
seq/<sample>_2.fq.gz
```

Place metadata file (`metadata.txt`) in the working directory (TSV format, first column = sample ID, must include Group column).

### 3. Run Pipeline

Edit working directory and parameters in `qiime2-16s-pipeline.sh`, then execute step by step:

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
├── qiime2-16s-pipeline.sh          # Main analysis pipeline
├── qiime2-16s-pipeline_install.sh  # Environment installation
├── examples/                       # Example analysis results
│   ├── metadata.txt                # Sample metadata reference
│   ├── manifest                    # Sample manifest reference
│   ├── export/                     # Final exported results (TSV/FASTA)
│   │   ├── feature-table.tsv       # ASV abundance table
│   │   ├── taxonomy.tsv            # Species annotation
│   │   ├── dna-sequences.fasta     # Representative sequences
│   │   ├── alpha-diversity.tsv     # Alpha diversity metrics
│   │   ├── distance-matrix.tsv     # Beta diversity distances
│   │   ├── ordination.txt          # PCoA coordinates
│   │   ├── stats.tsv               # Denoising statistics
│   │   └── tree.nwk                # Phylogenetic tree (Newick)
│   ├── qiime2/                     # QIIME2 visualizations (.qzv)
│   └── logs/                       # Pipeline run logs
├── .gitignore
└── README.md
```

## Example Data

Dataset: *Arabidopsis thaliana* rhizosphere microbiome (16S rRNA, V3-V4 region)
- 18 samples (WT/KO/OE groups, 6 replicates each)
- Platform: Illumina HiSeq 2500, PE250
- Published under CRA002352

View `.qzv` files online at: https://view.qiime2.org/

## Notes

- DADA2 parameters (`--p-trunc-len-f`, `--p-trunc-len-r`, `--p-max-ee`) should be adjusted based on your sequencing data quality
- `--p-sampling-depth` for core diversity metrics should be determined from alpha rarefaction curve and table.qzv
- The chimera removal uses `pooled` method (more stable than `consensus` for large datasets)

## License

[MIT](LICENSE)
