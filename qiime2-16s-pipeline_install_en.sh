#!/bin/bash
##########################################################
# Environment Installation Script
# Install environments for QIIME2 16S amplicon analysis
# Includes: qc_preprocess + qiime2-2025.7 (default) + qiime2-2023.2 (fallback) + silva-138.2 classifier
##########################################################

## 1. Miniconda installation (skip if already installed)
# wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh -b -f
# ~/miniconda3/condabin/conda init
# Restart terminal before proceeding

## 2. Install qc_preprocess environment (FastQC + cutadapt + MultiQC)
conda create -n qc_preprocess -c bioconda -c conda-forge fastqc cutadapt multiqc -y

# Verify
conda activate qc_preprocess
fastqc --version
cutadapt --version
multiqc --version

## 3. Install qiime2-2025.7 environment (recommended)

### Method: Conda online installation
### Official installation guide: https://library.qiime2.org/quickstart/qiime2
wget -c https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2025.7-py310-linux-conda.yml
conda env create -n qiime2-2025.7 --file qiime2-amplicon-2025.7-py310-linux-conda.yml

# Verify
conda activate qiime2-2025.7
qiime --version

## 4. Install qiime2-2023.2 environment (fallback)

### Option A: Pre-packaged installation (recommended, faster)
n=qiime2-2023.2
wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
mkdir -p ~/miniconda3/envs/${n}
tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
conda activate ${n}
conda unpack

### Option B: Conda online installation (fallback)
# wget -c https://data.qiime2.org/distro/amplicon/qiime2-2023.2-py310-linux-conda.yml
# conda env create -n qiime2-2023.2 --file qiime2-2023.2-py310-linux-conda.yml

# Verify
conda activate qiime2-2023.2
qiime --version

## 5. Download taxonomic classifier
mkdir -p classifiers

# ⚠ SILVA138.2 classifier is only compatible with qiime2 ≥ 2023.x
#   If using qiime2-2023.2 (fallback env), choose Option C (silva-138-99 classifier)

# Option A: Baidu Netdisk (manual download required)
echo "Manually download SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza and place it in classifiers/"
echo "Download link: https://pan.baidu.com/s/1yFOvSJXofc6H7AahHSJq2A?pwd=5gdz"

# Option B: SILVA official link
wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2025.7/SSU/full-length/uniform/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza \
    -O classifiers/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza

# Option C: Domestic mirror (if available)
# wget -c ftp://download.nmdc.cn/tools/amplicon/silva-138-99-nb-classifier.qza \
#     -O classifiers/silva-138-99-nb-classifier.qza
