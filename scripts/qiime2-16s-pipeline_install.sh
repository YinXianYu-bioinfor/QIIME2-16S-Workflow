#!/bin/bash
##########################################################
# 环境安装脚本
# 安装 QIIME2 16S 扩增子分析所需环境
# 包含: qc_preprocess + qiime2-2025.7 (默认) + qiime2-2023.2 (备选) + silva-138.2 分类器
##########################################################

## 1. Miniconda 安装 (如已安装可跳过)
# wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# bash Miniconda3-latest-Linux-x86_64.sh -b -f
# ~/miniconda3/condabin/conda init
# 重新打开终端后继续

## 2. 安装 qc_preprocess 环境 (FastQC + cutadapt + MultiQC)
conda create -n qc_preprocess -c bioconda -c conda-forge fastqc cutadapt multiqc -y

# 验证
conda activate qc_preprocess
fastqc --version
cutadapt --version
multiqc --version

## 3. 安装 qiime2-2025.7 环境 (推荐)

### 方法: Conda 在线安装
# 官网安装指南: https://library.qiime2.org/quickstart/qiime2
wget -c https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2025.7-py310-linux-conda.yml
conda env create -n qiime2-2025.7 --file qiime2-amplicon-2025.7-py310-linux-conda.yml

# 验证
conda activate qiime2-2025.7
qiime --version

## 4. 安装 qiime2-2023.2 环境 (备选)

### 方案 A: 预打包安装 (推荐, 速度快)
n=qiime2-2023.2
wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
mkdir -p ~/miniconda3/envs/${n}
tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
conda activate ${n}
conda unpack

### 方案 B: Conda 在线安装 (备选)
# wget -c https://data.qiime2.org/distro/amplicon/qiime2-2023.2-py310-linux-conda.yml
# conda env create -n qiime2-2023.2 --file qiime2-2023.2-py310-linux-conda.yml

# 验证
conda activate qiime2-2023.2
qiime --version

## 5. 下载物种分类器
mkdir -p classifiers

# ⚠ SILVA138.2 分类器仅兼容 qiime2 ≥ 2023 版本
#   若使用 qiime2-2023.2 (备选环境), 请选择方案 C (silva-138-99 分类器)

# 方案 A: 百度网盘 (需手动下载)
echo "手动下载 SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza 并放入 classifiers/"
echo "下载链接: https://pan.baidu.com/s/1yFOvSJXofc6H7AahHSJq2A?pwd=5gdz"

# 方案 B: SILVA 官方链接
wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2025.7/SSU/full-length/uniform/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza \
    -O classifiers/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza

# 方案 C: 国内镜像 (如可用)
# wget -c ftp://download.nmdc.cn/tools/amplicon/silva-138-99-nb-classifier.qza \
#     -O classifiers/silva-138-99-nb-classifier.qza
