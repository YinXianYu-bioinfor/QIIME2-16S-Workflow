#!/bin/bash
##########################################################
# 环境安装脚本
# 安装 QIIME2 16S 扩增子分析所需环境
# 包含: qc_preprocess + qiime2-2023.2 + silva-138-99 分类器
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

## 3. 安装 qiime2-2023.2 环境

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

## 4. 下载物种分类器
mkdir -p classifiers

# 方案 A: 百度网盘 (需手动下载)
echo "手动下载 silva-138-99-nb-classifier.qza 并放入 classifiers/"
echo "下载链接: https://pan.baidu.com/s/1FmThEjT_m7M-Zig0M34jxQ?pwd=42vv"

# 方案 B: QIIME2 官方链接 (下载较慢)
# wget -c https://data.qiime2.org/classifiers/silva/silva-138-99-nb-classifier.qza \
#     -O classifiers/silva-138-99-nb-classifier.qza

# 方案 C: 国内镜像 (如可用)
# wget -c ftp://download.nmdc.cn/tools/amplicon/silva-138-99-nb-classifier.qza \
#     -O classifiers/silva-138-99-nb-classifier.qza
