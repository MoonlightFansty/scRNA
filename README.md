# scRNA analysis
## 单细胞上游分析
### 前言
本部分介绍在Linux环境下10X单细胞上游分析环境构建，主要从软件安装和上游分析示例讲起，上游分析示例包括两个内容：
* 常规的SRA文件下载，转fastq格式，然后走`cellranger`流程
* 作者提供bam文件，先把bam转fastq，然后走`cellranger`流程

### 一、构建10x上游分析环境
**(1)conda创建分析环境及软件安装**
```
1. conda create -n scRNA python=3.8
2. conda activate scRNA
3. mamba install -y -c bioconda aspera-cli bwa samtools bedtools sambamba sra-tools bowtie2 samblaster fasterq-dump
```
**(2)CellRanger及参考基因组下载** \
2022年的最新版本是cellranger软件6.2.1，hg38和mm10参考基因组 \
这里`cellranger`的下载链接会失效，需要去官网获取：https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
>如果下载失败或下载速度过慢，可以在windows环境下载，在上传至服务器
```
## 1.1 cellranger软件 6.2.1
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1651612516&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTE2MTI1MTZ9fX1dfQ__&Signature=mveKaux6HIjgV5uTm~0D7stn-fy-akGwbRr3ENUe24VM6WMy5JAZ7uoBAREuoG9XpQu4aUsFrGpXZ07o9Q5NvWzRGItyE2j2LGKwNZ-Lw9shhrJufpqg0QqGqftDvXmvSma2vu7hSGt43YIDfW-Yq6qK3G4u-f~XmpKfU77dkFiSEDpmHLKRTkZNF4xCffKvAV7WcKoaAhrd2z0hgB4SPLgcOfs1yu-7u84NsD7Epgb8NUlTdHgQIH1huHw5RMgUdZ6eoWhMMlMN0RDbizYWnEm1OU31uxheiGAxz0~~QxbP4rgkeD~1yQ0kqCLsqMagYwuUryR8~ptG1pXmNt1log__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

md5sum cellranger-6.1.2.tar.gz
# 310d4453acacf0eec52e76aded14024c cellranger-6.1.2.tar
tar -xzvf cellranger-6.1.2.tar.gz
```
```
## 1.2 Mouse reference dataset required for Cell Ranger 小鼠参考基因组
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

md5sum refdata-gex-mm10-2020-A.tar.gz
# 886eeddde8731ffb58552d0bb81f533d refdata-gex-mm10-2020-A.tar.gz
tar -xzvf refdata-gex-mm10-2020-A.tar.gz
```
```
## 1.3 Human reference (GRCh38) dataset required for Cell Ranger 人类参考基因组
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

md5sum refdata-gex-GRCh38-2020-A.tar.gz
# dfd654de39bff23917471e7fcc7a00cd refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```
```
## 1.4 配置环境变量：PATH的路径设置对应的绝对路径
ln -s /home/liuhonghao/scRNA/cellranger-6.1.2/cellranger /home/liuhonghao/miniconda3/envs/scRNA/bin

cellranger
```
