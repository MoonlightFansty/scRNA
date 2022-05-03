# scRNA
## 单细胞上游分析
### 前言
本部分介绍在Linux环境下10X单细胞上游分析环境构建，主要从软件安装和上游分析示例讲起。上游分析示例包括两个内容：
* 常规的SRA文件下载，转fastq格式，然后走cellranger流程
* 作者提供bam文件，先把bam转fastq，然后走cellranger流程

### 一、构建10x上游分析环境
**（1）conda创建分析环境及软件安装** \
```
conda create -n scRNA python=3.8
conda activate scRNA
mamba install -y -c bioconda aspera-cli bwa samtools bedtools sambamba sra-tools bowtie2 samblaster fasterq-dump
```

**（2）CellRanger及参考基因组下载** \
2022年的最新版本是cellranger软件6.2.1，hg38和mm10参考基因组 \
这里cellranger的下载链接会失效，需要去官网获取：https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
