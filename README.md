# scRNA analysis
## Cellranger count
### 前言
本部分介绍在Linux环境下10X单细胞上游分析环境构建，主要从软件安装和上游分析示例讲起，上游分析示例包括两个内容：
* 常规的SRA文件下载，转fastq格式，然后走`cellranger`流程
* 作者提供bam文件，先把bam转fastq，然后走`cellranger`流程
