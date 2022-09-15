# Cellranger count

## 一、构建10x上游分析环境
**(1)conda创建分析环境及软件安装**
```
conda create -n scRNA python=3.8
conda activate scRNA
mamba install -y -c bioconda aspera-cli bwa samtools bedtools sambamba sra-tools bowtie2 samblaster fasterq-dump
```
**(2)CellRanger及参考基因组下载** \
2022年的最新版本是cellranger软件7.0.0，hg38和mm10参考基因组 \
这里`cellranger`的下载链接会失效，需要去官网获取：https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest \
如果下载失败或下载速度过慢，可以在windows环境下载，在上传至服务器
```
# cellranger软件 7.0.0
wget -O cellranger-7.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.0.tar.gz?Expires=1656706802&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTY3MDY4MDJ9fX1dfQ__&Signature=CbGIIFt1SKCrEmusxc55I6a6-F8QuTfSaMVDeFXXy2Bf5vpQegUNMaxPEaBNWtMJg1xZbTQLTjtBzMgDhYDQMnpRddZIlIBeYHVO4BeTkkntNtJqMfP98yd3PHluJgRHe-HfeMda~ZAhal2nNwHktwyjICMGqSvvUN9Ngj27WR2-34uX6z3oZQjHukd7~KAdK2IjEoNuegyl3a35pgo-ZiXGwX4VRj-ZlWiHaDiMwPa8XU4jvioRdBf15qS5RWG1rkC2Ysei1Ux6XWpdylASmtPAivzhEgpUKykzcfu39q~9MLefLgdA1vYH8U3j4xfxS4GJ8hzBnfOxskAb9G-ypA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

md5sum cellranger-7.0.0.tar.gz
# 30855cb96a097c9cab6b02bdb520423f cellranger-7.0.0.tar.gz

tar -xzvf cellranger-7.0.0.tar.gz
```
```
# Human reference (GRCh38) dataset required for Cell Ranger 人类参考基因组
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

md5sum refdata-gex-GRCh38-2020-A.tar.gz
# dfd654de39bff23917471e7fcc7a00cd refdata-gex-GRCh38-2020-A.tar.gz

tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```
```
# Mouse reference dataset required for Cell Ranger 小鼠参考基因组
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz

md5sum refdata-gex-mm10-2020-A.tar.gz
# 886eeddde8731ffb58552d0bb81f533d refdata-gex-mm10-2020-A.tar.gz

tar -xzvf refdata-gex-mm10-2020-A.tar.gz
```
```
# 配置环境变量：PATH的路径设置对应的绝对路径
source sourceme.bash
# ln -s /home/~/scRNA/cellranger-7.0.0/cellranger /home/~/miniconda3/envs/scRNA/bin

cellranger
```
![CellRanger](https://github.com/MoonlightFansty/scRNA/blob/main/Cellranger%20count/Figures/Cellranger-results.png)

## 二、单细胞SRR文件上游分析
创建分析环境的文件夹：
```
mkdir GEO/GSE117988
mkdir sra fastq cellranger
```
**(1)数据下载** \
示例数据GSE117988:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117988 \
Metadata是表型数据，Accession List提供了SRA ID \
![GSE117988-SRA](https://github.com/MoonlightFansty/scRNA/blob/main/Cellranger%20count/Figures/GSE117988-SRA.png)
```
cd GEO/GSE117988/sra

# Accession List
cat >GSE117988
SRR7722937
SRR7722938
SRR7722939
SRR7722940
SRR7722941
SRR7722942
# ctrl+D 结束文件

# 下载 SRA
cat GSE117988 | while read id; do (nohup prefetch $id &); done
cat GSE117988 | while read id; do(mv $id/*.sra ./); done
cat GSE117988 | while read id; do(rm -rf $id); done
```
**(2)SRA转fastq** \
常规的SRA转fastq文件，用的是fastq-dump软件，速度非常慢，4-5个小时才能处理完一个样本 \
这里用新办法fasterq-dump，2分钟完成一个样本
```
cd ../fastq

# while 循环
cat >fastq.sh
ls ../sra/SRR* | while read id; do (nohup fasterq-dump -O ./ --split-files -e 40 $id  --include-technical &); done

bash fastq.sh
```
```
# for 循环
for i in `ls ../sra/SRR*`
do
i=$i
echo "nohup fasterq-dump -O ./ --split-files -e 40 $i --include-technical &"
done >fastq.sh

bash fastq.sh
```
```
# 查看文件
ls -lh
```
![GSE117988-Fastq](https://github.com/MoonlightFansty/scRNA/blob/main/Cellranger%20count/Figures/GSE117988-Fastq.png)
```
# 压缩文件
ls SRR*fastq | while read id; do (nohup gzip $id &); done
```
**(3)CellRanger count流程** \
对10X的fq文件运行CellRanger的counts流程，先做一个测试： \
首先需要对fastq.gz文件改名字，SampleName_S1_L001_R1_001.fastq.gz
```
mv SRR7722937_1.fastq.gz SRR7722937_S1_L001_I1_001.fastq.gz
mv SRR7722937_2.fastq.gz SRR7722937_S1_L001_R1_001.fastq.gz
mv SRR7722937_3.fastq.gz SRR7722937_S1_L001_R2_001.fastq.gz
```
```
cd ../cellranger

# 运行cellranger
cellranger count 
--id=SRR7722937 
--transcriptome=/home/~/scRNA/refdata-gex-GRCh38-2020-A 
--fastqs=/home/~/scRNA/GEO/GSE117988/fastq 
--sample=SRR7722937 
--nosecondary 
--localmem=30
```
**可以使用shell脚本批量完成**
```
# 批量修改fastq文件名
cat ../sra/GSE117988 | while read i; do (nohup mv ${i}_1*.gz ${i}_S1_L001_I1_001.fastq.gz &); done
cat ../sra/GSE117988 | while read i; do (nohup mv ${i}_2*.gz ${i}_S1_L001_R1_001.fastq.gz &); done
cat ../sra/GSE117988 | while read i; do (nohup mv ${i}_3*.gz ${i}_S1_L001_R2_001.fastq.gz &); done
```
```
cd ../cellranger

# 批量运行cellranger
ref=/home/~/scRNA/refdata-gex-GRCh38-2020-A
cat ../sra/GSE117988 | while read id;
do
nohup cellranger count --id=$id \
--transcriptome=$ref \
--fastqs=/home/~/scRNA/GEO/GSE117988/fastq \
--sample=$id \
--nosecondary \
# --localcores=10 \ #设置核心数
--localmem=30 &
done
```
```
mkdir cellranger
cat ../sra/GSE117988 | while read id; do(mkdir ./cellranger/$id); done
cat ../sra/GSE117988 | while read id; do(mv $id/outs/raw_feature_bc_matrix/* ./cellranger/$id); done
```
```
cd cellranger
tar -zcvf cellranger.tar.gz SRR7722937 SRR7722938 SRR7722939 SRR7722940 SRR7722941 SRR7722942
```
**最主要的几个参数:**

* --id 指定输出文件夹的名字

* --transcriptome 指定参考基因组的路径

* --sample 指定需要处理的fastq文件的前缀

* --expect-cell 指定预期的细胞数目，默认参数是3000个

* --localcores 指定计算的核心数

* --mempercore 指定内存大小 GB

* --nosecondary 不需要进行降维聚类（因为后期会用R可视化）

* --chemistry，默认是自动识别chemistry，但是有些时候识别不出chemistry的时候，需要加入参数特别标明

使用cellranger count --help可查看更多参数

## 三、单细胞BAM文件上游分析
**bam转fastq，再走cellranger的流程非常耗费计算机资源和时间** \
**(1)下载bam数据：** \
示例数据来自：https://www.ebi.ac.uk/ena/browser/view/PRJNA727404?show=reads
```
cd GEO/GSE

# 构建下载list
cat >download_file
/vol1/run/SRR144/SRR14424777/740_possorted_genome_bam.bam
/vol1/run/SRR144/SRR14424778/725_possorted_genome_bam.bam
/vol1/run/SRR144/SRR14424779/713_possorted_genome_bam.bam
/vol1/run/SRR144/SRR14424780/119_possorted_genome_bam.bam
/vol1/run/SRR144/SRR14424781/114_possorted_genome_bam.bam
/vol1/run/SRR144/SRR14424782/106_possorted_genome_bam.bam
/vol1/run/SRR144/SRR14424783/104_possorted_genome_bam.bam
/vol1/run/SRR144/SRR14424784/095_possorted_genome_bam.bam

# 批量下载
nohup ascp -v -QT -l 300m -P33001 -k1 -i /home/data/ssy40/anaconda3/envs/10x/etc/asperaweb_id_dsa.openssh --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp --file-list download_file ./ & >download.log
```
```
ls -lh | cut -d" " -f 5-
 56G 5月   6  2021 095_possorted_genome_bam.bam
 23G 5月   9  2021 104_possorted_genome_bam.bam
 22G 5月  11  2021 106_possorted_genome_bam.bam
 51G 5月  22  2021 114_possorted_genome_bam.bam
 45G 5月   6  2021 119_possorted_genome_bam.bam
 22G 5月   8  2021 713_possorted_genome_bam.bam
 27G 5月   8  2021 725_possorted_genome_bam.bam
 53G 12月 11  2021 740_possorted_genome_bam.bam
449K 12月 10 13:45 wget-log
```
**(2)bam转fastq** \
参考官网流程：https://support.10xgenomics.com/docs/bamtofastq?src=pr&lss=none&cnm=&cid=NULL \
cellranger产生的bam文件里是带有barcode与UMI的，储存在tag标签里：
```
samtools view 104_possorted_genome_bam.bam | less -SN
samtools view 104_possorted_genome_bam.bam | head -3 | tr "\t" "\n" | cat -n
```
* CB、CR、CY表示barcode，一般是16个碱基

* UB、UR、UY表示UMI，一般是10个碱基

* R一般代表原始测序数据，Y代表质量分数，而B代表校正后的R，可能对应碱基质量分数太低等因素，一般来说R与B都是相同的

**cellranger bamtofastq:**
```
# 3.3 构建name list文件
cat >name.list
740
725
713
119
114
106
104
095
```
```
mkdir fastq_file
```
```
# 准备shell脚本
cat > bamtofastq.sh
cat name.list |while read id
do
cellranger bamtofastq --nthreads 30 --traceback ${id}_possorted_genome_bam.bam ./fastq_file/${id}
done

# 运行脚本
nohup bash bamtofastq.sh & 
```
```
# 必要时可批量kill任务
# ps -ef | grep bamtofastq | awk '{print $2}' | while read id;do kill $id;done  #批量Kill
```
**(3)cellranger count \
批量完成：**
```
# 看一下文件夹地址
find ~/Projects_gao/newhbvhc/fastq_file/*/*count*/ | grep [1-1000] | grep -v XX/bam 

# 输入剩余未完成的文件list
cat >other_file.list
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/095/
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/104/
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/106/
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/114/
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/119/
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/713/
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/725/
/home/data/ssy40/Projects_gao/newhbvhc/fastq_file/740/
```
```
# 批量shell代码
cat other_file.list |while read id
do
ref=/home/data/ssy40/cellrange_soft/10x_refernce/refdata-gex-GRCh38-2020-A
sample_name=${id:0-36:3}_out_file
echo "cellranger count --id=$sample_name \
--transcriptome=$ref \
--fastqs=$id \
--sample=bamtofastq \
--nosecondary \
--localmem=20 \
--localcores=30"
done>other_file.sh

# 运行
nohup bash other_file.sh>log_106.114.119.log 2>&1 &
```

## 四、Cellranger 结果
**一个样本产生的结果：** \
![Results](https://github.com/MoonlightFansty/scRNA/blob/main/Cellranger%20count/Figures/Cellranger-results.png) \
**最重要的结果在```outs/```文件夹下：** \
![outs]() \
**结果解读：**

* web_summary.html：必看，官方说明 summary HTML file ，包括许多QC指标，预估细胞数，比对率等

* metrics_summary.csv：CSV格式数据摘要，可以不看

* possorted_genome_bam.bam：比对文件，用于可视化比对的reads和重新创建FASTQ文件，可以不看

* possorted_genome_bam.bam.bai：索引文件

* ```filtered_gene_bc_matrices```：是重要的一个目录，下面又包含了```barcodes.tsv.gz、features.tsv.gz、matrix.mtx.gz```，是下游Seurat、Scater、Monocle等分析的输入文件，是经过Cell Ranger过滤后构建矩阵所需要的所有文件

* filtered_feature_bc_matrix.h5：过滤掉的barcode信息HDF5 format，可以不看

* raw_feature_bc_matrix：原始barcode信息，未过滤的可以用于构建矩阵的文件，可以不看

* raw_feature_bc_matrix.h5：原始barcode信息HDF5 format，可以不看

* analysis：数据分析目录，下面又包含聚类clustering（有graph-based & k-means）、差异分析diffexp、主成分线性降维分析pca、非线性降维tsne，因为我们自己会走Seurat流程，所以不用看

* molecule_info.h5：可用于整合多样本，使用```cellranger aggr```函数

* cloupe.cloupe：官方可视化工具Loupe Cell Browser 输入文件，无代码分析的情况下使用，会代码的同学通常用不到

