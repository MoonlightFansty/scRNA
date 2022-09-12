#### 其他还有什么方法剥离出肿瘤上皮（癌变上皮）呢，拷贝数推断值得学习
## 加载R包
library(seurat)
library(ggplot2)
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(stringr)
library(cowplot)
library(scales)
library(tibble)
library(gplots)
library(RColorBrewer)

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467",
  CNN = "chartreuse4",
  CNA = "orange")

epi <- readRDS("epi.RDS")
epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)


epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  print(DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

### 我们依然选择比较大的resolution，手动注释没有关系
### 控制在20个cluster左右即可
Idents(epi) <- epi@meta.data$SCT_snn_res.1


meta_tbl <- as_tibble(FetchData(epi, c("cell_id", "sample_id", "main_cell_type", "tissue_type", "patient_id")))

#瞄一眼
table(meta_tbl$main_cell_type)
table(meta_tbl$tissue_type)

cols.use <- c(
  p032 = "#6fc48c",
  p033 = "#6fcbcd",
  p034 = "#35c5f3",
  Normal = "steelblue",
  Tumor= "red",
  CNA = "red",
  CNN = "grey"
)

pids <- sort(unique(epi$patient_id))
# 瞄一眼
pids_with_tumor <- FetchData(epi, c("tissue_type", "patient_id")) %>% as_tibble() %>% distinct(tissue_type, patient_id) %>% filter(tissue_type == "Tumor") %>% pull(patient_id) %>% sort



##  准备infercnv的输入文件---------------------------------

## 1.准备原始count矩阵(raw count matrix, 简称rcm)
rcm <- epi@assays$RNA@counts
dim(rcm)

## 2.准备细胞注释文件(annotation file ,简称af)
## 结果文件中，我们的malignant指的是肿瘤组织取的样
## normal指的是正常组织取的样
## 理论上，肿瘤组织取得样品有癌上皮和正常上皮两种可能
af <- FetchData(epi, c("main_cell_type", "tissue_type", "patient_id", "cell_id")) %>%
  # filter(cell_id %in% cell_subset) %>%
  mutate(cell_type_infercnv = ifelse(tissue_type != "Normal", paste0("malignant_", patient_id), paste0("normal_", patient_id)),
         cell_id = str_replace_all(cell_id, ":", "_")) %>% 
  select(cell_id, cell_type_infercnv) %>% 
  as_tibble
dir.create('infercnv')
write.table(af, file = "./infercnv/annotation.txt", sep = "\t", col.names = F, row.names = F)



### 3.准备基因坐标文件gene order file ,gof
### 下载网站:https://data.broadinstitute.org/Trinity/CTAT/cnv/
library(readr)
gof <- read_tsv("./infercnv//gencode_v21_gen_pos.complete.txt", col_names = c("gene", "chr", "start", "end"))
### 处理下，以前有|
gof <- separate(gof,col='gene',into ='gene',sep = '\\|')
gof <- gof[!duplicated(gof$gene),]
common_genes <- intersect(gof$gene, rownames(rcm))
### 把原来的换掉
write.table(gof,file ='./infercnv/gencode_v21_gen_pos.complete_modi.txt',col.names = F,sep = '\t',row.names = F,quote = F)


## 4. 指定参考分组，选择normal为参考组！tumor组相对于normal组比较（实际上是相对值）
baseline_cell_types <- sort(unique(af$cell_type_infercnv)[!str_detect(unique(af$cell_type_infercnv), "malignant")])


#BiocManager::install('infercnv',update = F,ask=F)
library(dendextend)
library(infercnv)

# 构建cnv对象
cnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = rcm,
  annotations_file = "./infercnv/annotation.txt",
  gene_order_file = "./infercnv//gencode_v21_gen_pos.complete_modi.txt",
  ref_group_names = baseline_cell_types
) 

### 运行infercnv，运行速度比较缓慢
## 会过滤到很多低表达基因，序列长度标准化，log标准化,平滑，重中心化，聚类等
cnv_obj <- infercnv::run(
  cnv_obj, cutoff=0.1, out_dir="./infercnv/", 
  cluster_by_groups=T, denoise=T, HMM=F, 
  num_threads=parallel::detectCores()
)

write_rds(cnv_obj, "./infercnv/infercnv_results.rds")



## 可以加载运行结果，节约时间
cnv_obj <- readRDS("./infercnv/infercnv_results.rds")

expr <- cnv_obj@expr.data

colnames(af)=c('cell_id','group')
head(af)

gene <- rownames(expr)
gof=as.data.frame(gof)
rownames(gof)=gof$gene
# 取交集
sub_gof <-  gof[intersect(gene,gof$gene),]
expr=expr[intersect(gene,gof$gene),]


### 针对CNV推断的结果聚类
set.seed(123456)
### epi有20个cluster,此处我们也设20个好了
kmeans.result <- kmeans(t(expr), 20)
km <- data.frame(km_cluster=kmeans.result$cluster)
km$cell_id=rownames(km)

#合并km和分组信息
km=km%>%inner_join(af,by="cell_id") 
# 按照聚类结果排序
km_ordered=km[order(km$km_cluster),]
rownames(km_ordered)=km_ordered$cell_id
km_ordered$cell_id=NULL
km_ordered$km_cluster=factor(km_ordered$km_cluster) #将km_group转换为因子
head(km_ordered)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## 准备画复杂热图
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=c(RColorBrewer::brewer.pal(8, "Dark2")[1:8], RColorBrewer::brewer.pal(8, "RdGy")[1:8],
          RColorBrewer::brewer.pal(8, "Blues")[1:4])
names(color_v)=as.character(1:20)
left_anno <- rowAnnotation(df = km_ordered,col=list(group=c("normal_p032"="#0072b1","normal_p033"="#0072b3",
                                                            "normal_p034"="#0072b5",'malignant_p032'='#bc3c25',
                                                            'malignant_p033'='#bc3c27','malignant_p034'='#bc3c29'),km_cluster=color_v))

pdf("complexheatmap.pdf",width = 20,height =15)
ht = Heatmap(t(expr)[rownames(km_ordered),], 
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_gof$chr, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)

draw(ht)
dev.off()


##依托于正常细胞作为参考，我们可以推断8，16是正常的细胞
## 也就是说，8，16虽然也出现在肿瘤组织中，但他可能是肿瘤组织的正常细胞
km_ordered$malignant= ifelse(km_ordered$km_cluster == 8 |km_ordered$km_cluster == 16, 'Normal','Malignant')
km_ordered_add=km_ordered[colnames(epi),]

epi$malignant=km_ordered_add$malignant


epi$km_cluster=km_ordered_add$km_cluster

## 根据inferCNV的result来画图
DimPlot(epi,reduction = 'umap',split.by = 'tissue_type',group.by = 'malignant')
ggsave2("是否恶性（是否存在CNV）.pdf", path = "./result_3v3/analysis/epi_analysis/", width = 10,height=5)

## 
DimPlot(epi,reduction = 'umap',split.by = 'tissue_type',group.by = 'km_cluster',label = T)
ggsave2("km-cluster.pdf", path = "./result_3v3/analysis/epi_analysis/", width = 10,height=5)

### 结果覆盖原来的epi文件，你需要保存！！！！！！！
saveRDS(epi,file ='epi.RDS')