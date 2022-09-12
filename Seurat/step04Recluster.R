library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(readr)
library(stringr)
library(progeny)
library(scales)

theme_set(theme_cowplot())

# 设定配色
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")


# 加载数据
epi <- readRDS("seurat_object/Preprocessing/epi.RDS")
imm <- readRDS("seurat_object/Preprocessing/imm.RDS")
str <- readRDS("seurat_object/Preprocessing/str.RDS")


# 主要亚群重聚类
# 上皮重聚类
epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)
ggsave2("epi_Elbow.pdf", path = "figure/Recluster/epithelial/resolution/", width = 10, height = 5)

# for (i in c(10, 15, 20, 25)){
#   umaptest <- RunUMAP(epi, dims = 1:i, verbose = F)
#   print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, "dimensions")))
#   remove(umaptest)
# }

epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("resolution", i,".pdf"), path = "figure/Recluster/epithelial/resolution/", width = 10, height = 10, units = "cm")
}

# 我们依然选择比较大的resolution,手动注释没有关系
# 控制在20个cluster左右即可
Idents(epi) <- epi@meta.data$SCT_snn_res.1


# 免疫重聚类
imm <- RunPCA(imm)
ElbowPlot(imm,  ndims = 50)
ggsave2("imm_Elbow.pdf", path = "figure/Recluster/immune/resolution/", width = 10, height = 5)

# for (i in c(10, 15, 20, 25)){
#  umaptest <- RunUMAP(imm, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, " dimensions")))
#  remove(umaptest)
# }

imm <- RunUMAP(imm, dims = 1:20)
imm <- FindNeighbors(imm, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  imm <- FindClusters(imm, resolution = i)
  DimPlot(imm, reduction = "umap") + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("resolution", i,".pdf"), path = "figure/Recluster/immune/resolution/", width = 10, height = 10, units = "cm")
}

# 控制在20个cluster左右即可
Idents(imm) <- imm@meta.data$SCT_snn_res.0.5


# 基质重聚类
str <- RunPCA(str)
ElbowPlot(str, ndims = 50)
ggsave2("str_Elbow.pdf", path = "figure/Recluster/stromal/resolution/", width = 10, height = 5)

# for (i in c(5, 10, 15, 20, 25, 30)){
#  umaptest <- RunUMAP(str, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, " dimensions")))
#  print(DimPlot(umaptest, reduction = "umap", group.by = "tissue_type") + labs(title = paste0(i, " dimensions")))
#  remove(umaptest)
# }

str <- RunUMAP(str, dims = 1:20)
str <- FindNeighbors(str, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  str <- FindClusters(str, resolution = i)
  DimPlot(str, reduction = "umap") + labs(title = paste0("resolution: ", i))
  ggsave2(paste0("resolution", i,".pdf"), path = "figure/Recluster/stromal/resolution/", width = 10, height = 10, units = "cm")
}

# 细胞不多,可以少点
Idents(str) <- str@meta.data$SCT_snn_res.1


# 疾病特异性细胞亚群
# 正常和肿瘤的上皮对比
# 发现确实多了一群
DimPlot(epi, group.by = "SCT_snn_res.1", label = T, repel = T, split.by = "tissue_type")
ggsave2("normal-tumor.png", path = "figure/Recluster/epithelial/", width = 30, height = 15, units = "cm")

# 比较正常和肿瘤细胞比例柱形图
epi_clusters <- FetchData(epi, vars = c("SCT_snn_res.1", "tissue_type"))

count_tumor <- epi_clusters %>% filter(tissue_type == "Tumor") %>% count() %>% as.numeric()
count_normal <- epi_clusters %>% filter(tissue_type == "Normal") %>% count() %>% as.numeric()

# 每个cluster计数总数
epi_counts <- epi_clusters %>% group_by(tissue_type) %>% count(SCT_snn_res.1)

# 除以总count数,是朴素的比例计算
proportion_tumor <- epi_counts %>% filter(tissue_type == "Tumor") %>% mutate(proportion = n/count_tumor)
proportion_normal <- epi_counts %>% filter(tissue_type == "Normal") %>% mutate(proportion = n/count_normal)

# 组合,最后一行很有意思,认为哪群比例高就属于哪群
proportion_epi <- full_join(proportion_normal, proportion_tumor, by = "SCT_snn_res.1") %>% 
  mutate(proportion.x = ifelse(is.na(proportion.x), 0,  proportion.x)) %>%  
  mutate(proportion.y = ifelse(is.na(proportion.y), 0,  proportion.y)) %>%
  mutate(tissue_type.x = "Normal") %>%
  mutate(tissue_type.y = "Tumor") %>%
  mutate(cluster_type = ifelse(proportion.x > proportion.y, "Normal", "Tumor"))

# 每个细胞的详细归属cluster
cluster_type_data <- left_join(x = epi_clusters, y = proportion_epi, by = "SCT_snn_res.1")
rownames(cluster_type_data) <- rownames(epi_clusters)

# 加入metadata中
epi <- AddMetaData(epi, select(cluster_type_data, cluster_type))


# 柱状图
n1 <- select(proportion_epi, c(tissue_type.x, SCT_snn_res.1, proportion.x)) %>%
  mutate(tissue_type = tissue_type.x) %>% 
  mutate(proportion = proportion.x) %>%
  mutate(tissue_type.x = NULL) %>%
  mutate(proportion.x = NULL)
t1 <- select(proportion_epi, c(tissue_type.y, SCT_snn_res.1, proportion.y)) %>%
  mutate(tissue_type = tissue_type.y) %>% 
  mutate(proportion = proportion.y) %>%
  mutate(tissue_type.y = NULL) %>%
  mutate(proportion.y = NULL)

proportion_epi2 <- rbind(n1, t1)

ggplot(proportion_epi2, aes(fill = tissue_type, y = proportion, x = SCT_snn_res.1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = use_colors)

ggsave2("norma-tumor-ratio.pdf", path = "figure/Recluster/epithelial/", width = 40, height = 20, units = "cm")

# 根据下面的图,我们可以明确地识别出肿瘤上皮中哪些类群是新增加的
# 也就是异型增生的上皮细胞——癌细胞
DimPlot(epi, group.by = "cluster_type",split.by = 'tissue_type' ,cols = use_colors, pt.size = 0.1)
ggsave2("New_tumor_cells.png", path = "figure/Recluster/epithelial/", width = 8, height = 4)


# 详细细胞类型注释
# 先画一些landscape,了解一下再聚类的质量
# 我们需要了解到的是,不同细胞注释当然是要结合起来看的
# 另外,我们对免疫细胞加入一些阴性marker也很重要
# 髓系
myeloid_markers <- c("S100A12", "FCN1", "S100A8", "S100A9", "CD14", "CTSS","VCAN", "LYZ", 
                     "MARCO", "FCGR1A", "C1QA", "APOC1", "LGMN", "CTSB", "FCGR3A", 
                     "MAFB", "MAF", "CX3CR1", "ITGAM", "CSF1R",
                     "FABP4", "MCEMP1", 
                     "IL1B", "CXCL8", 
                     "APOE", "CD163", "C1QB", "C1QC", 
                     "FCER1A", "CD1C", "CLEC9A", 
                     "LILRA4", "CLEC4C", "JCHAIN", "IL3RA", "NRP1", 
                     "CLEC10A", "PTCRA", "CCR7", "LAMP3", 
                     "ITGAX", "CD68", "MKI67", "CDK1", "EPCAM")

# T细胞
tcell_nk_markers <- c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", 
                      "CD8A", "CCL5", "NCR1", "NKG7", "GNLY", "NCAM1", 
                      "KLRD1", "KLRB1", "CD69", "KLRG1", 
                      "MKI67", "CDK1", "EPCAM")

# B细胞
bcell_plasma_mast_markers <- c("MS4A1", "CD19", "CD79A", "JCHAIN", "IGHA1", 
                               "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", 
                               "IGKC", "IGLC2", "IGLC3", "CPA3", "KIT", "MS4A2", 
                               "GATA2",  "MKI67", "CDK1", "EPCAM")

# Dotplot是最关键的手动注释的工具
DotPlot(imm, features = myeloid_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_myeloid_markers.png", path = "figure/Recluster/immune/annotation/", width = 20, height = 30, units = "cm")

DotPlot(imm, features = tcell_nk_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_T_NK_markers.png", path = "figure/Recluster/immune/annotation/", width = 20, height = 20, units = "cm")

DotPlot(imm, features = bcell_plasma_mast_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_B_Plasma_markers.png", path = "figure/Recluster/immune/annotation/", width = 20, height = 20, units = "cm")

DimPlot(imm, group.by = "SCT_snn_res.0.5", label = T, split.by = "tissue_type")
ggsave2("DimPlot_immune_clusters.png", path = "figure/Recluster/immune/annotation/", width =10, height =5)

DimPlot(imm, group.by = "patient_id", cols = use_colors)
ggsave2("DimPlot_immune_patients.png", path = "figure/Recluster/immune/annotation/", width = 10, height = 5)

