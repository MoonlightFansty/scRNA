rm(list = ls())

library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)

# 读取病人信息
seu_obj <- readRDS('seurat_object/preprocessing/scRNA_main_annotated.RDS')
metatable <- read_excel("data/metadata/patients_metadata.xlsx")

# 从seurat对象里提取orig.ident(病人)
metadata <- FetchData(seu_obj, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample_id <- metadata$orig.ident
# dplyr中的left_join合并
metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
# 重新加行名
rownames(metadata) <- metadata$cell_id

# 通过AddMetaData直接补入15列
seu_obj <- AddMetaData(seu_obj, metadata = metadata)

# 细胞周期运算,简单过
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

score_cc <- function(seu_obj) {
  seu_obj <- CellCycleScoring(seu_obj, s.genes, g2m.genes)
  seu_obj@meta.data$CC.Diff <- seu_obj@meta.data$S.Score - seu_obj@meta.data$G2M.Score
  return(seu_obj)
}

seu_obj <- score_cc(seu_obj)

FeatureScatter(seu_obj, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
  coord_fixed(ratio = 1)
ggsave2("Cell Cycle Scoring.pdf", path = "figure/Preprocessing/QC/", width=7 ,height=6)
DimPlot(seu_obj, group.by = "Phase" ,reduction = "umap")
ggsave2("Cell Cycle.pdf", path = "figure/Preprocessing/QC/", width=7 ,height=6)

# 去除周期影响
# seu_obj <- ScaleData(seu_obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seu_obj))

saveRDS(seu_obj, file = "seurat_object/Preprocessing/all.RDS")


# 提取子集,重归一化
# Idents设置最主要的区分因子
Idents(seu_obj) <- seu_obj@meta.data$main_cell_type
epi <- subset(seu_obj, idents = "Epithelial")
imm <- subset(seu_obj, idents = "Immune")
str <- subset(seu_obj, idents = "Stromal")

# 因为总量变了，所以需要重新归一化，否则分布将被打破
epi <- ScaleData(epi)
imm <- ScaleData(imm)
str <- ScaleData(str)

saveRDS(epi, file = "seurat_object/Preprocessing/epi.RDS")
saveRDS(imm, file = "seurat_object/Preprocessing/imm.RDS")
saveRDS(str, file = "seurat_object/Preprocessing/str.RDS")

# 作图:不同注释方式的Umap
# 修改group.by
DimPlot(seu_obj, group.by = "tissue_type", cols = use_colors, pt.size = 0.1)
ggsave2("Tissue type.png", path = "figure/Preprocessing/", width=7 ,height=6)

# 区分病人但不区分N/T
DimPlot(seu_obj, group.by = "patient_id", cols = use_colors, pt.size = 0.1)
ggsave2("Patient.png", path = "figure/Preprocessing/",  width=7 ,height=6)

DimPlot(seu_obj, group.by = "main_cell_type", cols = use_colors, pt.size = 0.1)
ggsave2("Cell type.png", path = "figure/Preprocessing/",  width=7 ,height=6)

cell_types <- FetchData(seu_obj, vars = c("sample_id", "main_cell_type", "tissue_type")) %>% 
  mutate(main_cell_type = factor(main_cell_type, levels = c("Stromal", "Immune", "Epithelial"))) %>% 
  mutate(sample_id = factor(sample_id, levels = rev(c("p018t", "p019t", "p023t", "p024t", "p027t", "p028t", "p030t", "p031t", "p032t", "p033t", "p034t", "p018n", "p019n", "p027n", "p028n", "p029n", "p030n", "p031n", "p032n", "p033n", "p034n"))))

ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x = sample_id, fill = main_cell_type, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = use_colors) +
  coord_flip()
ggsave2("Cell distribution.pdf", path = "figure/Preprocessing/", width=5 ,height=5)

# 作图:各病人细胞分布图
# FetchData拿metadata数据出来
# 可以看到32号病人的上皮激增

# cell_types <- FetchData(seu_obj, vars = c("sample_id", "main_cell_type", "tissue_type"))
# ggplot(data = cell_types) + 
#   geom_bar(mapping = aes(x = sample_id, fill = main_cell_type, ), position = "fill", width = 0.75) +
#   scale_fill_manual(values = use_colors) +
#   coord_flip()
