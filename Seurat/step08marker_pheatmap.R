### 加载R包
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
library(stringr)
library(cowplot)
library(scales)
library(readr)
library(progeny)
library(gplots)
library(tibble)
library(grid)
library(rlang)

theme_set(theme_cowplot())

### 设计配色，细胞类型注释已经好了，所以根据的是上辑第五篇的注释结果
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
  Alveolar_Macrophages = "#6bAEd6",
  `Monocyte-derived macrophages`= "#fff500",
  Monocytes= "#FA9FB5",
  `Myeloid dendritic cells`= "#DD3497",
  Plasmacytoid_dendritic_cells= "#7A0177",
  `Conventional_T_cells`= "#c2e699",
  `Regulatory T cells`= "#006837",
  `CD8+_T_cells`= "#bcbddc",
  NK_cells= "#4a1486",
  B_cells= "#969696",
  Plasma_cells= "#636363")

# 差异基因
### 巨噬细胞差异基因

imm_macro <- subset(imm_anno, subset = cell_type_imm %in% c("Alveolar_Macrophages",
                                                            "Monocyte-derived macrophages"))

### 重scale
imm_macro <- ScaleData(imm_macro)

Idents(imm_macro) <- imm_macro@meta.data$cell_type_imm

macro_markers <- FindAllMarkers(imm_macro, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

# group_by和top_n组合
top_macro_markers <- macro_markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

#devtools::install_github("elliefewings/DoMultiBarHeatmap") 

#DoMultiBarHeatmap提供多层注释，但这个包安装较难，必须版本均不报warning错，量力而装
DoMultiBarHeatmap::DoMultiBarHeatmap(imm_macro, features = top_macro_markers$gene, group.by = "cell_type_imm", 
                                     additional.group.by = 'tissue_type',draw.lines = F) +scale_fill_viridis()

ggsave2("macro_marker_gene.pdf", path = "./result_3v3/", width = 10, height = 9)


### T细胞差异基因

#T细胞亚群
imm_T <- subset(imm_anno, subset = cell_type_imm %in% c("Conventional_T_cells",
                                                        "CD8+_T_cells",
                                                        "Regulatory T cells"))

imm_T <- ScaleData(imm_T)

Idents(imm_T) <- imm_T@meta.data$cell_type_imm

markers_T <- FindAllMarkers(imm_T, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

top_markers_T <- markers_T %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)

DoMultiBarHeatmap(imm_T, features = top_markers_T$gene, group.by = "cell_type_imm", additional.group.by = "tissue_type", draw.lines = F) +
  scale_fill_viridis()


ggsave2("T_cell_marker_gene.pdf", path = "./result_3v3/", width = 10, height = 9)