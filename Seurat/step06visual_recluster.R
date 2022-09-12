#epithelial analyses
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
  p034 = "#FD6467",
  AT1 = "#2B8CBE",
  AT2 = "#045A8D",
  Club = "#006D2C",
  DifferentiatingCiliated = "#31A354",
  Ciliated = "#74C476",
  Neuroendocrine = "#8856A7",
  lepidic = "#0B775E",
  acinar = "#74A089",
  `mucinuous (papillary)` = "#E2D200",
  `(micro)papillary` = "#CEAB07",
  solid = "#B40F20",
  sarcomatoid = "#5B1A18",
  CNN = "chartreuse4",
  CNA = "orange")


#load data

#epi_anno <- readRDS("seurat_objects/epi_anno.RDS")

epi_anno@meta.data$cell_type_epi <- factor(epi_anno@meta.data$cell_type_epi, levels = c("AT2",
                                                                                        "AT1",
                                                                                        "Club",
                                                                                        "Ciliated",
                                                                                        "Neuroendocrine",
                                                                                        "Tumor"))

#add inferCNV clone scores (for inferCNV, use seurat object epi_anno and the following code: https://github.com/bischofp/single_cell_lung_adenocarcinoma, computation time ~12h)

scna_scores <- read.delim("../data/inferCNV_output/infercnv_clone_scores_nsclc.tsv") %>% filter(tissue_type == "Tumor") %>% filter(!is.na(cna_clone))
rownames(scna_scores) <- scna_scores$cell_id
scna_scores <- scna_scores %>% select(cna_clone) %>% mutate(cna_clone = as.character(cna_clone))
epi_anno <- AddMetaData(epi_anno, scna_scores)

scna_data <- FetchData(epi_anno, c("tissue_type", "cna_clone"))
scna_data <- scna_data %>% mutate(cna_clone = ifelse(is.na(cna_clone), "CNN", cna_clone))
epi_anno <- AddMetaData(epi_anno, scna_data)



#some plots

DotPlot(subset(epi_anno, subset = "cluster_type", idents = "Normal"), features = c("ABCA3", "SFTPC", "AGER", "PDPN",  "KRT5", "TRP63", "NGFR", "SCGB1A1", "MUC5B", "FOXJ1", "TMEM190", "CHGA", "CALCA"), group.by = "cell_type_epi") + 
  coord_flip() + 
  scale_color_viridis()
#ggsave2("DotPlot_markergenes_epi_cell_type.emf", path = "../results", width = 11, height = 8, units = "cm")
ggsave2("Fig2B.png", path = "../results", width = 11, height = 8, units = "cm")


DimPlot(epi_anno, group.by = "cell_type_epi", cols = use_colors, pt.size = 0.5)
ggsave2("Fig2A_celltype.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(epi_anno, group.by = "patient_id", cols = use_colors, pt.size = 0.5)
ggsave2("Fig2A_patients.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(epi_anno, group.by = "tissue_type", cols = use_colors, pt.size = 0.5)
ggsave2("Fig2A_tissuetype.png", path = "../results", width = 15, height = 15, units = "cm")

epi_cell_counts <- FetchData(epi_anno, vars = c("tissue_type", "cell_type_epi", "cna_clone", "cluster_type")) %>%
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))

ggplot(data = epi_cell_counts, aes(x = tissue_type, fill = cell_type_epi)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  scale_y_reverse() +
  coord_flip()
ggsave2("Fig2A_barplot.pdf", path = "../results", width = 20, height = 5, units = "cm")

ggplot(data = epi_cell_counts, aes(x = tissue_type, fill = cna_clone)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  scale_y_reverse() +
  coord_flip()
ggsave2("SuppFig4_barplot_cna_clone.pdf", path = "../results", width = 20, height = 5, units = "cm")

ggplot(data = epi_cell_counts, aes(x = tissue_type, fill = cluster_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("cyan3", "darkorange2")) +
  coord_flip()
ggsave2("SuppFig4_barplot_cluster_type.pdf", path = "../results", width = 20, height = 5, units = "cm")

DimPlot(epi_anno, group.by = "cna_clone", cols = use_colors, pt.size = 0.5)
ggsave2("SuppFig4_umap_cna_clone.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(epi_anno, group.by = "cluster_type", cols = c("cyan3", "darkorange2"), pt.size = 0.5)
ggsave2("SuppFig4_umap_cluster_type.png", path = "../results", width = 15, height = 15, units = "cm")



# stromal analyses
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
  p034 = "#FD6467",
  Endothelial = "#FED976",
  Lymphaticendothelial = "salmon",
  Fibroblast = "#2166AC",
  Myofibroblast = "#5AAE61",
  Smoothmuscle = "#9970AB",
  Mesothelial = "#40004B")


#str_anno <- readRDS("seurat_objects/str_anno.RDS")

str_anno@meta.data$cell_type_str <- factor(str_anno@meta.data$cell_type_str, levels = c("Endothelial",
                                                                                        "Lymphaticendothelial",
                                                                                        "Fibroblast",
                                                                                        "Myofibroblast",
                                                                                        "Smoothmuscle",
                                                                                        "Mesothelial"))

DimPlot(str_anno, group.by = "tissue_type", cols = use_colors)
#ggsave2("DimPlot_str_Normal_Tumor.pdf", path = "output/fig3", width = 15, height = 15, units = "cm")

DimPlot(str_anno, group.by = "patient_id", cols = use_colors, pt.size = 0.5)
ggsave2("SuppFig1C_str_patients.pdf", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(str_anno, group.by = "cell_type_str", label = F, split.by = "tissue_type", cols = use_colors, pt.size = 0.5)
ggsave2("Fig3A_umap.pdf", path = "../results", width = 30, height = 15, units = "cm")

DotPlot(str_anno, features = c("WT1", "UPK3B", "MYH11", "PDGFRB", "ACTA2", "MYLK", "LUM", "PDGFRA", "CCL21", "PROX1", "PECAM1", "VWF"), group.by = "cell_type_str") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() + 
  scale_color_viridis()
ggsave2("Fig3B.pdf", path = "../results", width = 16, height = 12, units = "cm")

#immune analyses
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
  p034 = "#FD6467",
  Alveolar_Macrophages = "#6bAEd6",
  CD14_Macrophages= "#fff500",
  Monocytes= "#FA9FB5",
  Myeloid_Dendritic= "#DD3497",
  Plasmacytoid_Dendritic= "#7A0177",
  T_conv= "#c2e699",
  T_reg= "#006837",
  T_CD8= "#bcbddc",
  NK_cells= "#4a1486",
  B_cells= "#969696",
  Plasma= "#636363",
  Mast= "#252525")


#load data

#imm_anno <- readRDS("seurat_objects/imm_anno.RDS")

imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_imm, levels = c("Alveolar_Macrophages",
                                                                                         "CD14_Macrophages",
                                                                                         "Monocytes",
                                                                                         "Myeloid_Dendritic",
                                                                                         "Plasmacytoid_Dendritic",
                                                                                         "Mast",
                                                                                         "T_conv",
                                                                                         "T_reg",
                                                                                         "T_CD8",
                                                                                         "NK_cells",
                                                                                         "B_cells",
                                                                                         "Plasma"))


DimPlot(imm_anno, group.by = "tissue_type", cols = use_colors)
#ggsave2("DimPlot_imm_Normal_Tumor.pdf", path = "output/fig4", width = 15, height = 15, units = "cm")
#ggsave2("DimPlot_imm_Normal_Tumor.png", path = "output/fig4", width = 35, height = 15, units = "cm")

DimPlot(imm_anno, group.by = "patient_id", cols = use_colors, pt.size = 0.5)
#ggsave2("DimPlot_imm_patients.pdf", path = "output/fig4", width = 30, height = 15, units = "cm")
ggsave2("SuppFig1C_imm_patients.png", path = "figure/Recluster/immune/", width = 15, height = 15, units = "cm")

DimPlot(imm_anno, group.by = "cell_type_imm", split.by = "tissue_type", cols = use_colors, pt.size = 0.5)
#ggsave2("DimPlot_imm_celltype.pdf", path = "output/fig4", width = 35, height = 15, units = "cm")
ggsave2("Fig4A_umap.png", path = "figure/Recluster/immune/", width = 35, height = 15, units = "cm")

DotPlot(imm_anno, features = c("CD68", "LYZ", "FABP4", "MARCO", "LGMN", "CSF1R", "CD14", "S100A12", "FCN1", "CD1C", "FCER1A", "LILRA4", "IL3RA", "KIT", "GATA2", "CD3E", "CD4", "FOXP3", "IL2RA", "CD8A", "NKG7", "KLRD1", "MS4A1", "CD79A", "JCHAIN", "IGKC", "MKI67"), group.by = "cell_type_imm") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() + 
  scale_color_viridis()
ggsave2("Fig4B.pdf", path = "figure/Recluster/immune/", width = 20, height = 20, units = "cm")