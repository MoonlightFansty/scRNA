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


#载入数据
### 数据在阅读原文里有，当然更建议你自行从第一步制作下来
#saveRDS(imm_anno,file='imm_anno.RDS')
imm_anno <- readRDS("imm_anno.RDS")

imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_imm, levels = c("Alveolar_Macrophages","Monocyte-derived macrophages","Monocytes", "Myeloid dendritic cells","Plasmacytoid_dendritic_cells","Conventional_T_cells","Regulatory T cells","CD8+_T_cells", "NK_cells", "B_cells","Plasma_cells"))


# 分群分析############
# 免疫细胞分淋巴系和髓系

###淋巴细胞亚群
imm_lympho <- subset(imm_anno, subset = cell_type_imm %in% c("Conventional_T_cells",
                                                             "CD8+_T_cells",
                                                             "Regulatory T cells",
                                                             "NK_cells",
                                                             "B_cells",
                                                             "Plasma_cells"))
### 记住重新scale！
imm_lympho <- ScaleData(imm_lympho)


### 髓系亚群
imm_myelo <- subset(imm_anno, subset = cell_type_imm %in% c("Alveolar_Macrophages",
                                                            "Monocyte-derived macrophages",
                                                            "Monocytes",
                                                            "Myeloid dendritic cells",
                                                            "Plasmacytoid_dendritic_cells"))
### 记住重新scale
imm_myelo <- ScaleData(imm_myelo)

### FetchData快速获得附加信息
### 淋巴系
lympho_counts <- FetchData(imm_lympho, vars = c("tissue_type", "cell_type_imm", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))


# 细胞比例图
### count根据patient_id，和tissue_type来统计总数，相当于excel分类汇总了
lympho_counts_tbl <- lympho_counts %>%
  dplyr::count(cell_type_imm, patient_id, tissue_type)
write_csv(lympho_counts_tbl, path = "./result_3v3/analysis/lympho_counts_tbl.csv")


### 髓系也一样
myelo_counts <- FetchData(imm_myelo, vars = c("tissue_type", "cell_type_imm", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal"))) 

myelo_counts_tbl <- myelo_counts %>%
  dplyr::count(cell_type_imm, patient_id, tissue_type)
write_csv(myelo_counts_tbl, path = "./result_3v3/analysis/lmyelo_counts_tbl.csv")


ggplot(data = lympho_counts, aes(x = tissue_type, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("barplot_lymphoid.pdf", path = "./result_3v3/analysis/", width = 8, height = 2)


ggplot(data = myelo_counts, aes(x = tissue_type, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("barplot_myeloid.pdf", path = "./result_3v3/analysis/", width = 8, height = 2)


#病人的异质性可见一斑，这是只有细胞比例图才可以揭示的现象！
lympho_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("lymphoid_per_patient.pdf", path = "./result_3v3/analysis/", width = 8, height = 2.5)

myelo_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("barplot_myeloid_per_patient.pdf", path = "./result_3v3/analysis/", width = 8, height = 2.5)

lympho_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("normal_lymphoid.pdf", path = "./result_3v3/analysis/", width = 8, height =2.5)

myelo_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("normal_myeloid.pdf", path = "./result_3v3/analysis/", width =8, height =2.5)
