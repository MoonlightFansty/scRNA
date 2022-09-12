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



imm_anno <- readRDS("imm_anno.RDS")

imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_imm, levels = c("Alveolar_Macrophages","Monocyte-derived macrophages","Monocytes","Myeloid dendritic cells","Plasmacytoid_dendritic_cells","Conventional_T_cells","Regulatory T cells","CD8+_T_cells","NK_cells",
                                                                                         "B_cells", "Plasma_cells"))


# 分析点之一：对巨噬细胞分析M1,M2表型通路
### genelist文件在data文件夹中
###导入M1
m1m2_pws <- read_lines("./data/CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.gmt") %>%
  lapply(str_split, "\\t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)

### 导入M2
m1m2_pws <- append(m1m2_pws, read_lines("./data/CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.gmt") %>%
                     lapply(str_split, "\\t") %>% 
                     unlist(recursive = F) %>% 
                     lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
                     unlist(recursive = F))

### 运用addmodulescore函数，加入富集分数到metadata中
imm_anno <- AddModuleScore(object = imm_anno, features = m1m2_pws, name = c("m1up", "m1dn"), nbin = 12)


### 比较通路，当然这个没有显著性
VlnPlot(imm_anno, features = c("m1up1", "m1dn2"), 
        group.by = "cell_type_imm", pt.size = 0, idents = c("Alveolar_Macrophages",
                                                            "Monocyte-derived macrophages"), cols = use_colors)
ggsave2(filename = 'M1M2.pdf',path = './result_3v3/analysis/',width = 8,height = 7)


### 从下图我们看到单核来源的巨噬是M1表型的（略高），这个也提示我们不是所有的肿瘤相关巨噬都是M2表型的！

### 如果你跑了别的通路，如HALLMARK通路
#VlnPlot(imm_anno, features = c("HALLMARK_INFLAMMATORY_RESPONSE31",
#                               "HALLMARK_ALLOGRAFT_REJECTION46",
#                               "HALLMARK_INTERFERON_GAMMA_RESPONSE19",
#                               "HALLMARK_TNFA_SIGNALING_VIA_NFKB1",
#                               "m1up1", 
#                               "m1dn2"), 
#        group.by = "cell_type_imm", pt.size = 0, ncol = 3, idents = c("Alveolar_Macrophages",)


# 想设置显著性怎么办呢？
### 设置比较组
### 
comparisons <- list(c("Alveolar_Macrophages",
                      "Monocyte-derived macrophages"))
library(ggpubr)

### 难点在于我们的得分在metadata而非表达矩阵中，所以需要把它添加到表达矩阵中，再利用Vlnpubr包装函数
### 运行包装函数
Vlnpubr <- function(seo, gene_signature, file_name, test_sign,group_my,label  = "p.signif"){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(seo , features = signature,
            pt.size = 0.1, 
            group.by =  group_my, 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = label ) + NoLegend()
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  #file_name <- paste0(file_name, "_r.png")
  #ggsave(file_name, width = 14, height = 8)
}

### 我们的思路是新建一个scRNA，避免被破坏
scRNA=imm_anno
scRNA@assays$SCT@scale.data=rbind(scRNA@assays$SCT@scale.data,(scRNA$m1up1))
scRNA@assays$SCT@scale.data=rbind(scRNA@assays$SCT@scale.data,(scRNA$m1dn2))

scRNA_macro <- subset(scRNA, subset = cell_type_imm %in% c("Alveolar_Macrophages",
                                                           "Monocyte-derived macrophages"))
scRNA_macro=ScaleData(scRNA_macro)

### 出图
Vlnpubr(scRNA_macro,gene_signature = c('m1dn2','m1up1'), test_sign = comparisons,group_my = 'cell_type_imm')
ggsave2(filename = 'signa_M1M2.pdf',path = './result_3v3/analysis/',width = 8,height = 7)

### 我们终于加上了显著性

#分析点2：T耗竭和杀伤力分析

### 读取耗竭Marker
T_exhausted <- read_excel("./data/CD8_T_cells_exhausted.xlsx", skip = 1)
### 毒性marker
cytotoxicity <- c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7")

T_cell_markers <- list(T_exhausted$GeneSymbol, cytotoxicity)

imm_T <- AddModuleScore(imm_T, features = T_cell_markers, name = c("exhaustion", "cytotoxicity"), nbin = 12)

VlnPlot(imm_T, features = c("cytotoxicity2", "exhaustion1"), pt.size = 0, group.by = "cell_type_imm", cols = use_colors, idents = c("Conventional_T_cells",
                                                                                                                                    "CD8+_T_cells",
                                                                                                                                    "Regulatory T cells"), 
        ncol = 1)

ggsave2("T_cell_exhaustion_cytotxicity.pdf", path = "./result_3v3/analysis/", width = 4, height = 6)