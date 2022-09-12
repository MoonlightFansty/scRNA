# 上皮通路分析
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
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")

epi <- readRDS("epi.RDS")


# 1. hallmark50通路分析和progeny分析

library(readr)
broad_pws <- read_lines("./data/h.all.v6.2.symbols.gmt") %>%
  lapply(str_split, "\\t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)

epi <- AddModuleScore(object = epi, features = broad_pws, name = names(broad_pws))

# 简单可视化
FeaturePlot(epi,features = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB1',split.by = 'sample_id' )
### 当然也可以对imm_anno和后面的str
#imm_anno <- AddModuleScore(object = imm_anno, features = broad_pws, name = names(broad_pws))
#str_anno <- AddModuleScore(object = str_anno, features = broad_pws, name = names(broad_pws), nbin = 12)

# progeny signatures分析
### PROGENy主要是用来对一些经典的肿瘤通路进行分析
### 能做bulk和单细胞数据
### 运行很快
library(progeny)
epi <- progeny(epi, scale = F, organism="Human", top=500, perm=1, return_assay=T)
epi <- ScaleData(epi, assay = "progeny")

#imm_anno <- progeny(imm_anno, scale = F, organism="Human", top=500, perm=1, return_assay=T)
#imm_anno <- ScaleData(imm_anno, assay = "progeny")
#str_anno <- progeny(str_anno, scale = F, organism="Human", top=500, perm=1, return_assay=T)
#str_anno <- ScaleData(str_anno, assay = "progeny")


# 2.progeny病人异质性

### 提取我们认为的恶性上皮(InferCNV，具体看上一节)
### 之所以不用原文的(比例推断法，具体看上一节）是因为按照原文34病人没有什么肿瘤细胞，这显然是不科学的
epi_tumor <- subset(epi, subset = malignant%in% c('Malignant'))
epi_tumor <- ScaleData(epi_tumor)

progeny_scores <- as.data.frame(t(GetAssayData(epi_tumor, assay = "progeny", slot = "scale.data")))
progeny_scores$cell_id <- rownames(progeny_scores)
progeny_scores <- gather(progeny_scores, Pathway, Activity, -cell_id)

cells_clusters <- FetchData(epi_tumor, c("sample_id", "malignant")) %>% filter(str_detect(sample_id, "t"))
cells_clusters$cell_id <- rownames(cells_clusters)

progeny_scores <- inner_join(progeny_scores, cells_clusters)

summarized_progeny_scores <- progeny_scores %>% 
  group_by(Pathway, sample_id) %>% 
  summarise(avg = mean(Activity), std = sd(Activity)) %>%
  pivot_wider(id_cols = Pathway, names_from = sample_id, values_from = avg) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

heatmap.2(summarized_progeny_scores, trace = "none",density.info = "none", col = bluered(100))
### 完美且结果符合原文，异质性非常显著
ggsave2("progeny热图.pdf", path = "./result_3v3/analysis/epi_analysis/", width = 5, height = 5)


# 3. 有丝分裂活性评估
mitotic_activity <- FetchData(epi_tumor, c("tissue_type", "malignant", "Phase", "sample_id")) %>%
  mutate(sample_id = factor(sample_id, levels = c("p034t", "p033t", "p032t")))

ggplot(mitotic_activity, aes(x = sample_id, fill = Phase)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_manual(values = use_colors) +
  coord_flip()
ggsave2("有丝分裂活性评估.pdf", path = "./result_3v3/analysis/epi_analysis/", width = 5,height=3)