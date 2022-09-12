# Habermann et al.
# https://www.biorxiv.org/content/10.1101/753806v1

habermann_epi <- c("ABCA3", "SFTPB", "SFTPC", "AGER", "PDPN",  "KRT5", "TRP63", "NGFR", "SCGB1A1", "MUC5B", "KRT17", "FOXJ1", "TMEM190", "CAPS", "CHGA", "CALCA", "ASCL1", "PTPRC", "EPCAM")

habermann_imm <- c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3", "KIT", "MKI67", "CDK1", "EPCAM")

habermann_oth <- c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")


# Epithelial genes according to Habermann et al.
# manual annotation of normal cell types & detection of immune cell contaminated clusters
DotPlot(subset(epi, subset = cluster_type == "Normal"), features = habermann_epi, group.by = "SCT_snn_res.1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave2(filename = 'DotPlot_Normal_epi_markers.pdf',path = "figure/Recluster/epithelial/annotation", width =10, height =15)
# manual detection of immune cell contaminated clusters
DotPlot(subset(epi, subset = cluster_type == "Tumor"), features = habermann_epi, group.by = "SCT_snn_res.1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave2(filename = 'DotPlot_Tumor_epi_markers.pdf',path = "figure/Recluster/epithelial/annotation", width =10, height =15)

# Immune genes according to Habermann et al.
DotPlot(imm, features = habermann_imm, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2(filename = 'DotPlot_imm_markers.pdf',path = "figure/Recluster/immune/annotation", width =10, height =15)

# for (i in seq_along(habermann_imm)) {
#  plotlist <- list()
#  plotlist[1] <- FeaturePlot(imm, features = habermann_imm[i], sort.cell = T, combine = F)
#  plotlist[2] <- VlnPlot(imm, features = habermann_imm[i], pt.size = 0, combine = F)
#  print(CombinePlots(plots = plotlist))
# }

# Stromal genes according to Habermann et al.
DotPlot(str, features = habermann_oth, group.by = "SCT_snn_res.1") + coord_flip()
ggsave2(filename = 'DotPlot_str_markers.pdf',path = "figure/Recluster/stromal/annotation", width =10, height =15)

# for (i in seq_along(habermann_oth)) {
#  plotlist <- list()
#  plotlist[1] <- FeaturePlot(str, features = habermann_oth[i], sort.cell = T, combine = F)
#  plotlist[2] <- VlnPlot(str, features = habermann_oth[i], pt.size = 0, combine = F)
#  print(CombinePlots(plots = plotlist, ncol = 3))
# }


# Travaglini et al.
# https://www.biorxiv.org/content/10.1101/742320v1

# load marker gene lists
sheets <- paste0("Cluster ", c(1:58))
sheets <- sheets[-43]

signaturelist <- list()

for (i in seq_along(sheets)) {
  a <- read_excel("data/media-3.xlsx", sheet = sheets[[i]])
  a <- filter(a, a$...2 > 0.7 & a$...4 < 0.3)
  signaturelist <- c(signaturelist,a[1])
  remove(a)
}

# generate list with names of module scores in seurat object
names_of_modulescores <- c()
for (i in seq_along(signaturelist)){
  names_of_modulescores <- c(names_of_modulescores, paste0("T_", names(signaturelist[i]), i))
}

names_of_modulescores <- gsub(names_of_modulescores, pattern = " ", replacement = ".", fixed = TRUE)
names_of_modulescores <- gsub(names_of_modulescores, pattern = "+", replacement = ".", fixed = TRUE)
names_of_modulescores <- gsub(names_of_modulescores, pattern = "/", replacement = ".", fixed = TRUE)

# names_of_modulescores_unfiltered <- c()
# for (i in seq_along(signaturelist)){
#  names_of_modulescores_unfiltered <- c(names_of_modulescores_unfiltered, paste0(names(signaturelist[i]), "_unfiltered", i))
# }
# names_of_modulescores_unfiltered <- gsub(names_of_modulescores_unfiltered, pattern = " ", replacement = ".", fixed = TRUE)
# names_of_modulescores_unfiltered <- gsub(names_of_modulescores_unfiltered, pattern = "+", replacement = ".", fixed = TRUE)
# names_of_modulescores_unfiltered <- gsub(names_of_modulescores_unfiltered, pattern = "/", replacement = ".", fixed = TRUE)
# signature_list_updated <- list()
# for (i in seq_along(sheets)) {
#  signature_list_updated[[i]] <- checkGeneSymbols(signature_list[[i]])
# }

# calculate module scores for different subsets
epi <- AddModuleScore(object = epi, features = signaturelist, name = paste0("T_", names(signaturelist)))
imm <- AddModuleScore(object = imm, features = signaturelist, name = paste0("T_", names(signaturelist)))
str <- AddModuleScore(object = str, features = signaturelist, nbin = 12 , name = paste0("T_", names(signaturelist)))


# Vieira Braga et al.
# https://www.nature.com/articles/s41591-019-0468-5

# load marker gene lists
teichmann_signatures_epi <- read.csv("data/Fig1_DE_Lung_atlas_epithelial.csv")
teichmann_signatures_epi$gene <- as.character(teichmann_signatures_epi$gene)
teichmann_signatures_epi$cluster <- as.factor(teichmann_signatures_epi$cluster)
teichmann_epi <- levels(teichmann_signatures_epi$cluster)

teichmann_signatures_imm <- read.csv("data/Fig2_DE_Lung_atlas_immune.csv")
teichmann_signatures_imm$gene <- as.character(teichmann_signatures_imm$gene)
teichmann_signatures_imm$cluster <- as.factor(teichmann_signatures_imm$cluster)
teichmann_imm <- levels(teichmann_signatures_imm$cluster)

signaturelist2 <- list()

for (i in seq_along(teichmann_epi)) {
  signaturelist2 <- c(signaturelist2, teichmann_signatures_epi %>% filter(cluster == teichmann_epi[i], pct.2 < 0.3, avg_logFC > 0.7) %>% select(gene))
}

for (i in seq_along(teichmann_imm)) {
  signaturelist2 <- c(signaturelist2, teichmann_signatures_imm %>% filter(cluster == teichmann_imm[i], pct.2 < 0.3, avg_logFC > 0.7) %>% select(gene))
}

names(signaturelist2) <- gsub(c(teichmann_epi, teichmann_imm), pattern = "_", replacement = " ")

# generate list with names of module scores in seurat object
names_of_modulescores2 <- c()
for (i in seq_along(signaturelist2)){
  names_of_modulescores2 <- c(names_of_modulescores2, paste0("VB_", names(signaturelist2[i]), i))
}
names_of_modulescores2 <- gsub(names_of_modulescores2, pattern = " ", replacement = ".", fixed = TRUE)

# calculate module scores for different subsets
epi <- AddModuleScore(epi, features = signaturelist2, name = paste0("VB_", names(signaturelist2)))
imm <- AddModuleScore(imm, features = signaturelist2, name = paste0("VB_", names(signaturelist2)))
str <- AddModuleScore(str, features = signaturelist2, nbin = 12, name = paste0("VB_", names(signaturelist2)))


# cell type annotation and subsetting
# epithelial
annotation_curated_epi <- read_excel("data/annotation/annotation_epi.xlsx")
epi_anno <- epi
new_ids_epi <- annotation_curated_epi$cell_type_epi
names(new_ids_epi) <- levels(epi_anno)
epi_anno <- RenameIdents(epi_anno, new_ids_epi)
epi_anno@meta.data$cell_type_epi <- Idents(epi_anno)

epi_anno <- subset(epi_anno, subset = cell_type_epi != "Immune_contamination")
epi_anno <- ScaleData(epi_anno)


# immune
annotation_curated_imm <- read_excel("data/annotation/annotation_imm.xlsx")
imm_anno <- imm
new_ids_imm <- annotation_curated_imm$cell_type_imm
names(new_ids_imm) <- levels(imm_anno)
imm_anno <- RenameIdents(imm_anno, new_ids_imm)
imm_anno@meta.data$cell_type_imm <- Idents(imm_anno)

imm_anno <- subset(imm_anno, subset = cell_type_imm != "Epithelial_contamination")
imm_anno <- ScaleData(imm_anno)


# stromal
annotation_curated_str <- read_excel("data/annotation/annotation_str.xlsx")
str_anno <- str
new_ids_str <- annotation_curated_str$cell_type_str
names(new_ids_str) <- levels(str_anno)
str_anno <- RenameIdents(str_anno, new_ids_str)
str_anno@meta.data$cell_type_str <- Idents(str_anno)

str_anno <- subset(str_anno, subset = cell_type_str != "Immune/Epithelial contamination")
str_anno <- ScaleData(str_anno)


# Travaglini et al.
names_of_modulescores_original <- names(signaturelist)

# Epithelial
epi_type <- FetchData(epi_anno, vars = c(names_of_modulescores))

for(i in seq_along(names_of_modulescores_original)) {
  colnames(epi_type)[i] <- names_of_modulescores_original[i]
}

epi_type %>%
  merge(FetchData(epi_anno, vars = c("cluster_type", "cell_type_epi")), by = 0) %>%
  filter(cluster_type == "Normal") %>%
  group_by(cell_type_epi) %>%
  pivot_longer(cols = names_of_modulescores_original, names_to = "T_cell_type") %>%
  group_by(cell_type_epi, T_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = T_cell_type, y = cell_type_epi, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5B_epithelial.pdf", path = "../results", width = 30, height = 15, units = "cm")

# Immune
imm_type <- FetchData(imm_anno, vars = c(names_of_modulescores))

for(i in seq_along(names_of_modulescores_original)) {
  colnames(imm_type)[i] <- names_of_modulescores_original[i]
}

imm_type %>%
  merge(FetchData(imm_anno, vars = "cell_type_imm"), by = 0) %>%
  group_by(cell_type_imm) %>%
  pivot_longer(cols = names_of_modulescores_original, names_to = "T_cell_type") %>%
  group_by(cell_type_imm, T_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = T_cell_type, y = cell_type_imm, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5B_immune.pdf", path = "../results", width = 30, height = 30, units = "cm")

# Stromal
str_type <- FetchData(str_anno, vars = c(names_of_modulescores))

for(i in seq_along(names_of_modulescores_original)) {
  colnames(str_type)[i] <- names_of_modulescores_original[i]
}

str_type %>%
  merge(FetchData(str_anno, vars = "cell_type_str"), by = 0) %>%
  group_by(cell_type_str) %>%
  pivot_longer(cols = names_of_modulescores_original, names_to = "T_cell_type") %>%
  group_by(cell_type_str, T_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = T_cell_type, y = cell_type_str, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5B_stromal.pdf", path = "../results", width = 30, height = 20, units = "cm")


# Vieira Braga et al.
names_of_modulescores_original2 <- names(signaturelist2)

# Epithelial
epi_type <- FetchData(epi_anno, vars = c(names_of_modulescores2))

for(i in seq_along(names_of_modulescores_original2)) {
  colnames(epi_type)[i] <- names_of_modulescores_original2[i]
}

epi_type %>%
  merge(FetchData(epi_anno, vars = c("cluster_type", "cell_type_epi")), by = 0) %>%
  filter(cluster_type == "Normal") %>%
  group_by(cell_type_epi) %>%
  pivot_longer(cols = names_of_modulescores_original2, names_to = "VB_cell_type") %>%
  group_by(cell_type_epi, VB_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = VB_cell_type, y = cell_type_epi, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5A_epithelial.pdf", path = "../results", width = 20, height = 10, units = "cm")

# Immune
imm_type <- FetchData(imm_anno, vars = c(names_of_modulescores2))

for(i in seq_along(names_of_modulescores_original2)) {
  colnames(imm_type)[i] <- names_of_modulescores_original2[i]
}

imm_type %>%
  merge(FetchData(imm_anno, vars = "cell_type_imm"), by = 0) %>%
  group_by(cell_type_imm) %>%
  pivot_longer(cols = names_of_modulescores_original2, names_to = "VB_cell_type") %>%
  group_by(cell_type_imm, VB_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = VB_cell_type, y = cell_type_imm, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5A_immune.pdf", path = "../results", width = 20, height = 30, units = "cm")

# Stromal
str_type <- FetchData(str_anno, vars = c(names_of_modulescores2))

for(i in seq_along(names_of_modulescores_original2)) {
  colnames(str_type)[i] <- names_of_modulescores_original2[i]
}

str_type %>%
  merge(FetchData(str_anno, vars = "cell_type_str"), by = 0) %>%
  group_by(cell_type_str) %>%
  pivot_longer(cols = names_of_modulescores_original2, names_to = "VB_cell_type") %>%
  group_by(cell_type_str, VB_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = VB_cell_type, y = cell_type_str, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5A_stromal.pdf", path = "../results", width = 20, height = 20, units = "cm")


saveRDS(epi_anno, file = "seurat_objects/epi_anno.RDS")
saveRDS(imm_anno, file = "seurat_objects/imm_anno.RDS")
saveRDS(str_anno, file = "seurat_objects/str_anno.RDS")
