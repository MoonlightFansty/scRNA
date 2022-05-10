# 安装和加载所需包
BiocManager::install("Seurat") 
BiocManager::install("dplyr")
BiocManager::install("patchwork")
library(dplyr)
library(Seurat)
library(patchwork) 

# 导入示例数据
pbmc.data <- Read10X(data.dir = "cellranger/p018n/filtered_feature_bc_matrix")
# 创建Seurat对象
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "p018n", min.cells = 3, min.features = 200)
# 过滤检测少于200个基因的细胞（min.features = 200）和少于3个细胞检测出的基因（min.cells = 3）
pbmc

# #参数解释
# CreateSeuratObject(
#   counts, # 未标准化的数据，如原始计数或TPMs
#   project = "CreateSeuratObject", # 设置Seurat对象的项目名称
#   assay = "RNA", # 与初始输入数据对应的分析名称
#   names.field = 1, # 对于每个cell的初始标识类，从cell的名称中选择此字段。例如，如果cell在输入矩阵中被命名为BARCODE_CLUSTER_CELLTYPE，则设置名称。字段设置为3以将初始标识设置为CELLTYPE
#   names.delim = "_", # 对于每个cell的初始标识类，从cell的列名中选择此分隔符。例如，如果cell命名为bar - cluster - celltype，则将此设置为“-”，以便将cell名称分离到其组成部分中，以选择相关字段
#   meta.data = NULL, # 要添加到Seurat对象的其他单元级元数据。应该是data.frame，其中行是单元格名称，列是附加的元数据字段
#   ...
#   min.cells # 包含至少在这些细胞检测到的features
#   min.features # 包含至少检测到这些features的细胞
# )

# 查看这三个基因的前三十行矩阵
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# 向pbmc新增一列percent.mt数据
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 展示前5个细胞的QC指标
head(pbmc@meta.data, 5)

# 使用小提琴图可视化QC指标
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter通常用于可视化 feature-feature 相关性
# nCount_RNA 与percent.mt的相关性
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# nCount_RNA与nFeature_RNA的相关性
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 # 合并两图

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# 选取 2500 > nFeature_RNA >200 和percent.mt < 5的数据

# 数据标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 鉴定高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 查看最高变的10个基因
top10 <- head(VariableFeatures(pbmc), 10)

# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 数据缩放
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 线性降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# 可视化降维
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) # 1个PC 500个细胞
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) # 15个PC

# 确定数据的维度
# JackStraw
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

# Elbow plot
ElbowPlot(pbmc)

# 细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# dims = 1:10 即选取前10个主成分来分类细胞
# 查看前5个细胞的分类ID
head(Idents(pbmc), 5)

# 非线性降维（UMAP/tSNE）
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
# 显示在聚类标签
DimPlot(pbmc, reduction = "umap", label = TRUE)
# 使用TSNE聚类
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
# 显示在聚类标签
DimPlot(pbmc, reduction = "tsne", label = TRUE)

# 保存rds，用于后续分析
# saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# 找差异表达基因(聚类标志cluster biomarkers)
# cluster 1的标记基因
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# 找出区分cluster 5与cluster 0和cluster 3的所有标记
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

VlnPlot(pbmc, features = c('VWF', 'DEPP1'))

# you can plot raw counts as well
VlnPlot(pbmc, features = c('C1QB', 'RETN'), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

# 每个聚类前10个差异基因表达热图(如果小于10，则绘制所有标记)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# 加上注释
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
