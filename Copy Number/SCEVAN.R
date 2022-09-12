library(devtools)
install_github("miccec/yaGST")
# install_local('miccec-yaGST-56227df.tar.gz')
install_github("AntonioDeFalco/SCEVAN")
# install_local('AntonioDeFalco-SCEVAN-9a5e000.tar.gz')

library(SCEVAN)
library(Seurat)
suppressPackageStartupMessages(library(ggtree)) 
# count_mtx 是原始的reads的counts格式的表达量矩阵
count_mtx = sce@assays$RNA@counts
count_mtx[1:4,1:4]
phe=sce@meta.data
head(phe)
colnames(phe)

count_mtx[1:4,1:4]
results <- pipelineCNA(count_mtx,
                       sample = "MGH106",
                       par_cores = 4)
head(results)
table(results$class) 
save(results,file = 'MGH106_results.RData')