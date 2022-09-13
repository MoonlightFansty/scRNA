# main marker
# 2021 Clin Cancer Res  Spatially Distinct Reprogramming of the Tumor Microenvironment Based On Tumor Invasion in Diffuse-Type Gastric Cancers

# 2020 Dissecting transcriptional heterogeneity in primary gastric adenocarcinoma by single cell RNA sequencing
# 2019 Cell Reports Dissecting the Single-Cell Transcriptome Network Underlying Gastric Premalignant Lesions and Early Gastric Cancer

# antral basal gland mucous cells (GMCs, marked with MUC6 and PGC; Ep1),
# enterocytes (marked with FABP1 and FABP2; Ep3),
# pit mucous cells (PMC marked with TFF1 and MUC5AC; Ep4),
# chief cells (marked with PGA4; Ep6),
# enteroendocrine cells (marked with CHGA; Ep7; Supplementary Fig. S12A).
# Ep2 (“proliferative”) cells were characterized by the upregulation of TSPAN8 and S100A6, known as cancer promotion–related proliferative markers (37, 38).
# Ep5 cells showed the expression of PMC-related marker genes such as SOX4 and were annotated as “PMC-like”

library(ggplot2) 
genes_to_check = c('PTPRC', 
                   'MUC2' , 'ITLN1',
                   'FABP1' , 'APOA1',
                   'CEACAM5' , 'CEACAM6',
                   'EPCAM', 'KRT18', 'MUC1',
                   'MUC6' , 'TFF2',
                   'PGA4' , 'PGA3',
                   'MUC5AC' , 'TFF1','CHGA' , 'CHGB')
library(stringr)   
p_all_markers <- DotPlot(sce.all, features = genes_to_check )  + coord_flip()

p_all_markers
ggsave(plot=p_all_markers, 
       filename="check_gastric_marker.pdf")
# goblet cells (MUC2 and ITLN1) 
# enterocytes (FABP1 and APOA1),
#  tumor markers (CEACAM5 and CEACAM6
# epithelial cells (EPCAM, KRT18, and MUC1)
# antral basal gland mucous cells (GMCs, marked as MUC6 and TFF2), 
# pit mucous cells (PMCs, marked as MUC5AC and TFF1), 
# chief cells (PGA4 and PGA3), 
# enteroendocrine cells (CHGA and CHGB)

celltype[celltype$ClusterID %in% c( 0:5,7,9,11,12,14,15),2]='epi'  

celltype[celltype$ClusterID %in% c( 17),2]='goblet'    #  (MUC2 and ITLN1) 
celltype[celltype$ClusterID %in% c(11),2]='enterocytes'   # (FABP1 and APOA1)  
celltype[celltype$ClusterID %in% c( 7),2]='GMCs'  # antral basal gland mucous cells (
celltype[celltype$ClusterID %in% c(0,12),2]='PMCs'  # pit mucous cells
celltype[celltype$ClusterID %in% c(9,15),2]='enteroendocrine'  # (CHGA and CHGB) 

