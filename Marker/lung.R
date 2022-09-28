# main marker
# 2021 Science Advances Decoding the multicellular ecosystem of lung adenocarcinoma manifested as pulmonary subsolid nodules by single-cell RNA sequencing


# alveolar type I cell (AT1; AGER+)
# alveolar type II cell (AT2; SFTPA1)
# secretory club cell (Club; SCGB1A1+)
# basal airway epithelial cells (Basal; KRT17+)
# ciliated airway epithelial cells (Ciliated; TPPP3+)

# Pre-activated antiviral innate immunity in the upper airways controls early SARS-CoV-2 infection in children
library(ggplot2) 
genes_to_check =  c("SPRR3","GDPD3","SPRR1A","SPRR2A","RARRES2","TMPRSS11E",
                    "ASCL3","CFTR","FOXI2","1SG20","FOXI1",
                    "SAA4","SAA2","EFHC1","CCDC153","CCDC113","SAA1","CDC20B","FOXJ1",
                    "MYCL","FOXN4","CCNO",
                    "PIGR","BP1","MUC5A","VMO1","SCGB3A1","CYP2A13","CYP2B6","SCGB1A1",
                    "BCAM","KRT1","RT5","P63")

mainmarkers <- c('AGER', 'CLIC5', 'PDPN', # AT1
                 'LPCAT1', 'NAPSA', 'PGC', 'SFTPA1', 'SFTPA2', 'SFTPB', 'SFTPC', 'SLC34A2', # AT2
                 'KRT17', 'KRT5', 'KRT6A', # Basal
                 'AKAP14', 'ALDH3B1', 'ANKRD66', 'C11orf88', 'C11orf97', 'DNAI1', # Cilia
                 'PIGR', 'SCGB1A1', 'SCGB3A1', # Club
                 'CDH5', 'CLDN5', 'RAMP2', # EC
                 'C1R', 'COL1A2', 'DCN', # Fib
                 'AZGP1', 'CPE', 'TUBB2B', # NE
                 # immune
                 'CD19', 'CD79A', 'MS4A1', # B
                 'CD2', 'CD3D', 'CD3E', 'CD3G', 'CD4', 'TRAC', # CD4
                 'CD8A', 'CD8B', 'GZMK', # CD8
                 'C1orf54', 'LGALS2', # DC
                 'CD300E', 'CXCL8', 'EREG', 'S100A12', # Gran
                 'KIT', 'MS4A2', 'PTGS1', 'RGS13', # Mast
                 'CD68', 'FCGR1A', 'ITGAX', # M$
                 'GNLY', 'NKG7', # NK
                 'FOXP3', 'IL2RA', 'TNFRSF4' # Tregs
                 )
