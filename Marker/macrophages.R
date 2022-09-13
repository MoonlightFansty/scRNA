# Single cell characterization of the immune microenvironment of melanoma brain and leptomeningeal metastases
# 树突细胞细分亚群比较多，是cDC1，cDC2，cDC3，以及pDC，有一些文章里面也把上面的 cDC3 叫做是 mregDC

# 2020 NC  scRNA-seq of gastric tumor shows complex intercellular interaction with an alternative T cell exhaustion trajectory
# 单核就区分成为了：
# Mono_CD14 ， classical CD14+CD16-
# Mono_FCGR3A ，non-classical CD14-CD16+ monocytes,
# 巨噬细胞区分成为了：
# resident tissue macrophages (RTMs)  （高表达 THBS1 ）
# lipid-associated macrophages （高表达 APOE）

# Single-cell landscape of the ecosystem in early-relapse hepatocellular carcinoma

library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2', ## fibo 
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'Amy1' , 'Amy2a2', # Acinar_cells
                   'PECAM1', 'VWF',  ## endo 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )




