# 2020  Single-cell transcriptomic architecture and intercellular crosstalk of human intrahepatic cholangiocarcinoma

# normal epithelial
# cholangiocytes (546 cells, 1.7%, marked with FYXD2, TM4SF4, and ANXA4);
# hepatocytes (328 cells, 1.0%, marked with APOC3, FABP1, and APOA1);

# A single cell atlas of the human liver tumor microenvironment

library(ggplot2) 
genes_to_check =  c('FXYD3','CLDN4','CEACAM6','CEACAM5','ELF',
                    'CLDN10','SLC22A10','FETUB','LBP','HPR','LECT2',
                    'SERPINA10','CD5L','VCAM1','CETP','LILRB5','MARCO',
                    'SDC3','TREM2','GPNMB','CAPG','FCER1A','CD1C','CLEC10A',
                    'JAML','S100A9','FCN1','S100A8','FGR','XCR1','CLEC9A',
                    'IDO1','WDFY4','FLT3','CPNE3','GZMA','CD3E',
                    'KLRB1','NKG7','CD7','CCL5','IGLL5','FCRL5',
                    'TNFRSF17','DERL3','JCHAIN','MZB1','CPE',
                    'SLCO2A1','CLEC14A','TGM2','PODXL','VWA1',
                    'SOX18','PLVAP','CD34','ICAM2','RELN','CLEC4M',
                    'CLEC1B','CLEC4G','FCN2','OIT3','RERGL', 
                    'MYH11','ITGA8','PLN','ADIRF','OLFML2A','NDUFA4L2',
                    'RGS5','TPPP3','PLXDC1','FRZB','MMP11','CTHRC1','INHBA',
                    'HOPX','POSTN','LTBP2','CXCL12','PTGDS','MASP1','FBLN1',
                    'C7','HGF','TPX2','MKi67','UBE2C','ASPM','TOP2A','RRM2' )
