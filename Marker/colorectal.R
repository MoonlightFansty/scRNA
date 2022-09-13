# main marker
# 2021 EMBO Mol Med  Mitogen-activated protein kinase activity drives cell trajectories in colorectal cancer

# 2021 Cell  Differential pre-malignant programs and microenvironment chart distinct paths to malignancy in human colorectal polyps

# 恶性的癌症细胞
# ASC – adenoma specific cells,
# SSC – serrated specific cells
# 正常的上皮细胞，normal epithelial cells:
# 吸收细胞（ABS，absorptive cell）
# 隐窝顶部结肠细胞（CT，crypt top colonocytes）
# 内分泌细胞（EE enteroendocfine cells）
# 杯状细胞 （GOB，goblet cells）
# 干细胞（STM，stem cells）
# 过渡放大细胞（TAC，transit amplifing cells）
# 肠道簇细胞（TUF，tuft cells）

library(ggplot2) 
genes_to_check = c( "OTOP2","MEIS1","KRT20","GUCA2A",
                    "ALDOB","MSLN","MUC5AC","AQP5","TACSTD2",
                    "FSCN1","TFF2","ANXA1","ANXA10","REG4","MUC17",
                    "S100P","GSDMB","GSDMD","L18","RELB","MDK","RARA",
                    "RXRA","AHR","AGRN","PDX1","CLDN2","CD44","AXIN2",
                    "RNF43","TGFBI","EPHB2","TEAD2","CDX2","LGR5",
                    "OLFM4","ASCL2","PCNA","MKI67","ATOH1",
                    "MUC2","TFF3","CHGA","NEUROD1","POU2F3","SOX9")

