# main marker
# 2022 Genome Biology  Decoding the multicellular ecosystem of vena caval tumor thrombus in clear cell renal cell carcinoma by single‐cell RNA sequencing

# 2019 PNAS The single-cell transcriptomic landscape of early human diabetic nephropathy
# 2021 NC Single cell transcriptional and chromatin accessibility profiling redefine cellular heterogeneity in the adult human kidney

cg1=c('SLC34A1','LRP2','HAVCR1','CFH','SLC12A1','SLC12A3','SLC8A1','AQP2','SLC26A7','SLC26A4','NPHS2','EMCN','PIEZO2','COL1A2','PTPRC')
p <- DotPlot(sce.all, features = unique(cg1),
             assay='RNA'  )  + coord_flip()

p 
ggsave('check_GSE151302_markers.pdf' )


cg2=c('CUBN','LRP2','SLC34A1','SLC5A12','SLC5A2','ALDOB','CFH','SLC12A1','SLC12A3','SLC12A2','SLC8A1','AQP2','AQP6','SLC26A4','ATP6V0D2','NPHS1','NPHS2','PECAM1','FLT1','PDGFRB','PTPRC')
p <- DotPlot(sce.all, features = unique(cg2),
             assay='RNA'  )  + coord_flip()

p 
ggsave('check_GSE131882_markers.pdf' )


