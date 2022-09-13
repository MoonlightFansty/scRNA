# A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer

# Stromal cell diversity associated with immune evasion in human triple-negative breast cancer
# 2018 Profiling human breast epithelial cells using single cell RNA sequencing identifies cell diversity

# 2020 Aging-Associated Alterations in Mammary Epithelia and Stroma Revealed by Single-Cell RNA Sequencing
Myo=c("Krt17", "Krt14", "Krt5", "Acta2", "Myl9", "Mylk", "Myh11")
Lum=c("Krt19", "Krt18", "Krt8")
Hs=c("Prlr", "Cited1", "Pgr", "Prom1", "Esr1")  
AV=c("Mfge8", "Trf", "Csn3", "Wfdc18", "Elf5", "Ltf")
Lp=c("Kit", "Aldh1a3", "Cd14")
genes_to_check = list(
  Myo=Myo,
  Lum=Lum,
  Hs=Hs, 
  AV=AV,
  Lp=Lp ) 
genes_to_check = lapply(genes_to_check , str_to_upper)
p_all_markers=DotPlot(sce.all, 
                      features = genes_to_check,
                      scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
ggsave('markers_for_breast_by_celltyper.pdf')



