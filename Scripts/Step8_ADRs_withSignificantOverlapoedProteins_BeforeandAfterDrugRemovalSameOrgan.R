rm(list = ls())

ADRDP_Proteins = readRDS("Data/ADRDP_Proteins.rds")
ADRDP_Proteins_AfterDrugRemoval = readRDS("Data/ADRDP_Proteins_AfterDrugRemoval.rds")

N_proteins = matrix(0,length(ADRDP_Proteins_AfterDrugRemoval),3)
colnames(N_proteins) = c("N_before_remove", "N_after_remove", "N_intersection")
for(i in names(ADRDP_Proteins_AfterDrugRemoval)){
  
  ADRDP_Proteins_i = ADRDP_Proteins[[i]]
  N1 = length(ADRDP_Proteins_i)
  
  ADRDP_Proteins_AfterDrugRemoval_i = ADRDP_Proteins_AfterDrugRemoval[[i]]
  N2 = length(ADRDP_Proteins_AfterDrugRemoval_i)
  
  intersected_proteins = intersect(ADRDP_Proteins_i, ADRDP_Proteins_AfterDrugRemoval_i)
  N3 = length(intersected_proteins)
  N_proteins[r,] = c(N1, N2, N3)
}

rownames(N_proteins) = names(ADRDP_Proteins_AfterDrugRemoval)

Pval = c()
for(i in 1:nrow(N_proteins)){
  Nt = 12694
  N1 = N_proteins[i,1]
  N2 = N_proteins[i,2]
  K = N_proteins[i,3]
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  pval = fisher.test(mat, alternative = "greater")$p.value
  Pval = c(Pval, pval)
}
N_proteins = data.frame(N_proteins)
N_proteins$pval = Pval

qvalue = p.adjust(Pval, method = "BH")
N_proteins$qval = qvalue

N_proteins_sig = N_proteins[N_proteins$qval<0.05,]

ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan = rownames(N_proteins_sig)
saveRDS(ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan, "Data/ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan.rds")




