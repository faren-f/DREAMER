rm(list = ls())

drug_ADR = readRDS("data/preprocessed_graph/drug_ADR.rds")
ADRDP_Proteins = readRDS("data/ADRDP_Proteins.rds")
indication_proteins_before_HGTest = readRDS("data/Indication_Proteins.rds")
Intersected_ADRs = intersect(names(indication_proteins_before_HGTest), names(ADRDP_Proteins))

Intersected_proteins = list()
Pval = c()
Intersects = c()
for(i in Intersected_ADRs){
  
  ADRDP_Proteins_i = ADRDP_Proteins[[i]]
  indication_proteins_before_HGTest_i = indication_proteins_before_HGTest[[i]]
  
  Intersected_proteins[[i]] = intersect(ADRDP_Proteins_i, indication_proteins_before_HGTest_i)
  
  N1 = length(ADRDP_Proteins_i)
  N2 = length(indication_proteins_before_HGTest_i)
  K = length(Intersected_proteins[[i]])
  Intersects = c(Intersects,K)
  
  Nt = 12694
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  pval = fisher.test(mat, alternative = "greater")$p.value
  Pval = c(Pval, pval)
  
}
qval = p.adjust(Pval, method = "BH")
significant_ADRs = Intersected_ADRs[which(qval<0.05)]

saveRDS(significant_ADRs, "data/ADRs_withSignificantOvellap_ADRDPandIndication_Proteins.rds")

significant_ADRName_Indication = unique(drug_ADR[drug_ADR$meddra_id %in% significant_ADRs, 3])
saveRDS(significant_ADRName_Indication, "data/significant_ADRName_Indication.rds")



