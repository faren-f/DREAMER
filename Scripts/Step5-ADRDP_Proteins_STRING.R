rm(list = ls())

N_proteins = readRDS("Data/N_proteins.rds")
ADRDP_Proteins_before_HGTest = readRDS("Data/ADRDP_Proteins_before_HGTest.rds")
rownames(N_proteins) = names(ADRDP_Proteins_before_HGTest)

# number of proteins in STRING PPI network
N_nodes = 12694

Pval = c()
for(i in rownames(N_proteins)){
  
  N1 = as.numeric(N_proteins[i,1])
  N2 = as.numeric(N_proteins[i,2])
  K = as.numeric(N_proteins[i,3])
  
  Nt = N_nodes
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  pval = fisher.test(mat, alternative = "greater")$p.value
  Pval = c(Pval, pval)
}
sum(Pval<0.05)
qval = p.adjust(Pval, method = "BH")
sum(qval<0.05)


significant_ADRs = rownames(N_proteins)[which(qval<0.05)]
qval_significant_ADRs = qval[which(qval<0.05)]

Pval = data.frame(Pval, qval)
rownames(Pval) = rownames(N_proteins)
Pval = Pval[Pval$qval<0.05,]

saveRDS(Pval, "Pval_ADRDP_Proteins.rds")
saveRDS(significant_ADRs, "Significant_ADRs.rds")
saveRDS(qval_significant_ADRs, "qval_significant_ADRs.rds")

ADRDP_Proteins = ADRDP_Proteins_before_HGTest[significant_ADRs]
saveRDS(ADRDP_Proteins, "ADRDP_Proteins.rds")



