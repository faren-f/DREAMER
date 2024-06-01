rm(list = ls())

N_proteins = readRDS("Data/Validation/N_proteins.rds")
adj_Pr = readRDS("Data/Validation/adjustedProteins.rds")
rownames(N_proteins) = names(adj_Pr)
N_nodes = 12694

Pval = c()
for(i in rownames(N_proteins)){
  print(i)
  
  N1 = as.numeric(N_proteins[i,1])
  N2 = as.numeric(N_proteins[i,2])
  K = as.numeric(N_proteins[i,3])
  
  Nt = N_nodes
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  pval = fisher.test(mat, alternative = "greater")$p.value
  Pval = c(Pval, pval)
}
sum(Pval<0.05)
adj = p.adjust(Pval, method = "BH")
sum(adj<0.05)

adjusted_SEs = rownames(N_proteins)[which(adj<0.05)]
adjusted_SEs_pval_adj = adj[which(adj<0.05)]

saveRDS(adjusted_SEs, "Data/Validation/SENames_after_Correction.rds")
saveRDS(adjusted_SEs_pval_adj, "Data/Validation/pval_adjusted_SEs.rds")

adj_Pr_after_Correction = adj_Pr[adjusted_SEs]
saveRDS(adj_Pr_after_Correction, "Data/Validation/adjustedProteins_after_correction.rds")