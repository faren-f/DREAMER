rm(list = ls())

N_proteins = readRDS("data/N_proteins.rds")
ADRDP_Proteins_before_HG_test = readRDS("data/ADRDP_Proteins_before_HG_test.rds")
rownames(N_proteins) = names(ADRDP_Proteins_before_HG_test)
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

qval = p.adjust(Pval, method = "BH")

significant_ADRs = rownames(N_proteins)[which(qval<0.05)]
qval_significant_ADRs = qval[which(qval<0.05)]

saveRDS(significant_ADRs, "data/Significant_ADRs_Holdout.rds")
saveRDS(qval_significant_ADRs, "data/qval_significant_ADRs_Holdout.rds")

ADRDP_Proteins_Holdout = ADRDP_Proteins_before_HG_test[significant_ADRs]
saveRDS(ADRDP_Proteins_Holdout, "data/ADRDP_Proteins_Holdout.rds")

# Removing the score columns to make the process faster 

Scores = readRDS("data/ProteinScores_DiffusedFromDrugSideDiseasesSide_Holdout_Analysis.rds")

Positive_Negative_DrugsDiseases = list()
for(d in names(Scores)){
  Positive_Negative_DrugsDiseases[[d]] = Scores[[d]][4:7]
}

I = names(Positive_Negative_DrugsDiseases)[names(Positive_Negative_DrugsDiseases) %in% names(ADRDP_Proteins_Holdout)]
Positive_Negative_DrugsDiseases = Positive_Negative_DrugsDiseases[I]

saveRDS(Positive_Negative_DrugsDiseases, "data/Positive_Negative_DrugsDiseases_Holdout_Analysis.rds")

