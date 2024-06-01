rm(list = ls())
library(MASS)
library(igraph)


Scores = readRDS("Data/ProteinScores_DiffusedFromDrugSideDiseasesSide_Holdout_Analysis.rds")

PPI = read.csv("Data/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)
v = V(PPI_Graph)

ADRs = names(Scores)
ADRDP_Proteins_before_HG_test = list()
N = c()
for(s in ADRs){
  
  if(is.na(s)){
    next
  }
  
  Null_score_ADR_Protein = Scores[[s]][["Null_score_ADR_Protein"]]
  Null_score_DP_Protein = Scores[[s]][["Null_score_DP_Protein"]]
  actualScore = Scores[[s]][["Actual_scores"]]
  
  rownames(actualScore) = names(v)
  rownames(Null_score_ADR_Protein) = names(v)
  rownames(Null_score_DP_Protein) = names(v)
  
  pval_ADR_proteins = c()
  pval_DP_proteins = c()
  
  for(j in rownames(actualScore)){
    
    # P-value calculation 
    rate = fitdistr(Null_score_ADR_Protein[j,], "exponential")$estimate
    pval_ADR_proteins = c(pval_ADR_proteins, pexp(actualScore[j,1], rate, lower.tail = FALSE, log.p = FALSE))
    
    rate = fitdistr(Null_score_DP_Protein[j,], "exponential")$estimate
    pval_DP_proteins = c(pval_DP_proteins, pexp(actualScore[j,2], rate, lower.tail = FALSE, log.p = FALSE))
  }
  
  
  ## Correction
  qval_ADR_protein = p.adjust(pval_ADR_proteins, method = "BH")
  N1 = sum(qval_ADR_protein<0.05)
  
  qval_DP_protein = p.adjust(pval_DP_proteins, method = "BH")
  N2 = sum(qval_DP_protein<0.05)
  
  # And
  Sum_adj = sum(qval_ADR_protein < 0.05 & qval_DP_protein < 0.05)
  
  if(Sum_adj == 0){
    next
  }
  ADRDP_Proteins_before_HG_test[[s]] = rownames(Null_score_ADR_Protein)[which(qval_ADR_protein<0.05 & qval_DP_protein<0.05)]
  N = rbind(N, c(N1, N2, Sum_adj))
}
colnames(N) = c("N1", "N2", "K")

saveRDS(N, "Data/N_proteins.rds")
saveRDS(ADRDP_Proteins_before_HG_test, "Data/ADRDP_Proteins_before_HG_test.rds")

