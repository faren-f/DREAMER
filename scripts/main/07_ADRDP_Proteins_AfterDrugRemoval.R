rm(list = ls())
library(MASS)
library(igraph)

Scores = readRDS("data/ProteinScores_DiffusedFromDrugSideDiseasesSide_STRING_AfterDrugRemoval.rds")

ADRs_afterOrganRemoval = names(Scores)
saveRDS(ADRs_afterOrganRemoval, "data/ADRs_afterOrganRemoval.rds")

PPI = read.csv("data/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

v = V(PPI_Graph)

ADRDP_proteins = list()
N = c()
for(s in names(Scores)){
  
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
  
  for(i in rownames(actualScore)){
    
    # P-value calculation 
    rate = fitdistr(Null_score_ADR_Protein[i,], "exponential")$estimate
    pval_ADR_proteins = c(pval_ADR_proteins, pexp(actualScore[i,1], rate, lower.tail = FALSE, log.p = FALSE))
    
    rate = fitdistr(Null_score_DP_Protein[i,], "exponential")$estimate
    pval_DP_proteins = c(pval_DP_proteins, pexp(actualScore[i,2], rate, lower.tail = FALSE, log.p = FALSE))
  }
  
  
  ## Correction
  qval_ADR_protein = p.adjust(pval_ADR_proteins, method = "BH")
  N1 = sum(qval_ADR_protein<0.05)

  qval_DP_protein = p.adjust(pval_DP_proteins, method = "BH")
  N2 = sum(qval_DP_protein<0.05)
  
  # And
  Sum_adj = sum(qval_ADR_protein<0.05 & qval_DP_protein<0.05)
  
  if(Sum_adj == 0){
    next
  }
  ADRDP_proteins[[s]] = rownames(Null_score_ADR_Protein)[which(qval_ADR_protein<0.05 & qval_DP_protein<0.05)]
  N = rbind(N, c(N1, N2, Sum_adj))
  
}
colnames(N) = c("adj_Pr_ADR_Protein", "adj_Pr_DP_Protein", "adj_Pr_Overlap")

saveRDS(N,"data/N_proteins_after_Drug_Removal.rds")
saveRDS(ADRDP_proteins, "data/ADRDP_Proteins_AfterDrugRemoval.rds")



