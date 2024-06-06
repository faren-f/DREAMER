rm(list = ls())
library(MASS)
library(igraph)

Scores = readRDS("data/ProteinScores_DiffusedFromDrugSideDiseasesSide_STRING.rds")

PPI = read.csv("data/knowledge_graph/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

v = V(PPI_Graph)

Pvalues = list()
i = names(Scores)
for(i in names(Scores)){
  
  if(sum(is.na(Scores[[i]]))){
    next()
  }
  Null_score_ADR_Protein = Scores[[i]][["Null_score_ADR_Protein"]]
  Null_score_DP_Protein = Scores[[i]][["Null_score_DP_Protein"]]
  actualScore = Scores[[i]][["Actual_scores"]]
  
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
  qval_DP_protein = p.adjust(pval_DP_proteins, method = "BH")
  
  Pvalues[[i]][["ADR_Protein"]] = qval_ADR_protein
  Pvalues[[i]][["DP_Protein"]] = qval_DP_protein
}

saveRDS(Pvalues, "data/Pval_ADRProteins_DPProteins.rds")




