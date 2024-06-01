rm(list = ls())
library(MASS)
library(igraph)

Scores = readRDS("Data/ProteinScores_DiffusedFromIndicationSide_STRING.rds")
PPI = read.csv("Data/ProcessedData/Final_Network/PPI_STRING.csv")

PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)
v = V(PPI_Graph)

ADRs = names(Scores)

indication_proteins_before_HGTest = list()
N = c()
for(s in ADRs){
  
  if(is.na(s)){
    next
  }
  if(sum(is.na(Scores[[s]]))){
    next
  }
  
  Null_score_indication_protein = Scores[[s]][["Null_score_indication_protein"]]
  actualScore = Scores[[s]][["Actual_scores"]]
  actualScore = data.frame(actualScore)
  rownames(actualScore) = names(v)
  rownames(Null_score_indication_protein) = names(v)
  
  pval_indication_proteins = c()
  
  for(j in rownames(actualScore)){
    
    # P-value calculation: Parametric 
    rate = fitdistr(Null_score_indication_protein[j,], "exponential")$estimate
    pval_indication_proteins = c(pval_indication_proteins, pexp(actualScore[j,1], rate, lower.tail = FALSE, log.p = FALSE))
    
  }
  
  ## Correction
  qvalue_indication_protein = p.adjust(pval_indication_proteins, method = "BH")
  N = sum(qvalue_indication_protein<0.05)
  indication_proteins_before_HGTest[[s]] = rownames(Null_score_indication_protein)[which(qvalue_indication_protein<0.05)]
  
}

saveRDS(indication_proteins, "Data/Indication_Proteins.rds")





