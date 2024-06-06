rm(list = ls())
library(igraph)

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

drug_ADR = readRDS("data/preprocessed_graph/drug_ADR.rds")
drug_target = readRDS("data/preprocessed_graph/drug_target.rds")
disease_symptom = readRDS("data/preprocessed_graph/disease_symptom.rds")
disease_gene_protein = readRDS("data/preprocessed_graph/disease_gene_protein.rds")
drug_indication = readRDS("data/knowledge_graph/Drug_Indication.rds")
adj_Pr = readRDS("data/adjustedProteins_after_correction.rds")

ADRs = unique(drug_ADR$meddra_id)


#Physical PPI
PPI = read.csv("data/knowledge_graph/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)


clusterExport(cl, c("drug_ADR", "disease_gene_protein", 
                    "drug_target", "PPI_Graph",
                    "disease_symptom", "drug_indication"))
clusterEvalQ(cl,c(library(igraph)))

ParLoop = function(i){
  
  v = V(PPI_Graph)
  # indication
  D = drug_ADR[drug_ADR$meddra_id == i, 1]
  Indications_i = unique(drug_indication[drug_indication$drugbank_id %in% D, 2])
  
  G_all = disease_gene_protein[disease_gene_protein$mondo_id %in% Indications_i, 2]
  G = intersect(names(v), G_all)
  
  G_all = G_all[G_all%in% G]
  G_all = prop.table(table(G_all))
  
  
  if (nrow(G_all)!=0){
    
    initial_scores = rep(0, length(v))
    names(initial_scores) = names(v)
    initial_scores[names(G_all)] = G_all
    
    actual_score_indication = list(page_rank(
      PPI_Graph,
      directed = FALSE,
      damping = 0.7,
      personalized = initial_scores)$vector)[[1]]
    
    
    Null_score_indication = matrix(0,length(actual_score_indication), 1000)
    for(i1 in 1:1000){
      initial_scores = sample(initial_scores)
      null_score_disease = list(page_rank(
        PPI_Graph,
        directed = FALSE,
        damping = 0.7,
        personalized = initial_scores)$vector)[[1]]
      
      Null_score_indication[,i1] = null_score_disease
    }
    
    Scores = list()
    Scores[["Actual_scores"]] = actual_score_indication
    Scores[["Null_score_indication"]] = Null_score_indication
    
    return(Scores)
  }else{
    return(NA)
  }
}

scores = parLapply(cl, sapply(ADRs, list), ParLoop) 
saveRDS(scores, "data/ProteinScores_DiffusedFromIndications.rds")

stopCluster(cl)



