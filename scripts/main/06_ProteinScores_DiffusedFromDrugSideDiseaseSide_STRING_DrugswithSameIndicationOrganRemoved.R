# In this script, I want to find the scores of proteins that are diffused from drugSide and diseaseSide,
# but after removing the drugs that have the same organ indication as the given phenotype

rm(list = ls())
library(igraph)

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

drug_ADR = readRDS("data/preprocessed_graph/drug_ADR.rds")
drug_target = readRDS("data/preprocessed_graph/drug_target.rds")
disease_symptom = readRDS("data/preprocessed_graph/disease_symptom.rds")
disease_gene_protein = readRDS("data/preprocessed_graph/disease_gene_protein.rds")
ADR_Drugs_After_removal = readRDS("data/preprocessed_graph/ADR_Drugs_After_removal.rds")

PPI = read.csv("data/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)


clusterExport(cl, c("drug_ADR", "disease_gene_protein", "drug_target", 
                    "PPI_Graph", "disease_symptom", "ADR_Drugs_After_removal"))
clusterEvalQ(cl,c(library(igraph)))

ParLoop = function(i){
  
  ADR_Drugs_After_removal_i = ADR_Drugs_After_removal[[i]]
  drug_target_new = drug_target[drug_target$drugbank_id %in% ADR_Drugs_After_removal_i, 2]
  v = V(PPI_Graph)
  
  P = intersect(names(v), drug_target_new)
  drug_target_new = drug_target_new[drug_target_new %in% P]
  
  drug_target_new = prop.table(table(drug_target_new))
  
  if (length(drug_target_new)!=0){
    
    initial_scores = rep(0, length(v))
    names(initial_scores) = names(v)
    initial_scores[names(drug_target_new)] = drug_target_new
    
    actual_score_ADR_Protein = list(page_rank(
      PPI_Graph,
      directed = FALSE,
      damping = 0.7,
      personalized = initial_scores)$vector)[[1]]
    
    Null_score_ADR_Protein = matrix(0,length(actual_score_ADR_Protein), 1000)
    for(i1 in 1:1000){
      initial_scores = sample(initial_scores)
      null_score_ADR_Protein = list(page_rank(
        PPI_Graph,
        directed = FALSE,
        damping = 0.7,
        personalized = initial_scores)$vector)[[1]]
      Null_score_ADR_Protein[,i1] = null_score_ADR_Protein
    }
    # Disease
    Di = disease_symptom[disease_symptom$meddra_id == i,1]
    G_all = disease_gene_protein[disease_gene_protein$mondo_id %in% Di, 2]
    G = intersect(names(v), G_all)
    G_all = G_all[G_all%in% G]
    G_all = prop.table(table(G_all))
    
    if (nrow(G_all)!=0){
      
      initial_scores2 = rep(0, length(v))
      names(initial_scores2) = names(v)
      initial_scores2[names(G_all)] = G_all
      
      actual_score_DPProtein = list(page_rank(
        PPI_Graph,
        directed = FALSE,
        damping = 0.7,
        personalized = initial_scores2)$vector)[[1]]
      
      
      Null_score_DP_Protein = matrix(0,length(actual_score_DPProtein), 1000)
      for(i2 in 1:1000){
        initial_scores2 = sample(initial_scores2)
        null_score_DP_Protein = list(page_rank(
          PPI_Graph,
          directed = FALSE,
          damping = 0.7,
          personalized = initial_scores2)$vector)[[1]]
        
        Null_score_DP_Protein[,i2] = null_score_DP_Protein
      }
      Actual_scores = cbind(actual_score_ADR_Protein, actual_score_DPProtein)
      
      Scores = list()
      Scores[["Actual_scores"]] = Actual_scores
      Scores[["Null_score_ADR_Protein"]] = Null_score_ADR_Protein
      Scores[["Null_score_DP_Protein"]] = Null_score_DP_Protein
      
      return(Scores)
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}

ADRs = names(ADR_Drugs_After_removal)
scores = parLapply(cl, sapply(ADRs, list), ParLoop) 
saveRDS(scores, "data/ProteinScores_DiffusedFromDrugSideDiseasesSide_STRING_AfterDrugRemoval.rds")

stopCluster(cl)



