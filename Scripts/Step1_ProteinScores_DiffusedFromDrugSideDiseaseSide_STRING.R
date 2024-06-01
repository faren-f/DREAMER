rm(list = ls())
library(igraph)

# Read network
drug_ADR = readRDS("Data/drug_ADR.rds")
drug_target = readRDS("Data/drug_target.rds")
disease_symptom = readRDS("Data/disease_symptom.rds")
disease_gene_protein = readRDS("Data/disease_gene_protein.rds")

# all ADRs
ADRs = unique(drug_ADR$meddra_id)

#Read STRING PPI
PPI = read.csv("Data/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

Scores = list()
for(i in ADRs){
  #drugs that are linked to ADR_i
  Drugs_ADR_i = drug_ADR[drug_ADR$meddra_id == i, 1]
  
  # proteins that are linked to drugs that are linked to ADR_i 
  protein_drugs_ADR_i = drug_target[drug_target$drugbank_id %in% Drugs_ADR_i, 2]
  v = V(PPI_Graph)
  P = intersect(names(v), protein_drugs_ADR_i)
  protein_drugs_ADR_i = protein_drugs_ADR_i[protein_drugs_ADR_i%in% P]
  protein_drugs_ADR_i = prop.table(table(protein_drugs_ADR_i))

  if (length(protein_drugs_ADR_i)!=0){
    
    initial_scores = rep(0, length(v))
    names(initial_scores) = names(v)
    initial_scores[names(protein_drugs_ADR_i)] = protein_drugs_ADR_i
    
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
    # Disease that are connected to phenotype_i
    Di = disease_symptom[disease_symptom$meddra_id == i,1]
    
    G_all = disease_gene_protein[disease_gene_protein$mondo_id %in% Di, 2]
    G = intersect(names(v), G_all)
    G_all = G_all[G_all%in% G]
    G_all = prop.table(table(G_all))
    
    if (nrow(G_all)!=0){
      
      initial_scores2 = rep(0, length(v))
      names(initial_scores2) = names(v)
      initial_scores2[names(G_all)] = G_all
      
      actual_score_DP_Protein = list(page_rank(
        PPI_Graph,
        directed = FALSE,
        damping = 0.7,
        personalized = initial_scores2)$vector)[[1]]
      
      Null_score_DP_Protein = matrix(0,length(actual_score_DP_Protein), 1000)
      for(i2 in 1:1000){
        initial_scores2 = sample(initial_scores2)
        null_score_DP_Protein = list(page_rank(
          PPI_Graph,
          directed = FALSE,
          damping = 0.7,
          personalized = initial_scores2)$vector)[[1]]
        
        Null_score_DP_Protein[,i2] = null_score_DP_Protein
      }
      Actual_scores = cbind(actual_score_ADR_Protein, actual_score_DP_Protein)
      
      
      Scores[[i]][["Actual_scores"]] = Actual_scores
      Scores[[i]][["Null_score_ADR_Protein"]] = Null_score_ADR_Protein
      Scores[[i]][["Null_score_DP_Protein"]] = Null_score_DP_Protein
      
    }else{
      Scores[[i]] = NA
    }
  }else{
    Scores[[i]] = NA
  }
}

saveRDS(scores, "Data/ProteinScores_DiffusedFromDrugSideDiseasesSide_STRING.rds")
