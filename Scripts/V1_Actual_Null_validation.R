# In order to do Holdout validation:
# In this script, I want to find the actual and null diffusion scores of all the proteins diffusing from drug-side and diseases-side
# as well as the names of the discovery and validation drugs

rm(list = ls())
library(igraph)

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

drug_ADR = readRDS("Data/drug_ADR.rds")
drug_target = readRDS("Data/drug_target.rds")
disease_symptom = readRDS("Data/disease_symptom.rds")
disease_gene_protein = readRDS("Data/disease_gene_protein.rds")
ADRDP_Proteins = readRDS("Data/ADRDP_Proteins.rds")

ADRs = names(ADRDP_Proteins)
all_drugs = unique(drug_ADR$drugbank_id)
all_diseases = unique(disease_symptom$mondo_id)


PPI = read.csv("Data/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)


clusterExport(cl, c("drug_ADR", "disease_gene_protein", 
                    "drug_target", "PPI_Graph",
                    "disease_symptom", "all_drugs", "all_diseases"))
clusterEvalQ(cl,c(library(igraph)))

ParLoop = function(i){
  
  # drugs connected to ADR_i
  D = drug_ADR[drug_ADR$meddra_id == i, 1]
  
  # 20% of the drugs that are connected to ADR_i are randomly selected for validation and 80% for discovery 
  N_test_Dr = floor(0.2*length(unique(D)))
  D = sample(D)
  Drugs_val = D[1: N_test_Dr]
  Drugs_discovery = D[(N_test_Dr+1):length(unique(D))]
  
  rest_drugs = all_drugs[!all_drugs%in% D]
  
  Drugs_negativeSamples = c()
  for(r1 in 1:100){
    rest_drugs = sample(rest_drugs)
    Drugs_negativeSamples = rbind(Drugs_negativeSamples, rest_drugs[1: N_test_Dr])
  }
  
  ####Discovery
  protein_drugs_ADR_i = drug_target[drug_target$drugbank_id %in% Drugs_discovery, 2]
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
    
    
    # Disease
    Di = disease_symptom[disease_symptom$meddra_id == i,1]
    N_test_Di = floor(0.2*length(unique(Di)))
    Di = sample(Di)
    Diseases_test = Di[1:N_test_Di]
    Diseases_train = Di[(N_test_Di+1):length(unique(Di))]
    
    rest_diseases = all_diseases[!all_diseases%in% Di]
    
    Diseases_negativeSamples = c()
    for(r2 in 1:100){
      rest_diseases = sample(rest_diseases)
      Diseases_negativeSamples = rbind(Diseases_negativeSamples, rest_diseases[1: N_test_Di])
    }
    
    ####Training
    G_all = disease_gene_protein[disease_gene_protein$mondo_id %in% Diseases_train, 2]
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
      
      Data = list()
      Data[["Actual_scores"]] = Actual_scores
      Data[["Null_score_ADR_Protein"]] = Null_score_ADR_Protein
      Data[["Null_score_DP_Protein"]] = Null_score_DP_Protein
      
      Data[["Drugs_val"]] = Drugs_val
      Data[["Diseases_test"]] = Diseases_test
      
      Data[["Drugs_negativeSamples"]] = Drugs_negativeSamples
      Data[["Diseases_negativeSamples"]] = Diseases_negativeSamples
      
      
      return(Data)
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}

Data = parLapply(cl, sapply(se, list), ParLoop) 
saveRDS(Data, "Data/ProteinScores_DiffusedFromDrugSideDiseasesSide_Holdout_Analysis.rds")

stopCluster(cl)


