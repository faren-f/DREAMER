rm(list = ls())
#library(jsonlite)
library(igraph)
library(parallel)

no_cores = detectCores()
cl = makeCluster(2)

##############################
setwd("~/Desktop/ADR_Mechanism_Prediction/R/DREAMER3-STRING_DrugBank/")
drug_se = readRDS("Data/ProcessedData/Final_Network/drug_se.rds")
drug_target = readRDS("Data/ProcessedData/Final_Network/drug_target_db.rds")
disease_symptom = readRDS("Data/ProcessedData/Final_Network/disease_symptom.rds")
disease_gene_protein = readRDS("Data/ProcessedData/Final_Network/disease_gene_protein.rds")
#TrainTestSplits = fromJSON("Data/ProcessedData/Final_Network/TrainTestSplits_forEachADR.json")
TrainTestSplits = readRDS("Data/ProcessedData/Final_Network/TrainTestSplits_forEachADR.rds")


# setwd("/cosybio/project/Faren/DrugSiderAI/DREAMER3/")
# drug_se = readRDS("Data/drug_se.rds")
# drug_target = readRDS("Data/drug_target_db.rds")
# disease_symptom = readRDS("Data/disease_symptom.rds")
# disease_gene_protein = readRDS("Data/disease_gene_protein.rds")
# TrainTestSplits = fromJSON("Data/TrainTestSplits_forEachADR.json")

se = names(TrainTestSplits)

all_drugs = unique(drug_se$drugbank_id)
all_diseases = unique(disease_symptom$mondo_id)

#Physical PPI
PPI = read.csv("Data/ProcessedData/Final_Network/PPI_entrez_STRING800.csv")
#PPI = read.csv("Data/PPI_entrez_STRING800.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)


clusterExport(cl, c("drug_se", "disease_gene_protein", 
                    "drug_target", "PPI_Graph",
                    "disease_symptom", "all_drugs", 
                    "all_diseases", "TrainTestSplits"))
clusterEvalQ(cl,c(library(igraph)))

ParLoop = function(i){
  
  D = TrainTestSplits[[i]]$C1e
  Drugs_test = names(D)[D == "test"]
  Drugs_train = names(D)[D == "train"]
  
  rest_drugs = all_drugs[!all_drugs%in% names(D)]
  

  Drugs_negativeSamples = c()
  for(r1 in 1:100){
    rest_drugs = sample(rest_drugs)
    Drugs_negativeSamples = rbind(Drugs_negativeSamples, rest_drugs[1: length(Drugs_test)])
  }
  
  ####Training
  protein_drugs_se = drug_target[drug_target$drugbank_id %in% Drugs_train, 2]
  v = V(PPI_Graph)
  
  P = intersect(names(v), protein_drugs_se)
  protein_drugs_se = protein_drugs_se[protein_drugs_se%in% P]
  
  # in this loop I find the targets for each drug separately to see which 
  # protein targets are repeated in several drugs to give more score to them
  protein_drugs_se = prop.table(table(protein_drugs_se))
  #hist(protein_drugs_se)
  
  if (length(protein_drugs_se)!=0){
    
    initial_scores = rep(0, length(v))
    names(initial_scores) = names(v)
    initial_scores[names(protein_drugs_se)] = protein_drugs_se
    
    actual_score_drug = list(page_rank(
      PPI_Graph,
      directed = FALSE,
      damping = 0.7,
      personalized = initial_scores)$vector)[[1]]
    
    Null_score_drug = matrix(0,length(actual_score_drug), 1000)
    for(i1 in 1:1000){
      initial_scores = sample(initial_scores)
      null_score_drug = list(page_rank(
        PPI_Graph,
        directed = FALSE,
        damping = 0.7,
        personalized = initial_scores)$vector)[[1]]
      Null_score_drug[,i1] = null_score_drug
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
      
      actual_score_disease = list(page_rank(
        PPI_Graph,
        directed = FALSE,
        damping = 0.7,
        personalized = initial_scores2)$vector)[[1]]
      
      
      Null_score_disaese = matrix(0,length(actual_score_disease), 1000)
      for(i2 in 1:1000){
        initial_scores2 = sample(initial_scores2)
        null_score_disease = list(page_rank(
          PPI_Graph,
          directed = FALSE,
          damping = 0.7,
          personalized = initial_scores2)$vector)[[1]]
        
        Null_score_disaese[,i2] = null_score_disease
      }
      Actual_scores = cbind(actual_score_drug, actual_score_disease)
      
      Data = list()
      Data[["Actual_scores"]] = Actual_scores
      Data[["Null_score_drug"]] = Null_score_drug
      Data[["Null_score_disaese"]] = Null_score_disaese
      
      Data[["Drugs_test"]] = Drugs_test
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

Data = parLapply(cl, sapply(se[1:2], list), ParLoop) 
#saveRDS(Data, "Data/Validation/Data_DREAMER3.rds")

stopCluster(cl)
