rm(list = ls())
library(igraph)
library(jsonlite)

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-1)

##############################
setwd("~/Desktop/ADR_Mechanism_Prediction/R/DREAMER3-STRING_DrugBank/")
drug_se = readRDS("Data/ProcessedData/Final_Network/drug_se.rds")
drug_target = readRDS("Data/ProcessedData/Final_Network/drug_target_db.rds")
disease_symptom = readRDS("Data/ProcessedData/Final_Network/disease_symptom.rds")
disease_gene_protein = readRDS("Data/ProcessedData/Final_Network/disease_gene_protein.rds")
Adjusted_proteins_DREAMER = readRDS("Data/ProcessedData/Final_Result/adjustedProteins_after_correction.rds")
SMILES = read.csv("Data/ProcessedData/Final_Network/DrugBank_SMILES.csv")


se = names(Adjusted_proteins_DREAMER)

#STRING PPI
PPI = read.csv("Data/ProcessedData/Final_Network/PPI_entrez_STRING800.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

clusterExport(cl, c("drug_se", "disease_gene_protein", 
                    "drug_target", "PPI_Graph",
                    "disease_symptom", "SMILES"))
clusterEvalQ(cl,c(library(igraph)))

ParLoop = function(i){
  
  D = drug_se[drug_se$meddra_id == i, 1]
  Drug_genes = drug_target[drug_target$drugbank_id %in% D, 1:2]
  
  v = V(PPI_Graph)
  G_Dr = intersect(names(v), Drug_genes$entrez_id_drug)
  Drugs = unique(Drug_genes[Drug_genes$entrez_id_drug%in% G_Dr, 1])
  
  # Finding the SMILES of drugs for each ADR
  DrugBankid_SMILES = SMILES[SMILES$DrugBank.Id %in% Drugs,]
  
  # Disease
  # Di = disease_symptom[disease_symptom$meddra_id == i,1]
  # 
  # Disease_genes = disease_gene_protein[disease_gene_protein$mondo_id %in% Di, 1:2]
  # G_Di = intersect(names(v), Disease_genes$entrez_id_disaese)
  # Diseases = unique(Disease_genes[Disease_genes$entrez_id_disaese%in% G_Di, 1])
  
  #ADR_DrugDiseaseList = list()
  #ADR_DrugDiseaseList[["Drugs"]] = DrugBankid_SMILES
  #ADR_DrugDiseaseList[["Diseases"]] = Diseases

  return(DrugBankid_SMILES)
}

ADR_DrugDiseaseList = parLapply(cl, sapply(se, list), ParLoop) 

stopCluster(cl)

json_data = toJSON(ADR_DrugDiseaseList)
write(json_data, "Data/ProcessedData/Final_Network/DrugList_forEachADR.json")
ADR_DrugDiseaseList = fromJSON("Data/ProcessedData/Final_Network/DrugList_forEachADR.json")


