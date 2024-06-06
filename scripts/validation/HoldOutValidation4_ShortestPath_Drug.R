rm(list=ls())

library(igraph)
library(parallel)

no_cores = detectCores()
cl = makeCluster(no_cores-2)

Positive_Negative_DrugsDiseases = readRDS("data/Positive_Negative_DrugsDiseases_Holdout_Analysis.rds")
ADRDP_Proteins_Holdout = readRDS("data/ADRDP_Proteins_Holdout.rds")
ADRs = names(ADRDP_Proteins_Holdout)

drug_target = readRDS("data/drug_target.rds")
drug_target$entrez_id_drug = as.character(drug_target$entrez_id_drug)

disease_gene_protein = readRDS("data/disease_gene_protein.rds")
disease_gene_protein$entrez_id_disaese = as.character(disease_gene_protein$entrez_id_disaese)

PPI = read.csv("data/knowledge_graph/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

clusterExport(cl, c("Positive_Negative_DrugsDiseases", "ADRDP_Proteins_Holdout",
                    "drug_target", "disease_gene_protein", "PPI_Graph"))
clusterEvalQ(cl,c(library(igraph)))

ParLoop = function(i){
  
  Drugs_val = Positive_Negative_DrugsDiseases[[i]][["Drugs_val"]]
  shortest_path_for_all_drugs = c()
  for(j in Drugs_val){
    targets_Drugs_val_j = unique(drug_target[drug_target$drugbank_id %in% j, 2])
    ADRDP_Proteins_Holdout_ADR_i = ADRDP_Proteins_Holdout[[i]]
    
    shortest_path = c()
    for(l in targets_Drugs_val_j){
      length_sp = sapply(ADRDP_Proteins_Holdout_ADR_i, function(k) {
        path = shortest_paths(PPI_Graph, k, to = l, output = "vpath")
        path = length(path$vpath[[1]])
        path = ifelse(path>0, path-1, 7)
        return(path)
      })
      shortest_path = c(shortest_path, min(length_sp))
    }
    shortest_path_for_all_drugs = c(shortest_path_for_all_drugs, min(shortest_path))
  }
  
  ##############
  #Negative Samples
  Drugs_negativeSamples = Positive_Negative_DrugsDiseases[[i]][["Drugs_negativeSamples"]]
  shortest_path_for_NegSample_all_drugs_all_repeats = c()
  
  
  for(r in 1:nrow(Drugs_negativeSamples)){
    Drugs_negativeSamples_r = Drugs_negativeSamples[r,]
    
    shortest_path_for_NegSample_all_drugs = c()
    for(j2 in Drugs_negativeSamples_r){
      targets_Drugs_Negsample_j = unique(drug_target[drug_target$drugbank_id %in% j2, 2])
      ADRDP_Proteins_Holdout_ADR_i = ADRDP_Proteins_Holdout[[i]]
      
      shortest_path_NS = c()
      for(l in targets_Drugs_Negsample_j){
        length_sp = sapply(ADRDP_Proteins_Holdout_ADR_i, function(k) {
          path = shortest_paths(PPI_Graph, k, to = l, output = "vpath")
          path = length(path$vpath[[1]])
          path = ifelse(path>0, path-1, 7)
          return(path)
        })
        shortest_path_NS = c(shortest_path_NS, min(length_sp))
      }
      shortest_path_for_NegSample_all_drugs = c(shortest_path_for_NegSample_all_drugs, min(shortest_path_NS))
    }
    shortest_path_for_NegSample_all_drugs_all_repeats = c(shortest_path_for_NegSample_all_drugs_all_repeats, shortest_path_for_NegSample_all_drugs)
  }
  shortest_path = list()
  shortest_path[["pos"]] = shortest_path_for_all_drugs
  shortest_path[["neg"]] = shortest_path_for_NegSample_all_drugs_all_repeats
  
  return(shortest_path)
  
}

ShortestPath = parLapply(cl, sapply(ADRs, list), ParLoop) 

saveRDS(ShortestPath, "data/Shortest_path_drugs.rds")

stopCluster(cl)

