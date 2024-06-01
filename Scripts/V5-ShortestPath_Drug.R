rm(list=ls())

library(igraph)
library(parallel)
cl = makeCluster(40)

Data = readRDS("Data/Validation/Data_DREAMER3_withoutScores.rds")
adj_Proteins = readRDS("Data/Validation/adjustedProteins_after_correction.rds")
se = names(adj_Proteins)

drug_target = readRDS("Data/drug_target_db.rds")
drug_target$entrez_id_drug = as.character(drug_target$entrez_id_drug)

disease_gene_protein = readRDS("Data/disease_gene_protein.rds")
disease_gene_protein$entrez_id_disaese = as.character(disease_gene_protein$entrez_id_disaese)

PPI = read.csv("Data/PPI_entrez_STRING800.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

clusterExport(cl, c("Data", "adj_Proteins",
                    "drug_target", "disease_gene_protein", "PPI_Graph"))
clusterEvalQ(cl,c(library(igraph)))

ParLoop = function(i){
  
  Drugs_test = Data[[i]][["Drugs_test"]]
  shortest_path_for_all_drugs = c()
  for(j in Drugs_test){
    targets_Drugs_test_j = unique(drug_target[drug_target$drugbank_id %in% j, 2])
    identified_genes = adj_Proteins[[i]]
    
    shortest_path = c()
    for(l in targets_Drugs_test_j){
      length_sp = sapply(identified_genes, function(k) {
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
  Drugs_negativeSamples = Data[[i]][["Drugs_negativeSamples"]]
  shortest_path_for_NegSample_all_drugs_all_repeats = c()
  #for(r in 1:nrow(Drugs_negativeSamples)){
  for(r in 1:1){
    Drugs_negativeSamples_r = Drugs_negativeSamples[r,]
    
    shortest_path_for_NegSample_all_drugs = c()
    for(j2 in Drugs_negativeSamples_r){
      targets_Drugs_Negsample_j = unique(drug_target[drug_target$drugbank_id %in% j2, 2])
      identified_genes = adj_Proteins[[i]]
      
      shortest_path_NS = c()
      for(l in targets_Drugs_Negsample_j){
        length_sp = sapply(identified_genes, function(k) {
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

ShortestPath = parLapply(cl, sapply(se, list), ParLoop) 

saveRDS(ShortestPath, "Data/Validation/Shortest_path.rds")

stopCluster(cl)

