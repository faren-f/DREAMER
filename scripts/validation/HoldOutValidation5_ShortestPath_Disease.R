rm(list=ls())
library(igraph)

Positive_Negative_DrugsDiseases = readRDS("data/Positive_Negative_DrugsDiseases_Holdout_Analysis.rds")

ADRDP_Proteins_Holdout = readRDS("data/ADRDP_Proteins_Holdout.rds")
ADRs = names(ADRDP_Proteins_Holdout)

disease_gene = readRDS("data/disease_gene_protein.rds")
disease_gene$entrez_id_disaese = as.character(disease_gene$entrez_id_disaese)

PPI = read.csv("data/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

ShortestPath = list()

for(i in ADRs){
  
  Diseases_val = Data[[i]][["Diseases_val"]]
  shortest_path_for_all_diseases = c()
  for(j in Diseases_val){
    targets_Diseases_val_j = unique(disease_gene[disease_gene$mondo_id %in% j, 2])
    ADRDP_Proteins_Holdout_ADR_i = ADRDP_Proteins_Holdout[[i]]
    
    shortest_path = c()
    for(l in targets_Diseases_val_j){
      length_sp = sapply(ADRDP_Proteins_Holdout_ADR_i, function(k) {
        path = shortest_paths(PPI_Graph, k, to = l, output = "vpath")
        path = length(path$vpath[[1]])
        path = ifelse(path>0, path-1, 7)
        return(path)
      })
      shortest_path = c(shortest_path, min(length_sp))
    }
    shortest_path_for_all_diseases = c(shortest_path_for_all_diseases, min(shortest_path))
  }
  
  ##############
  #Negative Samples
  Diseases_negativeSamples = Data[[i]][["Diseases_negativeSamples"]]
  shortest_path_for_NegSample_all_Diseases_all_repeats = c()
  for(r in 1:nrow(Diseases_negativeSamples)){
    Diseases_negativeSamples_r = Diseases_negativeSamples[r,]
    
    shortest_path_for_NegSample_all_Diseases = c()
    for(j2 in Diseases_negativeSamples_r){
      targets_Diseases_Negsample_j = unique(disease_gene[disease_gene$mondo_id %in% j2, 2])
      ADRDP_Proteins_Holdout_ADR_i = ADRDP_Proteins_Holdout[[i]]
      
      v = V(PPI_Graph)
      targets_Diseases_Negsample_j = intersect(names(v), targets_Diseases_Negsample_j)
      
      
      shortest_path_NS = c()
      for(l in targets_Diseases_Negsample_j){
        length_sp = sapply(ADRDP_Proteins_Holdout_ADR_i, function(k) {
          path = shortest_paths(PPI_Graph, k, to = l, output = "vpath")
          path = length(path$vpath[[1]])
          path = ifelse(path>0, path-1, 7)
          return(path)
        })
        shortest_path_NS = c(shortest_path_NS, min(length_sp))
      }
      shortest_path_for_NegSample_all_Diseases = c(shortest_path_for_NegSample_all_Diseases, min(shortest_path_NS))
    }
    shortest_path_for_NegSample_all_Diseases_all_repeats = c(shortest_path_for_NegSample_all_Diseases_all_repeats, shortest_path_for_NegSample_all_Diseases)
  }
  ShortestPath[[i]][["pos"]] = shortest_path_for_all_diseases
  ShortestPath[[i]][["neg"]] = shortest_path_for_NegSample_all_Diseases_all_repeats
  
}


saveRDS(ShortestPath, "data/Shortest_path_diseases.rds")

