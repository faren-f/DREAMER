rm(list=ls())

library(igraph)

setwd("~/Desktop/ADR_Mechanism_Prediction/R/DREAMER3-STRING_DrugBank/")

#Data = readRDS("Data/Validation/Data_DREAMER2_withoutScores.rds")
Data = readRDS("Data/ProcessedData/Final_Result/Validation/Data_DREAMER3_withoutScores.rds")

#adj_Proteins = readRDS("Data/Validation/adjustedProteins_after_correction.rds")
adj_Proteins = readRDS("Data/ProcessedData/Final_Result/Validation/adjustedProteins_after_correction.rds")
se = names(adj_Proteins)

#disease_gene = readRDS("Data/disease_gene_protein.rds")
disease_gene = readRDS("Data/ProcessedData/Final_Network/disease_gene_protein.rds")
disease_gene$entrez_id_disaese = as.character(disease_gene$entrez_id_disaese)

#PPI = read.csv("Data/PPI_entrez.csv")
PPI = read.csv("Data/ProcessedData/Final_Network/PPI_entrez_STRING800.csv")

PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)


# clusterExport(cl, c("Data", "adj_Proteins",
#                     "disease_gene", "PPI_Graph"))
# clusterEvalQ(cl,c(library(igraph)))

rep = 0
ShortestPath = list()

for(i in se){
  print(rep); rep = rep+1
#ParLoop = function(i){
  Diseases_test = Data[[i]][["Diseases_test"]]
  shortest_path_for_all_diseases = c()
  for(j in Diseases_test){
    targets_Diseases_test_j = unique(disease_gene[disease_gene$mondo_id %in% j, 2])
    identified_genes = adj_Proteins[[i]]
    
    shortest_path = c()
    for(l in targets_Diseases_test_j){
      length_sp = sapply(identified_genes, function(k) {
        path = shortest_paths(PPI_Graph, k, to = l, output = "vpath")
        path = length(path$vpath[[1]])
        path = ifelse(path>0, path-1, 7)
        return(path)
      })
      shortest_path = c(shortest_path, min(length_sp))
    }
    #shortest_path_for_all_diseases = c(shortest_path_for_all_diseases, sum(shortest_path<=1)/length(shortest_path))
    shortest_path_for_all_diseases = c(shortest_path_for_all_diseases, min(shortest_path))
  }
  
  ##############
  #Negative Samples
  Diseases_negativeSamples = Data[[i]][["Diseases_negativeSamples"]]
  shortest_path_for_NegSample_all_Diseases_all_repeats = c()
  #for(r in 1:nrow(Diseases_negativeSamples)){
  for(r in 1:1){
    Diseases_negativeSamples_r = Diseases_negativeSamples[r,]
    
    shortest_path_for_NegSample_all_Diseases = c()
    for(j2 in Diseases_negativeSamples_r){
      targets_Diseases_Negsample_j = unique(disease_gene[disease_gene$mondo_id %in% j2, 2])
      identified_genes = adj_Proteins[[i]]
      
      v = V(PPI_Graph)
      targets_Diseases_Negsample_j = intersect(names(v), targets_Diseases_Negsample_j)
      
      
      shortest_path_NS = c()
      for(l in targets_Diseases_Negsample_j){
        length_sp = sapply(identified_genes, function(k) {
          path = shortest_paths(PPI_Graph, k, to = l, output = "vpath")
          path = length(path$vpath[[1]])
          path = ifelse(path>0, path-1, 7)
          return(path)
        })
        shortest_path_NS = c(shortest_path_NS, min(length_sp))
      }
      #shortest_path_for_NegSample_all_Diseases = c(shortest_path_for_NegSample_all_Diseases, sum(shortest_path_NS<=1)/length(shortest_path_NS))
      shortest_path_for_NegSample_all_Diseases = c(shortest_path_for_NegSample_all_Diseases, min(shortest_path_NS))
    }
    shortest_path_for_NegSample_all_Diseases_all_repeats = c(shortest_path_for_NegSample_all_Diseases_all_repeats, shortest_path_for_NegSample_all_Diseases)
  }
  #ShortestPath = list()
  ShortestPath[[i]][["pos"]] = shortest_path_for_all_diseases
  ShortestPath[[i]][["neg"]] = shortest_path_for_NegSample_all_Diseases_all_repeats
  
  
  #return(shortest_path)
  
}

#ShortestPath = parLapply(cl, sapply(se, list), ParLoop) 

#saveRDS(ShortestPath, "Data/Validation/Shortest_path.rds")
saveRDS(ShortestPath, "Data/ProcessedData/Final_Result/Validation/Shortest_path_diseases.rds")

#stopCluster(cl)
