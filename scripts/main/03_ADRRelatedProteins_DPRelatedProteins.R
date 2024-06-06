rm(list = ls())
library(igraph)

Pvalues = readRDS("data/Pval_ADRProteins_DPProteins.rds")

#Physical PPI
PPI = read.csv("data/knowledge_graph/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)

Pr_names = V(PPI_Graph)$name

adj_pval_ADRRelatedProteins_DPRelatedProteins = list()
ADRRelatedProteins_DPRelatedProteins = list()

for(i in names(Pvalues)){
  I_ADRProteins = which(Pvalues[[i]][["ADR_Protein"]]<0.05)
  I_DPProteins = which(Pvalues[[i]][["DP_Protein"]]<0.05)
  
  ADRRelatedProteins_DPRelatedProteins[[i]][["ADR_Protein"]] = Pr_names[I_ADRProteins]
  ADRRelatedProteins_DPRelatedProteins[[i]][["DP_Protein"]] = Pr_names[I_DPProteins]
  
  adj_pval_ADRRelatedProteins_DPRelatedProteins[[i]][["ADR_Protein"]] = Pvalues[[i]][["ADR_Protein"]][Pvalues[[i]][["ADR_Protein"]]<0.05]
  adj_pval_ADRRelatedProteins_DPRelatedProteins[[i]][["DP_Protein"]] = Pvalues[[i]][["DP_Protein"]][Pvalues[[i]][["DP_Protein"]]<0.05]
  
}
saveRDS(adj_pval_ADRRelatedProteins_DPRelatedProteins, "data/adj_pval_ADRRelatedProteins_DPRelatedProteins.rds")
saveRDS(ADRRelatedProteins_DPRelatedProteins, "data/ADRRelatedProteins_DPRelatedProteins.rds")




