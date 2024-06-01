rm(list = ls())

library(ggplot2)
library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

disease_symptom = readRDS("Data/disease_symptom.rds")
disease_gene_protein = readRDS("Data/disease_gene_protein.rds")

targets = as.character(unique(disease_gene_protein$entrez_id_disaese))
symptoms = unique(disease_symptom$meddra_id)

# number of all the diseases
Nt = length(unique(disease_gene_protein$mondo_id))

clusterExport(cl, c("disease_gene_protein", "disease_symptom", "targets", "symptoms","Nt"))
clusterEvalQ(cl,c())

ParLoop = function(i){
  
  #diseases connected to protein_i
  dt = disease_gene_protein[disease_gene_protein$entrez_id_disaese == i, 1]
  N1 = length(unique(dt))
  result = c()
  for(j in symptoms){
    ds = disease_symptom[disease_symptom$meddra_id == j, 1]
    N2 = length(unique(ds))
    K = length(intersect(dt, ds))
    mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
    Fisher_Test_1SG = fisher.test(mat, alternative = "greater")$p.value
    result = rbind(result, data.frame(entrez_id_target = i, SE_name = j, 
                                      N1 = N1, N2 = N2, Overlap = K, 
                                      Fisher_Test_1SG = Fisher_Test_1SG))
    
  }
  return(result)
}
result = parLapply(cl, sapply(targets, list), ParLoop) 

Result = data.frame()
for (k in targets){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)

qval = p.adjust(Result$Fisher_Test_1SG, "BH")
Result$qval = qval

DP_Proteins_baseline = Result[Result$qval<0.05,]
saveRDS(DP_Proteins_baseline, "Data/DP_Proteins_baseline.rds")


