rm(list = ls())

library(ggplot2)
library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

drug_ADR = readRDS("Data/drug_ADR.rds")
drug_target = readRDS("Data/drug_target.rds")

targets = as.character(unique(drug_target$entrez_id_drug))
ADRs = unique(drug_ADR$meddra_id)
# number of all the drugs
Nt = length(unique(drug_target$drugbank_id))

clusterExport(cl, c("drug_target", "drug_ADR", "targets", "ADRs","Nt"))
clusterEvalQ(cl,c())

ParLoop = function(i){
  
  #drugs connected to protein_i
  dt = drug_target[drug_target$entrez_id_drug == i, 1]
  N1 = length(unique(dt))
  result = c()
  for(j in ADRs){
    #drugs connected to ADR_j 
    ds = drug_ADR[drug_ADR$meddra_id == j, 1]
    N2 = length(unique(ds))
    K = length(intersect(dt, ds))
    mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
    Fisher_Test_1SG = fisher.test(mat, alternative = "greater")$p.value
    
    result = rbind(result, data.frame(Target_name = i, SE_name = j, 
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

ADR_Proteins_baseline = Result[Result$qval<0.05,]

saveRDS(ADR_Proteins_baseline, "Data/ADR_Proteins_baseline.rds")


