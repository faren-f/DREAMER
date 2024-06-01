rm(list = ls())
library(openxlsx)

drug_ADR = readRDS("Data/drug_ADR.rds")

Final_ADRs_Physical = readRDS("Data/Final_56_ADRs_Physical_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")
rownames(Final_ADRs_Physical) = Final_ADRs_Physical[,1]

Final_ADRs_STRING = readRDS("Data/Final_67_ADRs_STRING_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")
rownames(Final_ADRs_STRING) = Final_ADRs_STRING[,1]

intersected_ADRs = intersect(Final_ADRs_Physical[,1], Final_ADRs_STRING[,1])
intersected_ADRs_names = unique(drug_ADR[drug_ADR$meddra_id %in% intersected_ADRs, 3])


saveRDS(intersected_ADRs_names,"Data/Intersected_ADRs_Physical_STRING.rds")
saveRDS(intersected_ADRs,"Data/Intersected_ADRs_Physical_STRING.rds")

#########
# Phenotype identified in STRING but not Physical
ADR_Only_STRING = rownames(Final_ADRs_STRING)[!(Final_ADRs_STRING[,1] %in% intersected_ADRs)]
ADR_Only_Physical = rownames(Final_ADRs_Physical)[!(Final_ADRs_Physical[,1] %in% intersected_ADRs)]










