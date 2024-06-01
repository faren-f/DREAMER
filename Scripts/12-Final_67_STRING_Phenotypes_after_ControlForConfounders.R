rm(list = ls())
library(ggplot2)
library(dplyr)

drug_ADR = readRDS("Data/drug_ADR.rds")
conversion_table = read.csv("Data/Conversion_Proteins.csv")

ADRDP_Proteins = readRDS("Data/ADRDP_Proteins.rds")
ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan = readRDS("Data/ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan.rds")
Indication_Proteins = readRDS("Data/Indication_Proteins.rds")

sig_ADRs = intersect(ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan, names(ADRDP_Proteins))
sig_ADRs = intersect(names(Indication_Proteins), sig_ADRs)

# ADRs that have significant overlap with indication and should be removed
ADRs_withSignificantOverlapWithIndicationProteins = readRDS("Data/ADRs_withSignificantOvellap_ADRDPandIndication_Proteins.rds")

sig_ADRs = sig_ADRs[!sig_ADRs %in% ADRs_withSignificantOverlapWithIndicationProteins]

Final_ADRs = unique(drug_ADR[drug_ADR$meddra_id %in% sig_ADRs, 2:3])

# save the the name and ids of ADRs for which we obtained mechanism afetr controlling for confounders leveraging STRING PPI:
saveRDS(Final_ADRs, "Data/Final_67_ADRs_STRING_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")

# save the the name and ids of ADRs for which we obtained mechanism afetr controlling for confounders leveraging physical PPI:
#saveRDS(Final_ADRs, "Data/Final_56_ADRs_Physical_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")





