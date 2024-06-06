rm(list = ls())
library(openxlsx)

drug_ADR = readRDS("data/preprocessed_graph/drug_ADR.rds")

ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan = readRDS("data/ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan.rds")
ADRDP_Proteins = readRDS("data/ADRDP_Proteins.rds")

ADRDP_Proteins_AfterDrugRemoval = readRDS("data/ADRDP_Proteins_AfterDrugRemoval.rds")
ADRDP_Proteins_AfterDrugRemoval = ADRDP_Proteins_AfterDrugRemoval[ADRs_withSignificantOverlap_afterDrugRemovalSameOrgan]

ADRs_afterOrganRemoval = readRDS("data/ADRs_afterOrganRemoval.rds")
ADR_associated = ADRs_afterOrganRemoval[!(ADRs_afterOrganRemoval %in% intersect(names(ADRDP_Proteins_AfterDrugRemoval), ADRs_afterOrganRemoval))]
No_overlapping_drugs_after_removal = names(ADRDP_Proteins)[!(names(ADRDP_Proteins) %in% ADRs_afterOrganRemoval)]

ADR_table = data.frame(ADR_name = names(ADRDP_Proteins))
ADR_table[ADR_table$ADR_name %in% names(ADRDP_Proteins_AfterDrugRemoval),2] = "No association with organ"
ADR_table[ADR_table$ADR_name %in% ADR_associated ,2] = "associated with organ"
ADR_table[ADR_table$ADR_name %in% No_overlapping_drugs_after_removal ,2] = "ADRs without drugs after removal"
colnames(ADR_table)[2] = "Organ_Confounding_Control"


#########Indication-diffused Control
indication_proteins = readRDS("data/Indication_Proteins.rds")
significant_ADRName_Indication = readRDS("data/significant_ADRName_Indication.rds")
significant_ADRName_Indication = unique(drug_se[drug_se$se_name %in% significant_ADRName_Indication, 2])


# ADRs without at least one drug with at least one indication that has at least one associated gene listed in our dataset
ADR_table[!(ADR_table$ADR_name %in% names(indication_proteins)),3] = "ADRs without at least one indication that has at least one associated gene listed in our dataset"
ADR_table[ADR_table$ADR_name %in% significant_ADRName_Indication,3] = "associted with indications"

Non_associted_to_indications = names(indication_proteins)[!names(indication_proteins) %in% significant_ADRName_Indication]
ADR_table[ADR_table$ADR_name %in% Non_associted_to_indications,3] = "No association with indications-related proteins"
colnames(ADR_table)[3] = "Confounding_control_after_indication_diffused"

write.table(ADR_table, "data/TableS7_ADRs_associatedWith_Confounders_STRING.csv", 
            sep = ",", row.names = FALSE)
write.xlsx(ADR_table, file = "data/TableS7_ADRs_associatedWith_Confounders_STRING.xlsx")



