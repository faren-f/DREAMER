rm(list = ls())
library("readxl")

ADRDP_Proteins = readRDS("data/ADRDP_Proteins.rds")
ADR_Tissue = readRDS("data/knowledge_graph/ADR_Tissue.rds")
length(intersect(rownames(ADR_Tissue), names(ADRDP_Proteins)))

write.table(ADR_Tissue, "data/TableS6_ADR_Tissue.csv", 
            sep = ",", row.names = FALSE)

write.xlsx(ADR_Tissue, "data/TableS6_ADR_Tissue.xlsx")
