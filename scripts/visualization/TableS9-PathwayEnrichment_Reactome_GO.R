rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(ReactomePA)
library(openxlsx)

mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

Final_ADRs = readRDS("data/Final_67_ADRs_STRING_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")
ADRDP_Proteins = readRDS("data/ADRDP_Proteins.rds")
drug_ADR = readRDS("data/preprocessed_graph/drug_ADR.rds")
disease_symptom = readRDS("data/preprocessed_graph/disease_symptom.rds")

ADRs = Final_ADRs[,1]
ADRDP_Proteins = ADRDP_Proteins[ADRs]

Pathways_all_Reactome = c()
Pathways_all_GO = c()
for(i in names(ADRDP_Proteins)){
  ADR_name = unique(drug_ADR[drug_ADR$meddra_id %in% i,3])
  gene_ids = ADRDP_Proteins[[i]]
  
  if(length(gene_ids) == 0){
    print(i)
    next
  }
  pw_Reactome = enrichPathway(gene = gene_ids, 
                              organism = "human", 
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH", 
                              qvalueCutoff = 0.01, 
                              minGSSize = 10,
                              maxGSSize = 500, 
                              readable = FALSE)
  
  Pathway_Reactome = pw_Reactome@result
  Pathway_Reactome = Pathway_Reactome[1:5,c(1,2,5,6,7)]
  Pathway_Reactome$ADR_meddra_id = rep(i,5)
  Pathway_Reactome$ADR_name = rep(ADR_name,5)
  Pathways_all_Reactome = rbind(Pathways_all_Reactome, Pathway_Reactome)
  ###########
  
  GO_pw = enrichGO(gene = gene_ids,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)
  Pathway_GO = GO_pw@result
  Pathway_GO = Pathway_GO[1:5,c(1,2,5,6,7)]
  Pathway_GO$ADR_meddra_id = rep(i,5)
  Pathway_GO$ADR_name = rep(ADR_name,5)
  Pathways_all_GO = rbind(Pathways_all_GO, Pathway_GO)
}

write.table(Pathways_all_Reactome, "data/TableS9_Reactome.csv", sep = ",", row.names = FALSE)
write.table(Pathways_all_GO, "data/TableS9_GO.csv", sep = ",", row.names = FALSE)

write.xlsx(Pathways_all_Reactome, file = "data/TableS9_Reactome.xlsx")
write.xlsx(Pathways_all_GO, file = "data/TableS9_GO.xlsx")


