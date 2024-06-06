# In this script, I want to depict the dot plots of pathways enriched by Reactome and GO

rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)

library(biomaRt)
library(ReactomePA)
mart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

Final_ADRs = readRDS("Data/Final_67_ADRs_STRING_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")
ADRDP_Proteins = readRDS("Data/ADRDP_Proteins.rds")
drug_ADR = readRDS("Data/drug_ADR.rds")
disease_symptom = readRDS("Data/disease_symptom.rds")


Final_ADRs = unique(drug_ADR[drug_ADR$se_name %in% Final_ADRs,2])
ADRDP_Proteins = ADRDP_Proteins[Final_ADRs]

i = "meddra.10020957" # The meddra_id of the ADR for which we want to depict the dot plot for example hypochloremia

gene_ids = ADRDP_Proteins[[i]]
if(length(gene_ids) == 0){
  next
}
ADR_name = unique(drug_ADR[drug_ADR$meddra_id %in% i,3])

#########Reactome
reactom_pw = enrichPathway(gene = gene_ids, organism = "human", pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", qvalueCutoff = 0.01, minGSSize = 10,
                           maxGSSize = 500, readable = FALSE)

pdf(paste0("Figures/", ADR_name, "_Reactome.pdf"), width = 5, height = 4)
dotplot(reactom_pw, showCategory = 5)+
  theme(axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size=.2, color="lightgray"),
        panel.grid.minor.y = element_line(size=.2, color="lightgray"))+
  ggtitle(ADR_name)
dev.off()


#########GO

GO_pw = enrichGO(gene = gene_ids,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENTREZID",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 readable = TRUE)
pdf(paste0("Figures/", ADR_name, "_GO.pdf"), width = 5, height = 5)
dotplot(GO_pw, showCategory = 5)+ 
  theme(axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(size=.2, color="lightgray"),
        panel.grid.minor.y = element_line(size=.2, color="lightgray"))+
  ggtitle(ADR_name)

dev.off()


