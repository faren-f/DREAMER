rm(list = ls())
library(writexl)
library(igraph)

setwd("~/Desktop/DrugSiderAI/MyPaper/Codes/DREAMER/")
Conversion_Proteins = read.csv("data/knowledge_graph/Conversion_Proteins.csv")
meddra_se_sider = read.csv("data/knowledge_graph/meddra.tsv",sep = '\t', header = FALSE)
PPI = read.csv("data/knowledge_graph/PPI_STRING.csv")
PPI_Graph = graph_from_data_frame(PPI, directed = FALSE)
PPI_Graph = simplify(PPI_Graph, remove.loops = TRUE, remove.multiple = TRUE)
PPI_gene_ids = V(PPI_Graph)$name
PPI_geneNames_ids = Conversion_Proteins[Conversion_Proteins$entrez %in% PPI_gene_ids,]

###########################################
# Reading PubMed data
symptom_gene = read.csv("data/other_resources/literature_symptom_gene_associations.tsv", sep = "\t", header = TRUE)
symptom_gene = symptom_gene[,c(1,3)]

I = intersect(symptom_gene$Gene_Symbol, PPI_geneNames_ids[,1])
symptom_gene = symptom_gene[symptom_gene$Gene_Symbol%in% I,]

###########################################
# Read ADR-related proteins and DP-related proteins
ADRProtein_DPProtein = readRDS("data/ADRRelatedProteins_DPRelatedProteins.rds")

###########################################Baseline
# significant overlap association (SOA) Baseline method
SOA_ADRProteins = readRDS("data/ADR_Proteins_baseline.rds")
SOA_DPProteins = readRDS("data/DP_Proteins_baseline.rds")

#How many ADRs are common between PubMed literature-based and our network
all_SEs = names(ADRProtein_DPProtein)

UMLS_meddra = unique(meddra_se_sider[meddra_se_sider$V3 %in% gsub("meddra.","", all_SEs),c(1,3)])
colnames(UMLS_meddra) = c("umls_id", "meddra_id")
UMLS_meddra$meddra_id = paste0("meddra.", UMLS_meddra$meddra_id)

I = intersect(unique(symptom_gene$Symptom_CUI), UMLS_meddra$umls_id)
symptom_gene = symptom_gene[symptom_gene$Symptom_CUI %in% I,]
UMLS_meddra = UMLS_meddra[UMLS_meddra$umls_id %in% I,]

symptom_gene = merge(symptom_gene, UMLS_meddra, 
                            by.x = colnames(symptom_gene)[1],
                            by.y = colnames(UMLS_meddra)[1])

symptom_gene = merge(symptom_gene, PPI_geneNames_ids, 
                            by.x = colnames(symptom_gene)[2],
                            by.y = colnames(PPI_geneNames_ids)[1])

################################################################################
####Overlap SOA_DP-Proteins with PubMed
Phenotypes_SOA_DPProtein = c()

for(i in UMLS_meddra[,2]){

  G1_SOA = unique(SOA_DPProteins[SOA_DPProteins$SE_name %in% i, 1])
  G2 = unique(as.character(symptom_gene[symptom_gene$meddra_id %in% i,4]))
  
  I = intersect(G1_SOA, G2)
  Phenotypes_SOA_DPProtein = rbind(Phenotypes_SOA_DPProtein, c(length(G1_SOA),length(G2), length(I)))
}

rownames(Phenotypes_SOA_DPProtein) = UMLS_meddra[,2]
Phenotypes_SOA_DPProtein = data.frame(Phenotypes_SOA_DPProtein)
colnames(Phenotypes_SOA_DPProtein) = c("N_Pr_SOA", "N_Pr_PubMed", "N_Overlap")

Phenotypes_SOA_DPProtein$pval = rep(1, nrow(Phenotypes_SOA_DPProtein))
for(i in 1:nrow(Phenotypes_SOA_DPProtein)){
  Nt = length(PPI_gene_ids)
  N1 = Phenotypes_SOA_DPProtein[i,1]
  N2 = Phenotypes_SOA_DPProtein[i,2]
  K = Phenotypes_SOA_DPProtein[i,3]
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  Phenotypes_SOA_DPProtein$pval[i] = fisher.test(mat, alternative = "greater")$p.value
}

Phenotypes_SOA_DPProtein$qval = p.adjust(Phenotypes_SOA_DPProtein$pval, "BH")
Phenotypes_SOA_DPProtein = Phenotypes_SOA_DPProtein[Phenotypes_SOA_DPProtein$qval<0.05,]

################################################################################
####Overlap SOA_ADR-Proteins with PubMed

Phenotypes_SOA_ADRProtein = c()
for(i in UMLS_meddra[,2]){
  G1_SOA = unique(SOA_ADRProteins[SOA_ADRProteins$SE_name %in% i, 1])
  G2 = unique(as.character(symptom_gene[symptom_gene$meddra_id %in% i,4]))
  
  I = intersect(G1_SOA, G2)
  Phenotypes_SOA_ADRProtein = rbind(Phenotypes_SOA_ADRProtein, c(length(G1_SOA),length(G2), length(I)))
}

rownames(Phenotypes_SOA_ADRProtein) = UMLS_meddra[,2]
Phenotypes_SOA_ADRProtein = data.frame(Phenotypes_SOA_ADRProtein)
colnames(Phenotypes_SOA_ADRProtein) = c("N_Pr_SOA", "N_Pr_PubMed", "N_Overlap")

Phenotypes_SOA_ADRProtein$pval = rep(1, nrow(Phenotypes_SOA_ADRProtein))
for(i in 1:nrow(Phenotypes_SOA_ADRProtein)){
  Nt = length(PPI_gene_ids)
  N1 = Phenotypes_SOA_ADRProtein[i,1]
  N2 = Phenotypes_SOA_ADRProtein[i,2]
  K = Phenotypes_SOA_ADRProtein[i,3]
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  Phenotypes_SOA_ADRProtein$pval[i] = fisher.test(mat, alternative = "greater")$p.value
}
Phenotypes_SOA_ADRProtein$qval = p.adjust(Phenotypes_SOA_ADRProtein$pval, "BH")
Phenotypes_SOA_ADRProtein = Phenotypes_SOA_ADRProtein[Phenotypes_SOA_ADRProtein$qval<0.05,]

######################################
####Overlap Diffusion DP-Protein with PubMed

Phenotypes_diffusion_DPProtein = c()
r = 0
for(i in UMLS_meddra[,2]){
  print(r); r = r+1
  G1 = unique(as.character(ADRProtein_DPProtein[[i]][["Disease"]]))
  G2 = unique(as.character(symptom_gene[symptom_gene$meddra_id %in% i,4]))
  I = intersect(G1, G2)
  Phenotypes_diffusion_DPProtein = rbind(Phenotypes_diffusion_DPProtein, c(length(G1),length(G2), length(I)))
  
}
rownames(Phenotypes_diffusion_DPProtein) = UMLS_meddra[,2]
Phenotypes_diffusion_DPProtein = data.frame(Phenotypes_diffusion_DPProtein)
colnames(Phenotypes_diffusion_DPProtein) = c("N_Pr_DiffDisease", "N_Pr_PubMed", "N_Overlap")

Phenotypes_diffusion_DPProtein$pval = rep(1, nrow(Phenotypes_diffusion_DPProtein))
for(i in 1:nrow(Phenotypes_diffusion_DPProtein)){
  Nt = length(PPI_gene_ids)
  N1 = Phenotypes_diffusion_DPProtein[i,1]
  N2 = Phenotypes_diffusion_DPProtein[i,2]
  K = Phenotypes_diffusion_DPProtein[i,3]
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  Phenotypes_diffusion_DPProtein$pval[i] = fisher.test(mat, alternative = "greater")$p.value
}

Phenotypes_diffusion_DPProtein$qval = p.adjust(Phenotypes_diffusion_DPProtein$pval, "BH")
Phenotypes_diffusion_DPProtein = Phenotypes_diffusion_DPProtein[Phenotypes_diffusion_DPProtein$qval<0.05,]

######################################
####Overlap Diffusion ADR-Protein with PubMed

Phenotypes_diffusion_ADRProtein = c()
r = 0
for(i in UMLS_meddra[,2]){
  print(r); r = r+1
  
  G1 = unique(as.character(ADRProtein_DPProtein[[i]][["Drug"]]))
  G2 = unique(as.character(symptom_gene[symptom_gene$meddra_id %in% i,4]))
  I = intersect(G1, G2)
  Phenotypes_diffusion_ADRProtein = rbind(Phenotypes_diffusion_ADRProtein, c(length(G1),length(G2), length(I)))
}
rownames(Phenotypes_diffusion_ADRProtein) = UMLS_meddra[,2]
Phenotypes_diffusion_ADRProtein = data.frame(Phenotypes_diffusion_ADRProtein)
colnames(Phenotypes_diffusion_ADRProtein) = c("N_Pr_DiffDrug", "N_Pr_PubMed", "N_Overlap")

Phenotypes_diffusion_ADRProtein$pval = rep(1, nrow(Phenotypes_diffusion_ADRProtein))
for(i in 1:nrow(Phenotypes_diffusion_ADRProtein)){
  Nt = length(PPI_gene_ids)
  N1 = Phenotypes_diffusion_ADRProtein[i,1]
  N2 = Phenotypes_diffusion_ADRProtein[i,2]
  K = Phenotypes_diffusion_ADRProtein[i,3]
  mat = matrix(c(K, N1-K, N2-K, Nt-(N1+N2-K)), nrow = 2, byrow = TRUE)
  Phenotypes_diffusion_ADRProtein$pval[i] = fisher.test(mat, alternative = "greater")$p.value
}

Phenotypes_diffusion_ADRProtein$qval = p.adjust(Phenotypes_diffusion_ADRProtein$pval, "BH")
Phenotypes_diffusion_ADRProtein = Phenotypes_diffusion_ADRProtein[Phenotypes_diffusion_ADRProtein$qval<0.05,]

############## BarPlot
df1 = data.frame(methods = "ADR-proteins baseline", values = nrow(Phenotypes_SOA_ADRProtein))
df2 = data.frame(methods = "DP-proteins baseline", values = nrow(Phenotypes_SOA_DPProtein))
df3 = data.frame(methods = "ADR-proteins diffusion", values = nrow(Phenotypes_diffusion_ADRProtein))
df4 = data.frame(methods = "DP-proteins diffusion", values = nrow(Phenotypes_diffusion_DPProtein))

data = rbind(df1, df2, df3, df4)
data$methods = factor(data$methods, levels = c("ADR-proteins baseline", "DP-proteins baseline", "ADR-proteins diffusion", "DP-proteins diffusion"))

ggplot(data, aes(x=methods, y=values, fill=methods)) + 
  geom_col(position="dodge", width=0.5, color="black") +
  scale_fill_manual(values=c("#B71C1C","#B71C1C", "darkblue", "darkblue")) + 
  theme_classic() +scale_y_continuous(n.breaks = 10)+
  scale_x_discrete(guide = guide_axis(angle = 60))






