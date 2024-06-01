rm(list = ls())
library(igraph)
library(openxlsx)

#Read data
drug_ADR = readRDS("Data/drug_ADR.rds")
drug_target = readRDS("Data/drug_target.rds")
disease_symptom = readRDS("Data/disease_symptom.rds")
drug = read.csv("Data/DrugNodes.csv")
mondo_meddra_list = readRDS("Data/mondo_meddra_list.rds")
mondo_meddra_df = readRDS("Data/mondo_meddra_df.rds")

drug_indication = read.table("Data/direct_drug_indications.csv", sep = "-")[,2:1]
colnames(drug_indication) = c("drugbank_id", "mondo_id")
Conversion_Proteins = read.csv("Data/Conversion_Proteins_new.csv")

ADRDP_Proteins = readRDS("Data/ADRDP_Proteins.rds")
Final_ADRs = readRDS("Data/Final_67_ADRs_STRING_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")
ADRDP_Proteins = ADRDP_Proteins[Final_ADRs[,1]]

drug_name = unique(drug_ADR[,c(1,4)])

##########
# Finding the list of drugs that their indication is the same as the given ADR 

meddra_drug = list()
for(i in names(ADRDP_Proteins)){
  ADRDP_Proteins_ADR_i = ADRDP_Proteins[[i]]
  
  # all of the drugs that have target in the ADRDP_Proteins_ADR_i
  Dr_adj = unique(drug_target[drug_target$entrez_id_drug %in% ADRDP_Proteins_ADR_i, 1])
  
  # drugs that have the given ADR and have target in ADRDP_Proteins_ADR_i
  Dr_ADR = drug_ADR[drug_ADR$meddra_id %in% i, 1]
  Dr_ADR = intersect(Dr_ADR, Dr_adj)
  
  #Excluding all of the drugs that have ADR and target in ADRDP_Proteins_ADR_i from Dr_adj
  D_other = unique(Dr_adj[!(Dr_adj %in% intersect(Dr_adj, Dr_ADR))])
  
  #To find Dr_indication we should take some steps:
  #1) from D_other excluding drugs that do not have indication in the drug_indication and for the rest finding their indication
  D_other_indication = unique(drug_indication[drug_indication$drugbank_id %in% D_other,])
  
  ######
  length(D_other)
  length(unique(D_other_indication$drugbank_id))
  ######
  
  
  #2) from D_other_indication excluding disaeses that do not have equivalent meddra_id in the mondo_meddra_list
  if(nrow(D_other_indication)>0){
    Mondo_drug = D_other_indication[D_other_indication$mondo_id %in% names(mondo_meddra_list),]
    
    #######
    length(unique(D_other_indication$mondo_id))
    length(unique(Mondo_drug$mondo_id))
    #######
    
    # making a dataframe from list of indication_meddra_ids
    Mondo = names(mondo_meddra_list)[names(mondo_meddra_list) %in% D_other_indication$mondo_id]
    Meddra_Mondo_df_list = lapply(Mondo, function(x) {
      data.frame(
        mondo_id = x,
        meddra_id = mondo_meddra_list[[x]]
      )
    })
    
    Meddra_Mondo = do.call(rbind, Meddra_Mondo_df_list)
    Meddra_Mondo = unique(Meddra_Mondo)
    
    # constructing a dataframe that includes all the meddra_mondo_drug that 
    # drugs are among D_other, diseases are among the indications of that drugs and ADRs 
    # have equivalent in mondo
    meddra_mondo_drug = merge(Mondo_drug, Meddra_Mondo, by = "mondo_id")
    
    if(length(intersect(Meddra_Mondo$meddra_id, i))!= 0){
      meddra_drug[[i]] = meddra_mondo_drug[meddra_mondo_drug$meddra_id %in% i,2]
    }
  }
} 

#########################################
# meddra_drug: list of the drugs that have indication the same as ADR

indicative_drug_protein = c()
candidate_drug_protein = c()

for(i in names(meddra_drug)){
  
  ADR = unique(drug_ADR[drug_ADR$meddra_id %in% i, 3])
  pr_i_AND = ADRDP_Proteins[[i]]
  
  ## Drugs1:
  drugs_ADR_i = drug_ADR[drug_ADR$meddra_id %in% i,1]
  drugs_ADR_i_target_i = unique(drug_target[drug_target$drugbank_id %in% drugs_ADR_i & 
                                             drug_target$entrez_id_drug %in% 
                                             pr_i_AND,1:2])
  ## Drugs2:
  drugs_OtherwithIndicationSameAsADR = meddra_drug[[i]]
  drugs_OtherwithIndicationSameAsADR_target_i = unique(drug_target[drug_target$drugbank_id %in% 
                                                                    drugs_OtherwithIndicationSameAsADR & 
                                                                    drug_target$entrez_id_drug %in% 
                                                                    pr_i_AND,1:2])
  
  ## Proteins drug_ADR and Proteins Drug_indication
  P1 = unique(drugs_ADR_i_target_i$entrez_id_drug)
  P2 = unique(drugs_OtherwithIndicationSameAsADR_target_i$entrez_id_drug)
  List_drugs1 = unique(drugs_ADR_i_target_i$drugbank_id)
  List_drugs2 = unique(drugs_OtherwithIndicationSameAsADR_target_i$drugbank_id)
  drug_protein1 = merge(drugs_ADR_i_target_i, drug_name, "drugbank_id")
  drug_protein2 = merge(drugs_OtherwithIndicationSameAsADR_target_i, drug_name, "drugbank_id")
  
  # Yellow proteins that are targeted by some drug_ADR and Drug_indication
  K = intersect(drugs_ADR_i_target_i$entrez_id_drug, drugs_OtherwithIndicationSameAsADR_target_i$entrez_id_drug)
  
  # Proteins rest_drugs
  Dr_adj = unique(drug_target[drug_target$entrez_id_drug %in% pr_i_AND, 1])
  
  #Excluding all of the drugs that have ADR and target in ADRDP_Proteins_ADR_i from Dr_adj
  D_other = Dr_adj[!(Dr_adj %in% intersect(Dr_adj, drugs_ADR_i_target_i[,1]))]
  D_rest = D_other[!D_other %in% drugs_OtherwithIndicationSameAsADR_target_i[,1]]
  
  drugs_rest_i_target_i = unique(drug_target[drug_target$drugbank_id %in% D_rest & 
                                               drug_target$entrez_id_drug %in% pr_i_AND,1:2])
  
  drug_protein_rest = merge(drugs_rest_i_target_i, drug_name, "drugbank_id")
  List_drugs_rest = unique(drug_protein_rest$drugbank_id)

  List_drugs2_db_name = unique(drug_ADR[drug_ADR$drugbank_id %in% List_drugs2,c(1,4)])
  List_drugs_rest_db_name = unique(drug_ADR[drug_ADR$drugbank_id %in% List_drugs_rest,c(1,4)])
  
  indicative_drug_protein_i = c()
  candidate_drug_protein_i = c()
  
  indicative_drug_protein_i = merge(drug_protein2, Conversion_Proteins, by.x = "entrez_id_drug", by.y = "entrez")[,3:4]
  indicative_drug_protein_i = unique(indicative_drug_protein_i)
  
  candidate_drug_protein_i = merge(drug_protein_rest, Conversion_Proteins, by.x = "entrez_id_drug", by.y = "entrez")[,3:4]
  candidate_drug_protein_i = unique(candidate_drug_protein_i)
  
  indicative_drug_protein_i = cbind(indicative_drug_protein_i, rep(ADR, nrow(indicative_drug_protein_i)))
  colnames(indicative_drug_protein_i)[3] = "ADR"
  
  
  candidate_drug_protein_i = cbind(candidate_drug_protein_i, rep(ADR, nrow(candidate_drug_protein_i)))
  colnames(candidate_drug_protein_i)[3] = "ADR"
  
  indicative_drug_protein = rbind(indicative_drug_protein, indicative_drug_protein_i)
  candidate_drug_protein = rbind(candidate_drug_protein, candidate_drug_protein_i)
}

write.table(indicative_drug_protein, "Data/TableS10_indicative_drugs_targets.csv", sep = ",")
write.table(candidate_drug_protein, "Data/TableS11_candidate_drugs_targets.csv", sep = ",")

write.xlsx(indicative_drug_protein, file = "Data/TableS10_indicative_drugs_targets.xlsx")
write.xlsx(candidate_drug_protein, file = "Data/TableS11_candidate_drugs_targets.xlsx")


