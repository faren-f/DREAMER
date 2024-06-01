rm(list = ls())
library(openxlsx)
library(dplyr)

ADRDP_Proteins = readRDS("Data/ADRDP_Proteins.rds")
meddra_Sider = read.csv("Data/meddra.tsv",sep = '\t', header = FALSE)
conversion_table = read.csv("Data/Conversion_Proteins.csv")

meddra_Sider = meddra_Sider[,3:4]
meddra_Sider$V3 = paste0("meddra.", meddra_Sider$V3)
colnames(meddra_Sider) = c("meddra_id", "ADR_name")
  
adr_df <- lapply(names(ADRDP_Proteins), function(adr) {
  data.frame(meddra_id = adr, entrez = ADRDP_Proteins[[adr]], stringsAsFactors = FALSE)
}) %>% bind_rows()

adr_df = merge(adr_df, conversion_table, by = "entrez")
adr_df = merge(adr_df, meddra_Sider, by = "meddra_id")

adr_df = adr_df[,c(1,5,2,3)]
adr_df = unique(adr_df)
write.table(adr_df, "Data/TableS3_ADRDP_Proteins_DREAMER_STRING.csv", 
            sep = ",", row.names = FALSE)
write.xlsx(adr_df, file = "Data/TableS3_ADRDP_Proteins_DREAMER_STRING.xlsx")






