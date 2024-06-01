rm(list = ls())
library(jsonlite)
setwd("~/Desktop/ADR_Mechanism_Prediction/R/DREAMER3-STRING_DrugBank/")

TrainTestSplits = fromJSON("Data/ProcessedData/Final_Network/TrainTestSplits_forEachADR.json")

saveRDS(TrainTestSplits, "Data/ProcessedData/Final_Network/TrainTestSplits_forEachADR.rds")

