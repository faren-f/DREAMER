rm(list=ls())

library(ggplot2)
setwd("~/Desktop/ADR_Mechanism_Prediction/R/DREAMER3-STRING_DrugBank/")
ShortestPath = readRDS("Data/ProcessedData/Final_Result/Validation/Shortest_path_diseases.rds")
ShortestPath = readRDS("Data/ProcessedData/Final_Result/Validation/Shortest_path_drugs.rds")
###################
p1 = c()
p2 = c()
for(i in names(ShortestPath)){
  p1 = c(p1, ShortestPath[[i]][["pos"]])
  p2 = c(p2, ShortestPath[[i]][["neg"]])
}

min(p1); max(p1)
min(p2); max(p2)

wilcox.test(p1, p2, alternative = "less")
##############################################
Pvals = c()
Odds = c()
for(i in 0:(max(max(p1), max(p2))-1)){
  p_pos = p1; p_neg = p2
  p_pos = p_pos>i
  p_neg= p_neg>i
  Table = matrix(c(sum(!p_pos), sum(p_pos), sum(!p_neg), sum(p_neg)), nrow = 2, byrow = TRUE)
  F1 = fisher.test(Table, alternative = "greater")
  Pvals = c(Pvals, F1$p.value)
  Odds = c(Odds, F1$estimate[["odds ratio"]])
}

qval = p.adjust(Pvals, method = "BH")

# saveRDS(Pvals, "Data/ProcessedData/Final_Result/Validation/Pvals_atEachThreshol_Diseases_.rds")
# saveRDS(Odds, "Data/ProcessedData/Final_Result/Validation/Odds_atEachThreshold_Diseases.rds")

# saveRDS(Pvals, "Data/ProcessedData/Final_Result/Validation/Pvals_atEachThreshol_Drugs_.rds")
# saveRDS(Odds, "Data/ProcessedData/Final_Result/Validation/Odds_atEachThreshold_Drugs.rds")

barplot(Pvals)
#############Plot
data = data.frame(value = c(p1, p2), 
                  group = factor(c(rep("pos", length(p1)), rep("neg", length(p2)))))

ggplot(data, aes(x = group, y = value)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin plot of Positive and Negative Values",
       x = "Group",
       y = "Values")

boxplot(p1,p2)
wilcox.test(p1, p2, alternative = "less")




