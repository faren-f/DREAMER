rm(list=ls())
library(ggplot2)

ShortestPath_drug = readRDS("Data/Shortest_path_drugs.rds")
ShortestPath_disaese = readRDS("Data/Shortest_path_diseases.rds")
###################

p1 = c()
p2 = c()
for(i in names(ShortestPath)){
  p1 = c(p1, ShortestPath[[i]][["pos"]])
  p2 = c(p2, ShortestPath[[i]][["neg"]])
}


Pvals = c()
for(i in 0:(max(max(p1), max(p2))-1)){
  p_pos = p1; p_neg = p2
  p_pos = p_pos>i
  p_neg= p_neg>i
  Table = matrix(c(sum(!p_pos), sum(p_pos), sum(!p_neg), sum(p_neg)), nrow = 2, byrow = TRUE)
  F1 = fisher.test(Table, alternative = "greater")
  Pvals = c(Pvals, F1$p.value)
}
qval = p.adjust(Pvals, method = "BH")




