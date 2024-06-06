rm(list = ls())
library(ggplot2)
library(dplyr)

Pval_ADRDP_Proteins = readRDS("data/Pval_ADRDP_Proteins.rds")
Final_ADRs = readRDS("data/Final_67_ADRs_STRING_AfterControlFor_OrganRemoval_and_IndicationDiffused.rds")
intersected_ADRs = readRDS("data/Intersected_ADRs_Physical_STRING.rds")

Pval_ADRDP_Proteins = Pval_ADRDP_Proteins[Final_ADRs[,1],]
rownames(Pval_ADRDP_Proteins) = Final_ADRs[,2]

I = order(Pval_ADRDP_Proteins[,1], decreasing = TRUE)

Pval_ADRDP_Proteins_ranked = Pval_ADRDP_Proteins[I,]

ADR_Pval = data.frame(
  Effect = rownames(Pval_ADRDP_Proteins_ranked),
  P_Value = Pval_ADRDP_Proteins_ranked[,2])

#########################################################
pdf("data/FiguresS3_ADRs_sorted_by_pvals.pdf", width = 4, height = 7)

fill = ifelse(ADR_Pval[,1] %in% intersected_ADRs, "red", "royalblue")
ggplot() +
  geom_point(data = ADR_Pval, aes(x=reorder(Effect, P_Value, FUN = function(x) -x), y=-log10(P_Value)),
             size=3, color="black", fill = fill, stroke = .7, shape=21) +
  coord_flip() + # Flips the axes to have ADRs on the y-axis
  theme_minimal(base_size = 10) +
  labs(
    y = "-log10(p-value)", 
    x = "Phenotypes"
  ) +
  theme(
    panel.grid.major.x = element_blank(), # Remove vertical major grid lines
    panel.grid.minor.x = element_blank(), # Remove vertical minor grid lines
    panel.grid.major.y = element_line(size=.1, color="gray"),  # Keep horizontal lines
    panel.grid.minor.y = element_line(size=.1, color="gray"),  # Keep horizontal lines
    plot.title = element_text(face="bold", size=14),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12)
  )

dev.off()

