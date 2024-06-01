Data = readRDS("Data/Validation/Data_DREAMER3.rds")

Data1 = list()
for(d in names(Data)){
    Data1[[d]] = Data[[d]][4:7]
    }

adj_Pr_after_Correction = readRDS("Data/Validation/adjustedProteins_after_correction.rds")
length(adj_Pr_after_Correction)
length(names(Data1))
I = names(Data1)[names(Data1) %in% names(adj_Pr_after_Correction)]
length(I)
Data1 = Data1[I]
names(Data1)


saveRDS(Data1, "Data/Validation/Data_DREAMER3_withoutScores.rds")
