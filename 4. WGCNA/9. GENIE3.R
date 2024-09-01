library(GENIE3)
hubgenelist <- c("Ggta1","Lyz2","Fcgr2b","Tyrobp","AU020206","Cd300c2","C4b","S100a6","C1qa","Gfap","C1qc","Osmr","Trem2","Tgif1","C1qb","Adgre1","Apobec1","Csf3r","Cxcl16")
set.seed(123) # For reproducibility of results
ADdata <- read.csv("data/ADmerge_remove_TG.csv", header = TRUE, row.names = 1)
EPdata <- read.csv("data/EPmerge_remove_afterPCA.csv", header = TRUE, row.names = 1) 


modulelistEP <- read.csv("data/genes_and_modulesEP.csv", header = TRUE)
modulelistEP <- modulelistEP[modulelistEP$Module == "turquoise",]
modulelistEP <- modulelistEP$Genes

modulelistAD <- read.csv("data/genes_and_modulesTG.csv", header = TRUE)
modulelistAD <- modulelistAD[modulelistAD$Module == "turquoise",]
modulelistAD <- modulelistAD$Genes

ADdata <- ADdata[modulelistAD,]
EPdata <- EPdata[modulelistEP,]


count_normAD <- as.matrix(ADdata)
weightMatAD <- GENIE3(count_normAD, nCores = 7)
linkListAD <- getLinkList(weightMatAD,threshold = 0.01)
linkListtopAD <- linkListAD[linkListAD$regulatoryGene %in% hubgenelist & linkListAD$targetGene %in% hubgenelist,]
linkListtopAD <- linkListtopAD[1:30,]
write.table(linkListtopAD, file = "linkListAD.txt", sep = "\t", quote = FALSE, row.names = FALSE)



count_normEP <- as.matrix(EPdata)
weightMatEP <- GENIE3(count_normEP, nCores = 7)
linkListEP <- getLinkList(weightMatEP,threshold = 0.01)
linkListtopEP <- linkListEP[linkListEP$regulatoryGene %in% hubgenelist & linkListEP$targetGene %in% hubgenelist,]
linkListtopEP <- linkListtopEP[1:30,]
write.table(linkListtopEP, file = "linkListEP.txt", sep = "\t", quote = FALSE, row.names = FALSE)





