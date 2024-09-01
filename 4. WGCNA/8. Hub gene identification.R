library(WGCNA);
library(flashClust);
library(qvalue);
library(Hmisc);
library(impute);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
enableWGCNAThreads()
options(stringsAsFactors = FALSE);
set.seed(123)
collectGarbage();
load("data/1.expandtraitData.RData")
load("data/3.module definition.RData")
load("data/2.geneTree.RData")
load("data/5.geneModuleMembership.RData")
load("data/6.geneTraitCorColor.RData")
source("tutorialFunctions.R")  

#=====================================================================================
#Export Network and its connectivity
source("tutorialFunctions.R")   
colorsA1 = names(table(modulesA1)) 
for (co in colorsA1[colorsA1!="grey"]) 
  visantPrepOverall(modulesA1, co, expDataEP, colnames(expDataEP), 500, softPower, TRUE)

for (co in colorsA1[colorsA1!="grey"]) 
  visantPrepOverall(modulesA1, co, expDataTG, colnames(expDataTG), 500, softPower, TRUE)

for (co in colorsA1[colorsA1!="grey"]) 
  visantPrepOverall(modulesA1, co, expDataJ20, colnames(expDataTG), 500, softPower, TRUE)

#Hub genes specific to one network
datExprA12g = t(cbind(t(expDataEP),t(expDataTG)))
i1 = 1:dim(t(expDataEP))[[2]]; 
i2 = (1:dim(t(expDataTG))[[2]])+length(i1) 
for (co in colorsA1[colorsA1!="grey"]) 
  visantPrep(modulesA1, co, i1, i2, datExprA12g, colnames(expDataEP), 500, softPower, TRUE)

#Hub genes specific to one network
datExprA12g = t(cbind(t(expDataEP),t(expDataJ20)))
i1 = 1:dim(t(expDataEP))[[2]]; 
i2 = (1:dim(t(expDataJ20))[[2]])+length(i1) 
for (co in colorsA1[colorsA1!="grey"]) 
  visantPrep(modulesA1, co, i1, i2, datExprA12g, colnames(expDataEP), 500, softPower, TRUE)

colorsB2 = names(table(modulesB2_new))
for (co in colorsB2[colorsB2!="grey"]) 
  visantPrepOverall(modulesB2_new, co, expDataTG, colnames(expDataTG), 500, softPower, TRUE)
colorsB3 = names(table(modulesB3_new))
for (co in colorsB3[colorsB3!="grey"])
  visantPrepOverall(modulesB3_new, co, expDataJ20, colnames(expDataJ20), 500, softPower, TRUE)



#=====================================================================================
sigresultsEP <- read.csv("data/sigresultsEP.csv", header = TRUE, sep = ",")
sigresultsEP <- sigresultsEP[,2]
sigresultsTG <- read.csv("data/Tg4510_sig_interaction-ONLY_STATS.csv", header = TRUE, sep = ",")
sigresultsTG_filter <- sigresultsTG[sigresultsTG$GenotyperTg4510.Age_months8 > 1 | sigresultsTG$GenotyperTg4510.Age_months8 < -1,]
sigresultsTG_filter <- sigresultsTG_filter$X
shareDEG <- intersect(sigresultsEP, sigresultsTG_filter)

EPcolor <- rep("grey",length(colnames(expDataEP)))
for (i in 1:length(colnames(expDataEP))){
  if (colnames(expDataEP)[i] %in% sigresultsEP){
    EPcolor[i] <- "red"
  }
}

TGcolor <- rep("grey",length(colnames(expDataEP)))

for (i in 1:length(colnames(expDataEP))){
  if (colnames(expDataEP)[i] %in% sigresultsTG_filter){
    TGcolor[i] <- "red"
  }
}


pdf("plot/F12-RenewDEG_dendrogram_plotswith trait.pdf",height=5,width=15) 
plotDendroAndColors(geneTreeEP, data.frame(modulesA1,geneTraitCorColorEP,EPcolor), data.frame("Module colors","EP association","EP-DEG"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="EP dendrogram") 
plotDendroAndColors(geneTreeTG, data.frame(modulesB2_new,modulesA1,geneTraitCorColorEP,EPcolor,geneTraitCorColorTGInter,TGcolor), data.frame("Module colors","EP Module colors","EP association","EP-DEG","TG association","TG-DEG"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="TG dendrogram") 
plotDendroAndColors(geneTreeJ20, data.frame(modulesB3_new,modulesA1,geneTraitCorColorEP,EPcolor,geneTraitCorColorJ20Inter), data.frame("Module colors","EP Module colors","EP","EP-DEG","J20 association"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="J20 dendrogram")
dev.off()
#=====================================================================================
#Hub gene identification on EP module
path <- "Network/EP/"
file_list <- list.files(path = path, pattern = "*_connectivityOverall.csv", full.names = TRUE)

data_list <- list()
for(file in file_list) {
  data <- read.csv(file)
  data_list[[basename(file)]] <- data
}

process_data <- function(data) {
  n <- nrow(data)
  data <- data[1:(0.10 * n), ]
  return(data)
}

data_list <- lapply(data_list, process_data)

data_list <- lapply(data_list, function(x) {
  x <- x[["probes"]]
  return(x)
})

GS_MMEP <- read.csv("result/GS_MM_EP.csv")
GS_MMEP <- GS_MMEP[GS_MMEP$Module != "grey", ]
GS_MMEP$GS.Epileptic
hub_list <- list()
for (which.module in names(table(GS_MMEP$Module))) { 
  module_gene = GS_MMEP[GS_MMEP$Module == which.module, "Gene"]
  module_MM =  GS_MMEP[GS_MMEP$Module == which.module, paste("PC",which.module,".cor",sep = "")]
  module_GS =  GS_MMEP[GS_MMEP$Module == which.module, "GS.Epileptic"]
  data = data.frame("Gene" = module_gene, "MM" = module_MM, "GS" = module_GS) 
  data = data[abs(data$MM) > 0.9 & abs(data$GS) > 0.4, "Gene"] 
  hub_list[[which.module]] = data
}

names(data_list) <- names(hub_list)

Final_hub_list <- list()
for (module in names(data_list)) {
  x <- data_list[[module]]
  y <- hub_list[[module]]
  Final_hub_list[[module]] <- intersect(x, y)
  
}
max_length <- max(sapply(Final_hub_list, length))
Final_hub_list_df <- do.call(cbind, lapply(Final_hub_list, function(x) {
  length(x) <- max_length
  return(x)
}))
colnames(Final_hub_list_df) <- names(Final_hub_list)
write.csv(Final_hub_list_df, "result/Final_hub_listEP.csv", row.names = FALSE)
DEGintersect <- list()

for (module in names(Final_hub_list)) {
  x <- Final_hub_list[[module]]
  DEGintersect[[module]] <- intersect(x, sigresultsEP)
}
max_length <- max(sapply(DEGintersect, length))
DEGintersect_df <- do.call(cbind, lapply(DEGintersect, function(x) {
  length(x) <- max_length
  return(x)
}))
colnames(DEGintersect_df) <- names(DEGintersect)

write.csv(DEGintersect_df, "result/DEGintersectEP.csv", row.names = FALSE)
HubEP <- DEGintersect
#=====================================================================================
#Hub gene identification on TG module
path <- "Network/TG/"
file_list <- list.files(path = path, pattern = "*_connectivityOverall.csv", full.names = TRUE)
data_list <- list()
for(file in file_list) {
  data <- read.csv(file)
  data_list[[basename(file)]] <- data
}

process_data <- function(data) {
  n <- nrow(data)
  data <- data[1:(0.10 * n), ]
  return(data)
}

data_list <- lapply(data_list, process_data)

data_list <- lapply(data_list, function(x) {
  x <- x[["probes"]]
  return(x)
})

GS_MMTG <- read.csv("result/GS_MM_TG.csv")
GS_MMTG <- GS_MMTG[GS_MMTG$Module != "grey", ]
hub_list <- list()
for (which.module in names(table(GS_MMTG$Module))) { 
  module_gene = GS_MMTG[GS_MMTG$Module == which.module, "Gene"]
  module_MM =  GS_MMTG[GS_MMTG$Module == which.module, paste("PC",which.module,".cor",sep = "")]
  module_GS =  GS_MMTG[GS_MMTG$Module == which.module, "GS.Interaction"]
  data = data.frame("Gene" = module_gene, "MM" = module_MM, "GS" = module_GS) 
  data = data[abs(data$MM) > 0.9 & abs(data$GS) > 0.4, "Gene"] 
  hub_list[[which.module]] = data
}

names(data_list) <- names(hub_list)

Final_hub_list <- list()
for (module in names(data_list)) {
  x <- data_list[[module]]
  y <- hub_list[[module]]
  Final_hub_list[[module]] <- intersect(x, y)
  
}

max_length <- max(sapply(Final_hub_list, length))
Final_hub_list_df <- do.call(cbind, lapply(Final_hub_list, function(x) {
  length(x) <- max_length
  return(x)
}))
colnames(Final_hub_list_df) <- names(Final_hub_list)
write.csv(Final_hub_list_df, "result/Final_hub_listTG.csv", row.names = FALSE)
DEGintersect <- list()
for (module in names(Final_hub_list)) {
  x <- Final_hub_list[[module]]
  DEGintersect[[module]] <- intersect(x, sigresultsTG_filter)
}
max_length <- max(sapply(DEGintersect, length))
DEGintersect_df <- do.call(cbind, lapply(DEGintersect, function(x) {
  length(x) <- max_length
  return(x)
}))
colnames(DEGintersect_df) <- names(DEGintersect)

write.csv(DEGintersect_df, "result/DEGintersectTG.csv", row.names = FALSE)
HubTG <- DEGintersect

#=====================================================================================
#share hub of EP and TG
share_EP_TG_turquoise <- intersect(HubEP$turquoise, HubTG$turquoise)
write.csv(share_EP_TG_turquoise, "result/share_EP_TG_turquoise.csv", row.names = FALSE)
#=====================================================================================
#Hub gene identification on J20 module
path <- "Network/J20/"
file_list <- list.files(path = path, pattern = "*_connectivityOverall.csv", full.names = TRUE)
data_list <- list()
for(file in file_list) {
  data <- read.csv(file)
  data_list[[basename(file)]] <- data
}

process_data <- function(data) {
  n <- nrow(data)
  data <- data[1:(0.10 * n), ]
  return(data)
}

data_list <- lapply(data_list, process_data)

data_list <- lapply(data_list, function(x) {
  x <- x[["probes"]]
  return(x)
})

GS_MMJ20 <- read.csv("result/GS_MM_J20.csv")
GS_MMJ20 <- GS_MMJ20[GS_MMJ20$Module != "grey", ]
hub_list <- list()
for (which.module in names(table(GS_MMJ20$Module))) { 
  module_gene = GS_MMJ20[GS_MMJ20$Module == which.module, "Gene"]
  module_MM =  GS_MMJ20[GS_MMJ20$Module == which.module, paste("PC",which.module,".cor",sep = "")]
  module_GS =  GS_MMJ20[GS_MMJ20$Module == which.module, "GS.Interaction"]
  data = data.frame("Gene" = module_gene, "MM" = module_MM, "GS" = module_GS) 
  data = data[abs(data$MM) > 0.9 & abs(data$GS) > 0.4, "Gene"] 
  hub_list[[which.module]] = data
}

names(data_list) <- names(hub_list)

Final_hub_list <- list()
for (module in names(data_list)) {
  x <- data_list[[module]]
  y <- hub_list[[module]]
  Final_hub_list[[module]] <- intersect(x, y)
  
}

max_length <- max(sapply(Final_hub_list, length))
Final_hub_list_df <- do.call(cbind, lapply(Final_hub_list, function(x) {
  length(x) <- max_length
  return(x)
}))
colnames(Final_hub_list_df) <- names(Final_hub_list)
write.csv(Final_hub_list_df, "result/Final_hub_listJ20.csv", row.names = FALSE)
#=====================================================================================
