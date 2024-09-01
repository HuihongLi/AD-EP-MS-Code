library(WGCNA);
library(flashClust);
library(qvalue);
library(Hmisc);
library(impute);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
enableWGCNAThreads()
options(stringsAsFactors = FALSE);
set.seed(123)
collectGarbage();
load("data/1.expandtraitData.RData")
load("data/2.mColor.RData")
load("data/2.geneTree.RData")
load("data/3.module definition.RData")
source("tutorialFunctions.R")  

#=====================================================================================
# calculation of Eigengenes and module membership in EP module
PCs1A = moduleEigengenes(expDataEP,  colors=modulesA1)  
ME_1A = PCs1A$eigengenes 
Epileptic = datTraitsEP$Epileptic
ME_1Atrait = orderMEs(cbind(ME_1A, Epileptic))
pdf("plot/F8-ModuleEigengeneVisualizationsEP.pdf",height=12,width=15) 
plotEigengeneNetworks(ME_1Atrait, "", marDendro = c(0,4,1,2), marHeatmap = c(5,5,1,1), cex.lab = 0.8, xLabelsAngle
                      = 90)

MEexpress = as.data.frame(PCs1A$averageExpr)
MEexpress$Sample = 1:nrow(MEexpress)
MEexpressmelt <- melt(MEexpress, id.vars = "Sample")


background_data <- data.frame(
  xmin = c(0, 93),
  xmax = c(93, 180),
  ymin = -Inf, 
  ymax = Inf,  
  fill = c("Control", "Epileptic")  
)

for (which.module in names(table(modulesA1))) { 
  MEexpressmelt$color <- ifelse(substring(MEexpressmelt$variable, 3) == which.module, which.module, "#BEBEBE46")

  p <- ggplot(MEexpressmelt, aes(x = Sample, y = value, group = variable, color = color)) +
    geom_rect(data = background_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
              alpha = 0.3, inherit.aes = FALSE) + 
    geom_line(size = 1) +
    scale_color_identity() +
    scale_fill_manual(values = c("Control" = "lightblue", "Epileptic" = "lightgreen"), name = "Condition") + 
    labs(y = paste("averageExpression in module ", which.module), x = "") +
    theme_minimal() +
    theme(legend.position = "top")  
  
  print(p)
}
dev.off()




#=====================================================================================
load("data/1.expandtraitData.RData")
# calculation of Eigengenes and module membership in EP modulem
PCs2A = moduleEigengenes(expDataTG,  colors=modulesB2_new)  
ME_2A = PCs2A$eigengenes 
TG = datTraitsTG$Interaction
ME_2Atrait = orderMEs(cbind(ME_2A, TG))
pdf("plot/F8-ModuleEigengeneVisualizationsTG.pdf",height=12,width=15) 
plotEigengeneNetworks(ME_2Atrait, "", marDendro = c(0,4,1,2), marHeatmap = c(5,5,1,1), cex.lab = 0.8, xLabelsAngle
                      = 90)
datTraitsTG$Sample = rownames(ME_2Atrait)
order_index = order(datTraitsTG$Interaction)

datTraitsTG = datTraitsTG[order_index,]
rownames(datTraitsTG) = 1:nrow(datTraitsTG)
ME_2A = ME_2A[order_index,]

#ME_1A = ME_1A[order_index,]
#ME_1A$Sample = rownames(ME_1A)
#ME_1A$Epileptic = datTraitsEP$Epilepti

MEexpress = as.data.frame(PCs2A$averageExpr)
MEexpress = MEexpress[order_index,]

MEexpress$time = 1:nrow(MEexpress)
MEexpressmelt <- melt(MEexpress, id.vars = "time")


background_data <- data.frame(
  xmin = c(0, 29,38,45,50), 
  xmax = c(29,38,45,50,57), 
  ymin = -Inf, 
  ymax = Inf,  
  fill = c("Control", "TG2","TG4","TG6","TG8") 
)

for (which.module in names(table(modulesB2_new))) { 
  MEexpressmelt$color <- ifelse(substring(MEexpressmelt$variable, 3) == which.module, which.module, "#BEBEBE46")
  # 绘图
  p <- ggplot(MEexpressmelt, aes(x = time, y = value, group = variable, color = color)) +
    geom_rect(data = background_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
              alpha = 0.3, inherit.aes = FALSE) +
    geom_line(size = 1) +
    scale_color_identity() +
    scale_fill_manual(values = c("Control" = "lightblue", "TG2" = "lightgreen","TG4"= "lightcyan","TG6"="#F781BF", "TG8" = "#A65628"), name = "Condition") + 
    labs(y = paste("averageExpression in module ", which.module), x = "") +
    theme_minimal() +
    theme(legend.position = "top")  
  

  print(p)
}


dev.off()

#=====================================================================================
load("data/1.expandtraitData.RData")
# calculation of Eigengenes and module membership in EP modulem
PCs3A = moduleEigengenes(expDataJ20,  colors=modulesB3_new)  
ME_3A = PCs3A$eigengenes 
J20 = datTraitsJ20$Interaction
ME_3Atrait = orderMEs(cbind(ME_3A, J20))
pdf("plot/F8-ModuleEigengeneVisualizationsJ20.pdf",height=12,width=15) 
plotEigengeneNetworks(ME_3Atrait, "", marDendro = c(0,4,1,2), marHeatmap = c(5,5,1,1), cex.lab = 0.8, xLabelsAngle
                      = 90)
datTraitsJ20$Sample = rownames(ME_3Atrait)
order_index = order(datTraitsJ20$Interaction)

datTraitsJ20 = datTraitsJ20[order_index,]
rownames(datTraitsJ20) = 1:nrow(datTraitsJ20)
ME_3A = ME_3A[order_index,]

MEexpress = as.data.frame(PCs3A$averageExpr)
MEexpress = MEexpress[order_index,]

MEexpress$time = 1:nrow(MEexpress)
MEexpressmelt <- melt(MEexpress, id.vars = "time")
background_data <- data.frame(
  xmin = c(0, 30,38,45,52), 
  xmax = c(30,38,45,52,59), 
  ymin = -Inf, 
  ymax = Inf, 
  fill = c("Control", "J20_6","J20_8","J20_10","J20_12") 
)

for (which.module in names(table(modulesB3_new))) { 
  MEexpressmelt$color <- ifelse(substring(MEexpressmelt$variable, 3) == which.module, which.module, "#BEBEBE46")
  # 绘图
  p <- ggplot(MEexpressmelt, aes(x = time, y = value, group = variable, color = color)) +
    geom_rect(data = background_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), 
              alpha = 0.3, inherit.aes = FALSE) + 
    geom_line(size = 1) +
    scale_color_identity() +
    scale_fill_manual(values = c("Control" = "lightblue", "J20_6" = "lightgreen","J20_8"= "lightcyan","J20_10"="#F781BF", "J20_12" = "#A65628"), name = "Condition") + 
    labs(y = paste("averageExpression in module ", which.module), x = "") +
    theme_minimal() +
    theme(legend.position = "top")  
  

  print(p)
}


dev.off()

#=====================================================================================
