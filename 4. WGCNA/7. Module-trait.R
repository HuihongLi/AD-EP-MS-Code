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
load("data/3.module definition.RData")
source("tutorialFunctions.R")  

#=====================================================================================

pdf("plot/F11-module-trait-correlation.pdf",height=8,width=7) 
par(mfrow=c(1,1), mar=c(5, 9,4, 0) + 0.1, cex=1) 
#module-Trait association
PCs1A = moduleEigengenes(expDataEP,  colors=modulesA1)  
ME_1A = PCs1A$eigengenes 

datTraitsEP = datTraitsEP

MEs = orderMEs(ME_1A)
moduleTraitCor = cor(MEs, datTraitsEP, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(datTraitsEP[,1]))
moduleTraitCor_0.65 = data.frame(moduleTraitCor[abs(moduleTraitCor)>0.65,])

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

col = colorRampPalette(c("navy", "white", "red"))(50);

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraitsEP),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = col,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("EP Module-trait relationships"))

# Recalculate MEs with new color labels
PCs2A = moduleEigengenes(expDataTG,  colors=modulesB2_new)  
ME_2A = PCs2A$eigengenes 

colnames(datTraitsTG) = c("TG Interaction")

MEs = orderMEs(ME_2A)
moduleTraitCor = cor(MEs, datTraitsTG, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(datTraitsTG[,1]))
moduleTraitCor_0.65 = data.frame(moduleTraitCor[abs(moduleTraitCor)>0.65,])
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraitsTG),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = col,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("TG Module-trait relationships"))

# Recalculate MEs with new color labels
PCs3A = moduleEigengenes(expDataJ20,  colors=modulesB3_new)  
ME_3A = PCs3A$eigengenes 

colnames(datTraitsJ20) = c("J20 Interaction")

MEs = orderMEs(ME_3A)
moduleTraitCor = cor(MEs, datTraitsJ20, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, length(datTraitsJ20[,1]))
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraitsJ20),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = col,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("J20 Module-trait relationships"))
dev.off()

