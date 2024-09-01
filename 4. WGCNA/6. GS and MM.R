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
PCs1A = moduleEigengenes(expDataEP,  colors=modulesA1)  
ME_1A = PCs1A$eigengenes 
colorsA1 = names(table(modulesA1)) 
geneModuleMembership1 = as.data.frame(cor(expDataEP, ME_1A, use = "p"));
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep="");  
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(expDataEP)[[2]]);  
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep=""); 
Gene = rownames(t(expDataEP)) 
kMEtable1  = cbind(Gene,Gene,modulesA1) 
for (i in 1:length(colorsA1)) 
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i]) 
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), 
                                                    colnames(MMPvalue1)))) 

Epileptic = as.data.frame(datTraitsEP$Epileptic);
names(Epileptic) = "Epileptic"
geneTraitSignificanceEP = as.data.frame(cor(expDataEP, Epileptic, use = "p"));
GSPvalueEP = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceEP), length(expDataEP[,1])));
names(geneTraitSignificanceEP) = paste("GS.", names(Epileptic), sep="");
names(GSPvalueEP) = paste("p.GS.", names(Epileptic), sep="")

kMEtable1 = cbind(kMEtable1, geneTraitSignificanceEP, GSPvalueEP)
write.csv(kMEtable1,"result/GS_MM_EP.csv",row.names=FALSE) 


pdf("plot/F9-GS.MM_EP.pdf",height=5,width=5) 
par(bg = "white")  

for (which.module in names(table(modulesA1))) { 
  module = which.module
  column = match(module, colorsA1);
  moduleGenes = modulesA1==module;
  verboseScatterplot(abs(geneModuleMembership1[moduleGenes, column]),
                     abs(geneTraitSignificanceEP[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for EP",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
#=====================================================================================
PCs2A = moduleEigengenes(expDataTG,  colors=modulesB2_new)  
ME_2A = PCs2A$eigengenes 
colorsB2 = names(table(modulesB2_new))
geneModuleMembership2 = as.data.frame(cor(expDataTG, ME_2A, use = "p"));
colnames(geneModuleMembership2)=paste("PC",colorsB2,".cor",sep="");
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(expDataTG)[[2]]);
colnames(MMPvalue2)=paste("PC",colorsB2,".pval",sep="");
Gene = rownames(t(expDataTG))
kMEtable2  = cbind(Gene,Gene,modulesB2_new)
for (i in 1:length(colorsB2)) 
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership2), 
                                                    colnames(MMPvalue2))))
Interaction = as.data.frame(datTraitsTG$Interaction)
names(Interaction) = "Interaction"
geneTraitSignificanceTG = as.data.frame(cor(expDataTG, Interaction, use = "p"));
GSPvalueTG = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceTG), length(expDataTG[,1])));
names(geneTraitSignificanceTG) = paste("GS.", names(Interaction), sep="")
names(GSPvalueTG) = paste("p.GS.", names(Interaction), sep="")
kMEtable2 = cbind(kMEtable2, geneTraitSignificanceTG, GSPvalueTG)
write.csv(kMEtable2,"result/GS_MM_TG.csv",row.names=FALSE)

pdf("plot/F9-GS.MM_TG.pdf",height=5,width=5)
par(bg = "white")  

for (which.module in names(table(modulesB2_new))) { 
  module = which.module
  column = match(module, colorsB2);
  moduleGenes = modulesB2_new==module;
  verboseScatterplot(abs(geneModuleMembership2[moduleGenes, column]),
                     abs(geneTraitSignificanceTG[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for TG",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
#=====================================================================================
PCs3A = moduleEigengenes(expDataJ20,  colors=modulesB3_new)
ME_3A = PCs3A$eigengenes
colorsB3 = names(table(modulesB3_new))
geneModuleMembership3 = as.data.frame(cor(expDataJ20, ME_3A, use = "p"));
colnames(geneModuleMembership3)=paste("PC",colorsB3,".cor",sep="");
MMPvalue3=corPvalueStudent(as.matrix(geneModuleMembership3),dim(expDataJ20)[[2]]);
colnames(MMPvalue3)=paste("PC",colorsB3,".pval",sep="");
Gene = rownames(t(expDataJ20))
kMEtable3  = cbind(Gene,Gene,modulesB3_new)
for (i in 1:length(colorsB3)) 
  kMEtable3 = cbind(kMEtable3, geneModuleMembership3[,i], MMPvalue3[,i])
colnames(kMEtable3)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership3), 
                                                    colnames(MMPvalue3))))

Interaction = as.data.frame(datTraitsJ20$Interaction)
names(Interaction) = "Interaction"
geneTraitSignificanceJ20 = as.data.frame(cor(expDataJ20, Interaction, use = "p"));
GSPvalueJ20 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificanceJ20), length(expDataJ20[,1])));
names(geneTraitSignificanceJ20) = paste("GS.", names(Interaction), sep="")
names(GSPvalueJ20) = paste("p.GS.", names(Interaction), sep="")
kMEtable3 = cbind(kMEtable3, geneTraitSignificanceJ20, GSPvalueJ20)
write.csv(kMEtable3,"result/GS_MM_J20.csv",row.names=FALSE)

pdf("plot/F9-GS.MM_J20.pdf",height=5,width=5)
par(bg = "white")  

for (which.module in names(table(modulesB3_new))) { 
  module = which.module
  column = match(module, colorsB3);
  moduleGenes = modulesB3_new==module;
  verboseScatterplot(abs(geneModuleMembership3[moduleGenes, column]),
                     abs(geneTraitSignificanceJ20[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for J20",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
#=====================================================================================
pdf("plot/F9-GS_EP_TG.pdf",height=5,width=5)
par(bg = "white")  

for (which.module in names(table(modulesA1))) { 
  module = which.module
  column = match(module, colorsB3);
  moduleGenes = modulesA1==module;
  verboseScatterplot(abs(geneTraitSignificanceEP[moduleGenes, 1]),
                     abs(geneTraitSignificanceTG[moduleGenes, 1]),
                     xlab = paste("Gene significance for EP", module, "module"),
                     ylab = "Gene significance for TG",
                     main = paste(" Gene significance EP vs. gene significance TG\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
pdf("plot/F9-GS_EP_J20.pdf",height=5,width=5)
par(bg = "white")  

for (which.module in names(table(modulesA1))) { 
  module = which.module
  column = match(module, colorsB3);
  moduleGenes = modulesA1==module;
  verboseScatterplot(abs(geneTraitSignificanceEP[moduleGenes, 1]),
                     abs(geneTraitSignificanceJ20[moduleGenes, 1]),
                     xlab = paste("Gene significance for EP", module, "module"),
                     ylab = "Gene significance for J20",
                     main = paste("Gene significance EP vs. gene significance J20\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
#=====================================================================================

geneTraitCorColorEP <- numbers2colors(geneTraitSignificanceEP)
geneTraitCorColorTGInter <- numbers2colors(geneTraitSignificanceTG)
geneTraitCorColorJ20Inter <- numbers2colors(geneTraitSignificanceJ20)

pdf("plot/F10-Renew_dendrogram_plotswith trait.pdf",height=5,width=15) 
plotDendroAndColors(geneTreeEP, data.frame(modulesA1,geneTraitCorColorEP), data.frame("Module colors","EP"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="EP dendrogram") 
plotDendroAndColors(geneTreeTG, data.frame(modulesB2_new,modulesA1,geneTraitCorColorEP,geneTraitCorColorTGInter), data.frame("Module colors","EP Module colors","EP","TG Interaction"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="TG dendrogram") 
plotDendroAndColors(geneTreeJ20, data.frame(modulesB3_new,modulesA1,geneTraitCorColorEP,geneTraitCorColorJ20Inter), data.frame("Module colors","EP Module colors","EP","J20 Interaction"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="J20 dendrogram")
dev.off()
save(geneModuleMembership1,geneModuleMembership2,geneModuleMembership3, file = "data/5.geneModuleMembership.RData")
save(geneTraitCorColorEP,geneTraitCorColorTGInter,geneTraitCorColorJ20Inter, file = "data/6.geneTraitCorColor.RData")