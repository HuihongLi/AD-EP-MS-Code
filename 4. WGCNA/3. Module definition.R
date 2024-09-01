library(WGCNA);
library(flashClust);
library(qvalue);
library(Hmisc);
library(impute);
enableWGCNAThreads()
options(stringsAsFactors = FALSE);
set.seed(123)
collectGarbage();
load("data/1.expandtraitData.RData")
load("data/2.mColor.RData")
load("data/2.geneTree.RData")
source("tutorialFunctions.R")  
#=====================================================================================
# calculation of Eigengenes and module membership in EP module
modulesA1 =  mColorhEP[,3]
# calculation of Eigengenes
PCs1A = moduleEigengenes(expDataEP,  colors=modulesA1)  
ME_1A = PCs1A$eigengenes 
colorsA1 = names(table(modulesA1)) 
#calculation of geneModuleMembership1 using Eigengenes  
geneModuleMembership1 = signedKME(expDataEP, ME_1A) 
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep="");  
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(expDataEP)[[2]]);  
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep=""); 
Gene = rownames(t(expDataEP)) 
kMEtable1  = cbind(Gene,Gene,modulesA1) 
for (i in 1:length(colorsA1)) 
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i]) 
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), 
                                                    colnames(MMPvalue1)))) 
write.csv(kMEtable1,"result/MM_EP.csv",row.names=FALSE) 
#=====================================================================================

#resigned the modules color according to modulesA1
modulesB2 =  mColorhTG[,3]
modulesB2_new = matchModules(Gene, modulesA1, Gene, modulesB2) 
dataframemodulesA1 = t(rbind(Gene,modulesA1))
enrichmentsB2A1 = userListEnrichment(Gene,modulesB2_new,"result/MM_EP.csv","EP","result/enrichmentTG_EP.csv") # show the overlap of module

modulesB3 =  mColorhJ20[,3]
modulesB3_new = matchModules(Gene, modulesA1, Gene, modulesB3) 
enrichmentsB3A1 = userListEnrichment(Gene,modulesB3_new,"result/MM_EP.csv","EP","result/enrichmentJ20_EP.csv") 

pdf("plot/F6-Renew_dendrogram_plots.pdf",height=5,width=15) 
plotDendroAndColors(geneTreeEP, data.frame(modulesA1), data.frame("Module colors"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="EP dendrogram") 
plotDendroAndColors(geneTreeTG, data.frame(modulesB2_new,modulesA1), data.frame("Module colors","EP-Module colors"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="TG dendrogram")   
plotDendroAndColors(geneTreeJ20, data.frame(modulesB3_new,modulesA1), data.frame("Module colors","EP-Module colors"), dendroLabels=FALSE, hang=0.03,  
                    addGuide=TRUE, guideHang=0.05, main="J20 dendrogram")   
dev.off()


genes_and_modulesJ20 <- data.frame(Gene,modulesB3_new)
names(genes_and_modulesJ20) <- c("Genes", "Modules")
write.csv(genes_and_modulesJ20, "result/genes_and_modulesJ20.csv", row.names = FALSE)

genes_and_modulesTG <- data.frame(Gene,modulesB2_new)
names(genes_and_modulesTG) <- c("Genes", "Modules")
write.csv(genes_and_modulesTG, "result/genes_and_modulesTG.csv", row.names = FALSE)

genes_and_modulesEP <- data.frame(Gene,modulesA1)
names(genes_and_modulesEP) <- c("Genes", "Modules")
write.csv(genes_and_modulesEP, "result/genes_and_modulesEP.csv", row.names = FALSE)
#=====================================================================================

library(ggplot2)
library(gridExtra) 
genenumberEP <- data.frame(table(modulesA1))
genenumberTG <- data.frame(table(modulesB2_new))
genenumberJ20 <- data.frame(table(modulesB3_new))

color_vector1 <- names(table(modulesA1))
color_vector2 <- names(table(modulesB2_new))
color_vector3 <- names(table(modulesB3_new))

genenumberEP <- genenumberEP[order(genenumberEP$Freq, decreasing = TRUE), ]
genenumberTG <- genenumberTG[order(genenumberTG$Freq, decreasing = TRUE), ]
genenumberJ20 <- genenumberJ20[order(genenumberJ20$Freq, decreasing = TRUE), ]

genenumberEP <- genenumberEP[genenumberEP$modulesA1 != "grey", ]
genenumberTG <- genenumberTG[genenumberTG$modulesB2_new != "grey", ]
genenumberJ20 <- genenumberJ20[genenumberJ20$modulesB3_new != "grey", ]


p1 <- ggplot(genenumberEP, aes(x=reorder(modulesA1, Freq), y=Freq)) +
  geom_bar(stat="identity", fill=color_vector1[match(genenumberEP$modulesA1, color_vector1)]) +
  coord_flip() +
  theme_minimal() + 
  xlab("Module") + 
  ylab("Number of genes") + 
  ggtitle("EP")

p2 <- ggplot(genenumberTG, aes(x=reorder(modulesB2_new, Freq), y=Freq)) +
  geom_bar(stat="identity", fill=color_vector2[match(genenumberTG$modulesB2_new, color_vector2)]) +
  coord_flip() +
  theme_minimal() + 
  xlab("Module") + 
  ylab("Number of genes") + 
  ggtitle("TG")

p3 <- ggplot(genenumberJ20, aes(x=reorder(modulesB3_new, Freq), y=Freq)) +
  geom_bar(stat="identity", fill=color_vector3[match(genenumberJ20$modulesB3_new, color_vector3)]) +
  coord_flip() +
  theme_minimal() + 
  xlab("Module") + 
  ylab("Number of genes") + 
  ggtitle("J20")

p3 <- grid.arrange(p1, p2, p3, ncol=3)
ggsave("Plot/F6.1-genenummodule.pdf", p3, width = 12, height = 7)

save(modulesA1, modulesB2_new, modulesB3_new, file = "data/3.module definition.RData")


