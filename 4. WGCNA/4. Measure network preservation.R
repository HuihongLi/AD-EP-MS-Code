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

#=====================================================================================
modulesA1 =  mColorhEP[,3]
# (This step will take ~10-30 minutes) 
multiExpr  = list(A1=list(data=expDataEP),A2=list(data=expDataTG)) 
multiColor = list(A1 = modulesA1) 
#mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed", 
#                      nPermutations=100,maxGoldModuleSize=2000,maxModuleSize=2000) 
#save mp
#save(mp,file="data/modulePreservationEP_TG.RData")
load("data/modulePreservationEP_TG.RData")
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2 
stats[order(-stats[,2]),c(1:2)] 
stats2 = mp[["preservation"]][["observed"]][["ref.A1"]][["inColumnsAlsoPresentIn.A2"]]
write.csv(stats2[order(-stats2[,2]),c(1:2)] ,"result/modulePreservationEP_TG_Rank.csv",row.names=TRUE)
write.csv(stats[order(-stats[,2]),c(1:2)] ,"result/modulePreservationEP_TG.csv",row.names=TRUE)

multiExpr  = list(A1=list(data=expDataEP),A2=list(data=expDataJ20)) 
multiColor = list(A1 = modulesA1) 
#mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed", 
#                      nPermutations=100,maxGoldModuleSize=2000,maxModuleSize=2000) 
#save mp
#save(mp,file="data/modulePreservationEP_J20.RData")
load("data/modulePreservationEP_J20.RData")
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2 
stats[order(-stats[,2]),c(1:2)] 
stats2 = mp[["preservation"]][["observed"]][["ref.A1"]][["inColumnsAlsoPresentIn.A2"]]
write.csv(stats2[order(-stats2[,2]),c(1:2)] ,"result/modulePreservationEP_J20_Rank.csv",row.names=TRUE)
write.csv(stats[order(-stats[,2]),c(1:2)] ,"result/modulePreservationEP_J20.csv",row.names=TRUE)

#=====================================================================================
#genes in between-network comparisons.   
colorsA1 = names(table(modulesA1)) 
PCs1A = moduleEigengenes(expDataEP,  colors=modulesA1)  
ME_1A = PCs1A$eigengenes 

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


#  Calculate  MEs of modulesA1 in TG and J20
PCs2A = moduleEigengenes(expDataTG,  colors=modulesA1)  
ME_2A = PCs2A$eigengenes 
geneModuleMembership2 = signedKME(expDataTG, ME_2A) 
colnames(geneModuleMembership2)=paste("PC",colorsA1,".cor",sep="");  
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(expDataTG)[[2]]);  
colnames(MMPvalue2)=paste("PC",colorsA1,".pval",sep=""); 
kMEtable2  = cbind(Gene,Gene,modulesA1) 
for (i in 1:length(colorsA1)) 
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i]) 
colnames(kMEtable2)=colnames(kMEtable1) 
write.csv(kMEtable2,"result/MM_EPinTG.csv",row.names=FALSE)


PCs3A = moduleEigengenes(expDataJ20,  colors=modulesA1)  
ME_3A = PCs3A$eigengenes 
geneModuleMembership3 = signedKME(expDataJ20, ME_3A) 
colnames(geneModuleMembership3)=paste("PC",colorsA1,".cor",sep="");  
MMPvalue3=corPvalueStudent(as.matrix(geneModuleMembership3),dim(expDataJ20)[[2]]);  
colnames(MMPvalue3)=paste("PC",colorsA1,".pval",sep=""); 
kMEtable3  = cbind(Gene,Gene,modulesA1) 
for (i in 1:length(colorsA1)) 
  kMEtable3 = cbind(kMEtable3, geneModuleMembership3[,i], MMPvalue3[,i]) 
colnames(kMEtable3)=colnames(kMEtable1) 
write.csv(kMEtable3,"result/MM_EPinJ20.csv",row.names=FALSE)

#=====================================================================================
pdf("plot/F7_1_all_kMETG_vs_kMEEP.pdf",height=8,width=8) 
par(mfrow=c(1,1), mar=c(6, 6, 6, 6) + 0.1, cex=1) 
for (c in 1:length(colorsA1)){ 
  verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsA1[c], 
                     xlab="kME in TG",ylab="kME in EP") 
}; dev.off() 

pdf("plot/F7_2_inModule_kMETG_vs_kMEEP.pdf",height=8,width=8) 
par(mfrow=c(1,1), mar=c(6, 6, 6, 6) + 0.1, cex=1) 
for (c in 1:length(colorsA1)){ 
  inMod = modulesA1== colorsA1[c] 
  verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsA1[c], 
                     xlab="kME in TG",ylab="kME in EP") 
}; dev.off() 

pdf("plot/F7_3_all_kMEJ20_vs_kMEEP.pdf",height=8,width=8) 
par(mfrow=c(1,1), mar=c(6, 6, 6, 6) + 0.1, cex=1) 
for (c in 1:length(colorsA1)){ 
  verboseScatterplot(geneModuleMembership3[,c],geneModuleMembership1[,c],main=colorsA1[c], 
                     xlab="kME in J20",ylab="kME in EP") 
}; dev.off() 

pdf("plot/F7_4_inModule_kMEJ20_vs_kMEEP.pdf",height=8,width=8) 
par(mfrow=c(1,1), mar=c(6, 6, 6, 6) + 0.1, cex=1) 
for (c in 1:length(colorsA1)){ 
  inMod = modulesA1== colorsA1[c] 
  verboseScatterplot(geneModuleMembership3[inMod,c],geneModuleMembership1[inMod,c],main=colorsA1[c], 
                     xlab="kME in J20",ylab="kME in EP") 
}; dev.off() 

#=====================================================================================


topGenesKME1 = NULL 
for (c in 1:length(colorsA1)){ 
  kMErank1    = rank(-geneModuleMembership1[,c]) 
  kMErank2    = rank(-geneModuleMembership2[,c]) 
  maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max)) 
  topGenesKME1 = cbind(topGenesKME1,Gene[maxKMErank<=50]) 
}; colnames(topGenesKME1) = colorsA1 
topGenesKME1
write.csv(topGenesKME1,"result/topGenesKME1.csv",row.names=FALSE)
topGenesKME2 = NULL 
for (c in 1:length(colorsA1)){ 
  kMErank1    = rank(-geneModuleMembership1[,c]) 
  kMErank2    = rank(-geneModuleMembership3[,c]) 
  maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max)) 
  topGenesKME2 = cbind(topGenesKME2,Gene[maxKMErank<=50]) 
}; colnames(topGenesKME2) = colorsA1 
topGenesKME2 
write.csv(topGenesKME2,"result/topGenesKME2.csv",row.names=FALSE)

#=====================================================================================
library(ggplot2)
library(patchwork)


moduleper1 <- read.csv("result/modulePreservationEP_TG.csv")
moduleper2 <- read.csv("result/modulePreservationEP_J20.csv")
moduleper1_rank <- read.csv("result/modulePreservationEP_TG_Rank.csv")
moduleper2_rank <- read.csv("result/modulePreservationEP_J20_Rank.csv")


custom_theme <- theme(
  plot.background = element_blank(),  
  panel.background = element_blank(),  
  panel.border = element_rect(color = "black", fill = NA), 
  plot.border = element_blank(),       
  panel.grid = element_blank(),        
  axis.line = element_line(color = "black") 
)


p1 <- ggplot(moduleper1, aes(x = module_size, y = Zsummary.pres, color = module_name)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray") +
  labs(title = "Zsummary.pres of EP module in TG", x = "Module Size", y = "Zsummary.pres") +
  scale_color_identity() +
  custom_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


p2 <- ggplot(moduleper1_rank, aes(x = moduleSize, y = medianRank.pres, color = X)) +
  geom_point(size = 3) +
  labs(title = "MedianRank.pres of EP module in TG", x = "Module Size", y = "MedianRank.pres") +
  scale_color_identity() +
  custom_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p3 <- ggplot(moduleper2, aes(x = module_size, y = Zsummary.pres, color = module_name)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray") +
  labs(title = "Zsummary.pres of EP module in J20", x = "Module Size", y = "Zsummary.pres") +
  scale_color_identity() +
  custom_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


p4 <- ggplot(moduleper2_rank, aes(x = moduleSize, y = medianRank.pres, color = X)) +
  geom_point(size = 3) +
  labs(title = "MedianRank.pres of EP module in J20", x = "Module Size", y = "MedianRank.pres") +
  scale_color_identity() +
  custom_theme +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


combined_plot <- (p2 | p1) / (p4 | p3) 


ggsave("plot/F7_modulePreservationEP.pdf", plot = combined_plot, height = 16, width = 16)
