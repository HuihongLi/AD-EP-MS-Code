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
#=====================================================================================
# Network construction and module detection
adjacencyJ20 = adjacency(expDataJ20,power=softPower,type="signed"); 
diag(adjacencyJ20)=0 
#load("data/TOMJ20.RData")
dissTOMJ20  = 1-TOMsimilarity(adjacencyJ20, TOMType="signed") 
save(dissTOMJ20, file = "data/TOMJ20.RData")
geneTreeJ20  = flashClust(as.dist(dissTOMJ20), method="average") 


adjacencyTG = adjacency(expDataTG,power=softPower,type="signed");
diag(adjacencyTG)=0
#load("data/TOMTG.RData")
dissTOMTG  = 1-TOMsimilarity(adjacencyTG, TOMType="signed")
save(dissTOMTG, file = "data/TOMTG.RData")
geneTreeTG  = flashClust(as.dist(dissTOMTG), method="average")

adjacencyEP = adjacency(expDataEP,power=softPower,type="signed");
diag(adjacencyEP)=0
#load("data/TOMEP.RData")
dissTOMEP  = 1-TOMsimilarity(adjacencyEP, TOMType="signed")
save(dissTOMEP, file = "data/TOMEP.RData")
geneTreeEP  = flashClust(as.dist(dissTOMEP), method="average")

pdf("plot/F4-dendrogram.pdf",height=6,width=12) 
plot(geneTreeJ20,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity J20", 
     labels=FALSE,hang=0.04); 
plot(geneTreeTG,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity TG", 
     labels=FALSE,hang=0.04); 
plot(geneTreeEP,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity EP", 
     labels=FALSE,hang=0.04);  
dev.off() 

#=====================================================================================
#deepsplit choose
mColorhEP=NULL 
for (ds in 0:3){ 
  tree = cutreeHybrid(dendro = geneTreeEP, pamStage=FALSE, 
                      minClusterSize = (30-3*ds), cutHeight = 0.99,  
                      deepSplit = ds, distM = dissTOMEP) 
  mColorhEP=cbind(mColorhEP,labels2colors(tree$labels)); 
} 

mColorhTG=NULL
for (ds in 0:3){ 
  tree = cutreeHybrid(dendro = geneTreeTG, pamStage=FALSE, 
                      minClusterSize = (30-3*ds), cutHeight = 0.99,  
                      deepSplit = ds, distM = dissTOMTG) 
  mColorhTG=cbind(mColorhTG,labels2colors(tree$labels)); 
} 

mColorhJ20=NULL 
for (ds in 0:3){ 
  tree = cutreeHybrid(dendro = geneTreeJ20, pamStage=FALSE, 
                      minClusterSize = (30-3*ds), cutHeight = 0.99,  
                      deepSplit = ds, distM = dissTOMJ20) 
  mColorhJ20=cbind(mColorhJ20,labels2colors(tree$labels)); 
} 

save(mColorhEP, mColorhTG, mColorhJ20, file = "data/2.mColor.RData")
save(geneTreeJ20,geneTreeTG,geneTreeEP, file = "data/2.geneTree.RData")

pdf("plot/F5-Module_choices.pdf", height=10,width=25);  
plotDendroAndColors(geneTreeEP, mColorhEP, paste("dpSplt =",0:3), main = "Gene dendrogram and module colors EP",dendroLabels=FALSE); 
plotDendroAndColors(geneTreeTG, mColorhTG, paste("dpSplt =",0:3), main = "Gene dendrogram and module colors TG",dendroLabels=FALSE); 
plotDendroAndColors(geneTreeJ20, mColorhJ20, paste("dpSplt =",0:3), main = "Gene dendrogram and module colors J20",dendroLabels=FALSE); 
dev.off() 
#=====================================================================================