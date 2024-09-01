library(WGCNA);
library(flashClust);
library(qvalue);
library(Hmisc);
library(impute);
enableWGCNAThreads()
options(stringsAsFactors = FALSE);
collectGarbage();
set.seed(123)

#=====================================================================================

femDataJ20 = read.csv("data/norm_J20.csv", header = TRUE, row.names = 1);
femDataTG = read.csv("data/norm_TG.csv", header = TRUE, row.names = 1);
femDataEP = read.csv("data/norm_EP.csv", header = TRUE, row.names = 1);

m.mad <- apply(femDataJ20,1,mad)
femDataJ20 <- femDataJ20[which(m.mad > 
                                 max(quantile(m.mad, probs=seq(0, 1, 0.15))[2],0.01)),]

m.mad <- apply(femDataTG,1,mad)
femDataTG <- femDataTG[which(m.mad > 
                               max(quantile(m.mad, probs=seq(0, 1, 0.15))[2],0.01)),]

m.mad <- apply(femDataEP,1,mad)
femDataEP <- femDataEP[which(m.mad > 
                               max(quantile(m.mad, probs=seq(0, 1, 0.15))[2],0.01)),]


# common gene
commonProbesA = intersect(rownames(femDataEP),rownames(femDataJ20))
commonProbesA = intersect(commonProbesA,rownames(femDataTG))
femDataEP = femDataEP[commonProbesA,]
femDataJ20 = femDataJ20[commonProbesA,]
femDataTG = femDataTG[commonProbesA,]

#=====================================================================================


expDataJ20 = as.data.frame(t(femDataJ20));
expDataTG = as.data.frame(t(femDataTG));
expDataEP = as.data.frame(t(femDataEP));


#=====================================================================================

gsg = goodSamplesGenes(expDataJ20, verbose = 3);
gsg$allOK
gsg = goodSamplesGenes(expDataTG, verbose = 3);
gsg$allOK
gsg = goodSamplesGenes(expDataEP, verbose = 3);
gsg$allOK


#=====================================================================================
sampleTreeJ20 = hclust(dist(expDataJ20), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTreeJ20, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 36, col = "red");
clust = cutreeStatic(sampleTreeJ20, cutHeight = 36, minSize = 10)
table(clust)
keepSamplesJ20 = (clust==1)
expDataJ20 = expDataJ20[keepSamplesJ20, ]

sampleTreeTG = hclust(dist(expDataTG), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTreeTG, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 40, col = "red");
clust = cutreeStatic(sampleTreeTG, cutHeight = 40, minSize = 10)
table(clust)
keepSamplesTG = (clust==1)
expDataTG = expDataTG[keepSamplesTG, ]

sampleTreeEP = hclust(dist(expDataEP), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTreeEP, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 45, col = "red");
clust = cutreeStatic(sampleTreeEP, cutHeight = 45, minSize = 10)
table(clust)
keepSamplesEP = (clust==1)
expDataEP = expDataEP[keepSamplesEP, ]


#=====================================================================================
pdf(file = "Plot/F1-sampleClustering.pdf", width = 12, height = 9);
traitDataEP = read.csv("data/metaEP.csv");
traitData3 = traitDataEP
traitData3$Outlier = as.numeric(!keepSamplesEP)
traitColors2 = numbers2colors(traitData3$Epileptic, signed = FALSE) 
traitColors1 = numbers2colors(traitData3$Outlier, signed = FALSE)
datColorsEP = data.frame(traitColors1, traitColors2)
names(datColorsEP) = c("Outlier", "Epileptic")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTreeEP, groupLabels = names(datColorsEP), 
                    colors = datColorsEP, 
                    main = "Sample dendrogram and trait heatmap of EP dataset")

traitDataTG = read.csv("data/metaTG.csv");
traitData2 = traitDataTG

traitData2$Outlier = as.numeric(!keepSamplesTG)
traitColors2 = numbers2colors(traitData2$Interaction, signed = FALSE) 
traitColors1 = numbers2colors(traitData2$Outlier, signed = FALSE)
datColorsTG = data.frame(traitColors1, traitColors2)
names(datColorsTG) = c("Outlier", "Interaction")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTreeTG, groupLabels = names(datColorsTG), 
                    colors = datColorsTG, 
                    main = "Sample dendrogram and trait heatmap of TG dataset")

traitDataJ2O = read.csv("data/metaJ20.csv");
traitData1 = traitDataJ2O

traitData1$Outlier = as.numeric(!keepSamplesJ20)
traitColors2 = numbers2colors(traitData1$Interaction, signed = FALSE) 
traitColors1 = numbers2colors(traitData1$Outlier, signed = FALSE)
datColorsJ20 = data.frame(traitColors1, traitColors2)
names(datColorsJ20) = c("Outlier", "Interaction")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTreeJ20, groupLabels = names(datColorsJ20), 
                    colors = datColorsJ20, 
                    main = "Sample dendrogram and trait heatmap of J20 dataset")

dev.off()


#=====================================================================================
#Trait dataframe

SamplesJ20 = rownames(expDataJ20);
traitRows = match(SamplesJ20, traitDataJ2O$Sample);
datTraitsJ20 = as.data.frame(traitDataJ2O[traitRows, -1]);
colnames(datTraitsJ20) = colnames(traitDataJ2O)[2];

SamplesTG = rownames(expDataTG);
traitRows = match(SamplesTG, traitDataTG$Sample);
datTraitsTG = as.data.frame(traitDataTG[traitRows, -1]);
colnames(datTraitsTG) = colnames(traitDataTG)[2];

SamplesEP = rownames(expDataEP);
traitRows = match(SamplesEP, traitDataEP$Sample);
datTraitsEP = as.data.frame(traitDataEP[traitRows, -1]);
colnames(datTraitsEP) = colnames(traitDataEP)[2];

collectGarbage();

#=====================================================================================

powers = c(c(3:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expDataJ20, powerVector = powers, verbose = 5)
sftgraph = data.frame(soft = sft$fitIndices[,1], scale = -sign(sft$fitIndices[,3])*sft$fitIndices[,2],conn = sft$fitIndices[,5],group="J20")

sft = pickSoftThreshold(expDataTG, powerVector = powers, verbose = 5)
sftgraph = rbind(sftgraph,data.frame(soft = sft$fitIndices[,1], scale = -sign(sft$fitIndices[,3])*sft$fitIndices[,2],conn = sft$fitIndices[,5],group="TG"))

sft = pickSoftThreshold(expDataEP, powerVector = powers, verbose = 5)
sftgraph = rbind(sftgraph,data.frame(soft = sft$fitIndices[,1], scale = -sign(sft$fitIndices[,3])*sft$fitIndices[,2],conn = sft$fitIndices[,5],group="EP"))

library(ggplot2)
p1 <- ggplot(sftgraph, aes(x=soft, y=scale, color=group)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept=0.85) +
  scale_y_continuous(breaks=seq(0, 1, 0.05)) +
  ggtitle("Scale independence") +
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  geom_text(aes(label=soft), hjust=0, vjust=1.2, show.legend = FALSE) +
  theme(legend.position = "none")
p2 <- ggplot(sftgraph, aes(x=soft, y=conn, color=group)) +
  geom_point() +
  theme_bw() +
  ggtitle("Mean connectivity") +
  xlab("Soft Threshold (power)") +
  ylab("Mean connectivity") +
  geom_text(aes(label=soft), hjust=0, vjust=1.2, show.legend = FALSE) + 
  labs(color="Group") +
  theme(legend.position = c(0.98, 0.98), 
        legend.justification = c("right", "top")) 
library(gridExtra)
p3 <-grid.arrange(p1, p2, ncol=2)
ggsave("Plot/F2-soft-thresholding.pdf", p3, width = 12, height = 7)

#=====================================================================================
#check compraibility of networks

softPower = 9;
rankExprJ20= rank(rowMeans(t(expDataJ20)))
rankExprTG= rank(rowMeans(t(expDataEP)))
rankExprEP= rank(rowMeans(t(expDataTG)))
rankConnJ20= rank(softConnectivity(expDataJ20,type="signed",power=softPower)) 
rankConnEP= rank(softConnectivity(expDataEP,type="signed",power=softPower)) 
rankConnTG= rank(softConnectivity(expDataTG,type="signed",power=softPower)) 
pdf("plot/F3-compara-generalNetworkProperties.pdf", height=12, width=12)
par(mfrow=c(2,2)) 
verboseScatterplot(rankExprJ20,rankExprEP, xlab="Ranked Expression J20",  
                   ylab="Ranked Expression EP") 
verboseScatterplot(rankConnJ20,rankConnEP, xlab="Ranked Connectivity J20",  
                   ylab="Ranked Connectivity EP") 
verboseScatterplot(rankExprTG,rankExprEP, xlab="Ranked Expression TG",  
                   ylab="Ranked Expression EP") 
verboseScatterplot(rankConnTG,rankConnEP, xlab="Ranked Connectivity TG",  
                   ylab="Ranked Connectivity EP") 
dev.off()


save(expDataJ20, expDataTG, expDataEP,softPower,datTraitsJ20,datTraitsTG,datTraitsEP, file = "data/1.expandtraitData.RData")
