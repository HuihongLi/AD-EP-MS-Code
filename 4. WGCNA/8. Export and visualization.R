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
load("data/TOMEP.RData")
load("data/TOMTG.RData")
load("data/TOMJ20.RData")
#=====================================================================================

# TURQUOISE
TOM = -dissTOMEP + 1
probes = colnames(expDataEP)
module = "turquoise"
# Takes gene names in the module and matches them 
inModule = is.finite(match(modulesA1, module))
# Find the probes in the module of interest
modProbes = probes[inModule]
# Select gene names which are turquoise
#modGenes_turquoise = geneInfoALL_turquoise$probeIDs
modTom = TOM[inModule,inModule]
dimnames(modTom) = list(modProbes)
dim(modTom)
# Calculate treshold -> top 50 nodes (genes) for adjacency
diag(modTom) <- 0
a <- sort(modTom, decreasing=TRUE)
percent <- 1
repeat {
  limit <- nrow(modTom) * percent
  thresh <- a[limit]
  b <- apply(modTom >= thresh, 1, any)
  if (nrow(modTom[b, b]) <= 50) {
    break
  } else {
    percent <- percent - 0.005
  }
}
# treshold corresponds to treshold for weight. We selected top 50 nodes (genes) values for weight
# The top 50 nodes (genes) for network interconnectedness (or proximity, measured as topological overlap measure or TOM)
thresh # 0.3776317

# Export the network into edge and node list files for Cytoscape, set the threshold for how dense/sparse you want the network to be. Threshold is the weighting where the edges below will be removed from the network.
cyt = exportNetworkToCytoscape(modTom,
                               edgeFile = paste("CytoEdgeEP", paste(module, collapse = "-"), ".txt", sep = ""),
                               nodeFile = paste("CytoNodeEP",paste(module,collapse="-"),".txt",sep=""),
                               weighted = TRUE,
                               threshold = thresh,
                               nodeNames = modProbes,
                               altNodeNames = modProbes,
                               nodeAttr = modulesA1[inModule]) ################################################### thickeness -> in the software!

# TURQUOISE in TG
TOM = -dissTOMTG + 1

#TOM = TOMsimilarityFromExpr(expDataEP, power=8) # power = soft-thresholding power for network construction.
#save(TOM, file = "Tg4510_RNAseq_WGCNA_TOM.RData")

#load(file = "Tg4510_RNAseq_WGCNA_TOM.RData")

# Select module which you're interested in from the TOM block
probes = colnames(expDataTG)

module = "turquoise"
# Takes gene names in the module and matches them 
inModule = is.finite(match(modulesA1, module))
# Find the probes in the module of interest
modProbes = probes[inModule]
# Select gene names which are turquoise
#modGenes_turquoise = geneInfoALL_turquoise$probeIDs
modTom = TOM[inModule,inModule]
dimnames(modTom) = list(modProbes)
dim(modTom)
# Calculate treshold -> top 50 nodes (genes) for adjacency
diag(modTom) <- 0
a <- sort(modTom, decreasing=TRUE)
percent <- 1
repeat {
  limit <- nrow(modTom) * percent
  thresh <- a[limit]
  b <- apply(modTom >= thresh, 1, any)
  if (nrow(modTom[b, b]) <= 50) {
    break
  } else {
    percent <- percent - 0.005
  }
}
# treshold corresponds to treshold for weight. We selected top 50 nodes (genes) values for weight
# The top 50 nodes (genes) for network interconnectedness (or proximity, measured as topological overlap measure or TOM)
thresh # 0.3776317

# Export the network into edge and node list files for Cytoscape, set the threshold for how dense/sparse you want the network to be. Threshold is the weighting where the edges below will be removed from the network.
cyt = exportNetworkToCytoscape(modTom,
                               edgeFile = paste("CytoEdgeAD", paste(module, collapse = "-"), ".txt", sep = ""),
                               nodeFile = paste("CytoNodeAD",paste(module,collapse="-"),".txt",sep=""),
                               weighted = TRUE,
                               threshold = thresh,
                               nodeNames = modProbes,
                               altNodeNames = modProbes,
                               nodeAttr = modulesA1[inModule]) ################################################### thickeness -> in the software!
