library(FactoMineR)
library(factoextra)
library(DESeq2)
library(ggplot2)
library(ggrepel)


count_matrix <- read.csv("EPmerge_remove.csv", row.names = 1)
count_matrix_filtered <- count_matrix[rowSums(count_matrix != 0) > 0, ]
metadata <- read.csv("EPmeta.csv", header = TRUE, row.names = "Sample", stringsAsFactors = TRUE, check.names = FALSE)
metadata$Group <- factor(metadata$Group)
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered, colData = metadata, design = ~Group)


vst.deseq <- vst(dds, blind = TRUE)  # transformed data
vst1 <- as.data.frame(assay(vst.deseq)) # corrected to assay function


pca_resultsvst1 <- PCA(t(vst1), graph = FALSE)


pca_scores <- as.data.frame(pca_resultsvst1$ind$coord)  
mahalanobis_distance <- mahalanobis(pca_scores, colMeans(pca_scores), cov(pca_scores))
cutoff <- qchisq(0.975, df = ncol(pca_scores))  
outliers <- which(mahalanobis_distance > cutoff)
pca_scores$outlier <- ifelse(1:nrow(pca_scores) %in% outliers, "Outlier", "Normal")


ggplot(pca_scores, aes(x = Dim.1, y = Dim.2, color = outlier)) +
  geom_point(aes(shape = outlier), size = 3) +
  scale_color_manual(values = c("Normal" = "gray", "Outlier" = "red")) +
  scale_shape_manual(values = c("Normal" = 1, "Outlier" = 19)) +
  geom_text_repel(data = subset(pca_scores, outlier == "Outlier"), aes(label = rownames(pca_scores[outliers,])), size = 3) +
  ggtitle("PCA Analysis on EP dataset - Outliers highlighted") +
  theme_minimal()

removelist <- names(outliers)

count1 <- count_matrix_filtered[ , !(names(count_matrix_filtered) %in% removelist)]
write.csv(count1, "EPmerge_remove_afterPCA.csv", row.names = TRUE)



write.csv(outliers, "PCA_outliers.csv")

write.csv(mahalanobis_distance, "PCA_mahalanobis_distance.csv")




