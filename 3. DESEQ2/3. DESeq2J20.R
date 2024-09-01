


# sessionInfo()

library("DESeq2")
library("RUVSeq")
library("EDASeq")
library("genefilter")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("geneplotter")


# column data / phenotypes
coldata <- read.csv("data/metaJ20.csv", row.names=1, stringsAsFactors=FALSE)
colnames(coldata)
rownames(coldata)

coldata$Age_months <- as.factor(coldata$Age_months)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$Genotype <- relevel(coldata$Genotype, "WT")
levels(coldata$Genotype)

countdata <- read.csv("data/ADmerge_remove_J20.csv", header=TRUE, row.names=1)
rownames(countdata)
colnames(countdata)

# Convert to matrix
countdata <- as.matrix(countdata) # count matrix raw reads
head(countdata)

# design: ~condition + time + condition:time (time = age; condition = genotype)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~Genotype + Age_months + Genotype:Age_months)
# We have 2 factors: age (with values 6, 8, 10, and 12) and genotype (with values WT and J20)

# Pre-filtering the data set
# Remove the rows that have no or nearly no information about the amount of gene expression
# Removing rows of the DESeqDataSet that have no counts, or only a single count across all samples
nrow(dds)
dds <- dds[ rowSums(counts(dds)>= 20) > 6, ] # filtering for non-expressed and lowly expressed genes (sum of expression accross all samples at least 6 because final n = 6-8 mice per group)
nrow(dds)

##### Differential expression analysis

### BULK: WALD TEST ###
dds_Wald <- DESeq(dds, test="Wald")
# The Wald test (also called the Wald Chi-Squared Test) is a way to find out if explanatory variables in a model are significant.
# “Significant” means that they add something to the model; variables that add nothing can be deleted without affecting the model in any meaningful way.
# The test can be used for a multitude of different models including those with binary variables or continuous variables.
# The Wald test is a rough approximation of the Likelihood Ratio Test. However, you can run it with a single model (the LR test requires at least two).
# It is also more broadly applicable than the LRT: often, you can run a Wald in situations where no other test can be run.
# For large values of n, the Wald test is roughly equivalent to the t-test; both tests will reject the same values for large sample sizes.
# The Wald, LRT and Lagrange multiplier tests are all equivalent as sample sizes approach infinity (called “asymptotically equivalent”).
# However, samples of a finite size, especially smaller samples, are likely to give very different results.

# replacing outliers and refitting for 9 genes
dds_Wald # class: DESeqDataSet

resultsNames(dds_Wald)

# I want the coefficients and p-values for:
# Genotype_J20_vs_WT -> effect of genotype
# Age_months_8_vs_6, Age_months_10_vs_6, Age_months_12_vs_6 -> each time point compared to age 6 (effect of age)
# GenotypeJ20.Age_months12, GenotypeJ20.Age_months10, GenotypeJ20.Age_months12 -> example: is GenotypeJ20.Age_months12 different than GenotypeJ20.Age_months10?

# Building the results table
# Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula
res_Wald <- results(dds_Wald)
res_Wald
results_Wald_frame <- data.frame(res_Wald, stringsAsFactors=FALSE)
write.csv(results_Wald_frame, file = "J20_res_Wald.csv") # WRONG P-VALUE (Wald test p-value: GenotypeJ20.Age months8)

stats_Wald <- mcols(dds_Wald)
rownames(stats_Wald) <- rownames(dds_Wald)
head(stats_Wald)
stats_Wald_frame <- data.frame(stats_Wald, stringsAsFactors=FALSE)

# DESeq preserves the original counts in counts(dds) saving the replacement counts as a matrix named replaceCounts in assays(dds).
replacedcounts_Wald <- assays(dds_Wald)[["replaceCounts"]] # 
head(replacedcounts_Wald) # class(): matrix


### PROGRESSIVE CHANGES (genes which show a diff expression profile over time across genotype): LRT TEST ###
dds_LRT <- DESeq(dds, reduced=~Genotype + Age_months, test="LRT") # DESeq: standard differential expression analysis steps are wrapped in this single function
# LRT: Likelihood-ratio test
# The Likelihood-Ratio test (sometimes called the likelihood-ratio chi-squared test) is a hypothesis test that helps you choose the “best” model between two nested models.
# “Nested models” means that one is a special case of the other.

# replacing outliers and refitting for 9 genes
dds_LRT # class: DESeqDataSet

resultsNames(dds_LRT)

# Building the results table
# Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula
res_LRT <- results(dds_LRT)
res_LRT # pvalue = interaction
results_LRT_frame <- data.frame(res_LRT, stringsAsFactors=FALSE)
write.csv(results_LRT_frame, file = "J20_res_LRT.csv")

stats_LRT <- mcols(dds_LRT)
head(stats_LRT)
stats_LRT_frame <- data.frame(stats_LRT, stringsAsFactors=FALSE)


### EFFECT OF AGE: MANUAL LRT ###
full <- stats::model.matrix.default(~Genotype*Age_months, data = as.data.frame(colData(dds)))
head(full)
reduced <- full[,-c(3:5)]
head(reduced)

dds_age <- DESeq(dds, full=full, reduced=reduced, test="LRT")

# replacing outliers and refitting for 9 genes
dds_age # class: DESeqDataSet

resultsNames(dds_age)

# Building the results table
# Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula
res_age <- results(dds_age)
res_age # pvalue = what I am interested (???)
results_age_frame <- data.frame(res_age, stringsAsFactors=FALSE)
write.csv(results_age_frame, file = "J20_res_age.csv")

stats_age <- mcols(dds_age)
head(stats_age)
stats_age_frame <- data.frame(stats_age, stringsAsFactors=FALSE)

save(countdata, coldata, dds, file="DEseq2_object_counts_J20.RData") # saves as objects; to open use load()
save(dds_Wald, res_Wald, stats_Wald, file="DEseq2_Wald_test_results_J20.RData")
save(dds_LRT, res_LRT, stats_LRT, file="DEseq2_LRT_test_results_J20.RData")
save(dds_age, res_age, stats_age, file="DEseq2_modified_LRT_test_results_J20.RData")

# Stats Tables (final)
FDR_adj_genotype <- p.adjust(stats_Wald[,"WaldPvalue_Genotype_J20_vs_WT"], method = "fdr") # we have calculated the FDR-adj pvalue ourselves because 1) they were not available for all analysis,and 2) from DESeq some were "NA"
FDR_adj_age <- p.adjust(stats_age[,"LRTPvalue"], method = "fdr")
FDR_adj_LRT <- p.adjust(stats_LRT[,"LRTPvalue"], method = "fdr")

stats_table <-cbind(FDR_adj_genotype, as.data.frame(stats_Wald[,c("WaldPvalue_Genotype_J20_vs_WT", "Genotype_J20_vs_WT")]), 
                    FDR_adj_age,
                    res_age[,"pvalue"],	
                    as.data.frame(stats_Wald[,c("WaldPvalue_Age_months_8_vs_6","Age_months_8_vs_6","WaldPvalue_Age_months_10_vs_6","Age_months_10_vs_6",
                                                "WaldPvalue_Age_months_12_vs_6","Age_months_12_vs_6")]), 
                    FDR_adj_LRT,
                    res_LRT[,"pvalue"],  
                    as.data.frame(stats_Wald[,c("WaldPvalue_GenotypeJ20.Age_months8","GenotypeJ20.Age_months8","WaldPvalue_GenotypeJ20.Age_months10","GenotypeJ20.Age_months10",
                                                "WaldPvalue_GenotypeJ20.Age_months12","GenotypeJ20.Age_months12")]))
rownames(stats_table) <- rownames(res_LRT)
write.csv(stats_table, file = "J20_stats_table.csv")

sig_genotype <- stats_table[which(stats_table[,"FDR_adj_genotype"]<0.05),]
write.csv(sig_genotype, file = "J20_sig_genotype.csv")

sig_age <- stats_table[which(stats_table[,"FDR_adj_age"]<0.05),]
write.csv(sig_age, file = "J20_sig_age.csv")

sig_LRT <- stats_table[which(stats_table[,"FDR_adj_LRT"]<0.05),]
write.csv(sig_LRT, file = "J20_sig_LRT.csv")

### New Table (Log2FC, etc)
## Genotype
table_publication_genotype <- cbind(FDR_adj_genotype, as.data.frame(stats_Wald[,c("WaldPvalue_Genotype_J20_vs_WT", "Genotype_J20_vs_WT", "WaldStatistic_Genotype_J20_vs_WT")]))
dim(table_publication_genotype)
rownames(table_publication_genotype) <- rownames(res_Wald)
head(table_publication_genotype)

genes_J20_genotype <- rownames(sig_genotype)
length(genes_J20_genotype)

table_publication_genotype <- table_publication_genotype[genes_J20_genotype,]
head(table_publication_genotype)

write.csv(table_publication_genotype, file = "J20_sig_genotype_STATS.csv")

# binomial test
threshold <- 0.05
ntotal <- nrow(table_publication_genotype) # number of sig genes I identified
ndownregulated <- nrow(table_publication_genotype[which(table_publication_genotype[,"Genotype_J20_vs_WT"]<0),])
test <- binom.test(ndownregulated, ntotal, 0.5)
pvalue <- test[3]
pvalue <- signif(as.numeric(pvalue), 3)


## Progression (LRT)
table_publication_interaction <- cbind(FDR_adj_LRT, as.data.frame(stats_LRT[,c("GenotypeJ20.Age_months8", "GenotypeJ20.Age_months10", "GenotypeJ20.Age_months12", "LRTStatistic", "LRTPvalue")]))
dim(table_publication_interaction)
rownames(table_publication_interaction) <- rownames(res_LRT)
head(table_publication_interaction)

genes_J20_interaction <- rownames(sig_LRT)
length(genes_J20_interaction)

table_publication_interaction <- table_publication_interaction[genes_J20_interaction,]
head(table_publication_interaction)

write.csv(table_publication_interaction, file = "J20_sig_interaction_STATS.csv")


### Filter genes that show an age effect, from the interaction list from DESeq2
common <- intersect(rownames(sig_age), rownames(sig_LRT))
interaction_only <- table_publication_interaction[!rownames(table_publication_interaction) %in% common, ]
write.csv(interaction_only, file = "J20_sig_interaction-ONLY_STATS.csv")


