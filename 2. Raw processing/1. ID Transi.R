count <- read.table("EP raw.txt", header = TRUE, sep = "\t", row.names = 1)
colnames(count) <- gsub("\\.bam", "", colnames(count))

count <- count[, -c(1:5)]
library(dplyr)
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
my_ensmusg_ids <- rownames(count)
genes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
               filters='ensembl_gene_id', 
               values=my_ensmusg_ids, 
               mart=mart)

genes <- genes[genes$external_gene_name != "",]
count$Geneid <- rownames(count)

merged_data <- merge(count, genes, by.x='Geneid', by.y='ensembl_gene_id', all.x=TRUE)

merged_data$Geneid <- ifelse(is.na(merged_data$external_gene_name), merged_data$Geneid, merged_data$external_gene_name)

final_data <- merged_data[, -which(names(merged_data) %in% 'external_gene_name')]

sum(is.na(final_data$Geneid))
sum(final_data$Geneid == "")

library(dplyr)

final_data <- final_data %>% group_by(Geneid) %>% summarise_all(sum)
write.csv(final_data, "EPmerge.csv", row.names = FALSE)



count <- read.table("AD raw.txt", header = TRUE, sep = "\t", row.names = 1)
colnames(count) <- gsub("\\.bam", "", colnames(count))
count <- count[, -c(1:5)]
my_ensmusg_ids <- rownames(count)
genes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
               filters='ensembl_gene_id', 
               values=my_ensmusg_ids, 
               mart=mart)

genes <- genes[genes$external_gene_name != "",]
count$Geneid <- rownames(count)

merged_data <- merge(count, genes, by.x='Geneid', by.y='ensembl_gene_id', all.x=TRUE)

merged_data$Geneid <- ifelse(is.na(merged_data$external_gene_name), merged_data$Geneid, merged_data$external_gene_name)

final_data <- merged_data[, -which(names(merged_data) %in% 'external_gene_name')]

sum(is.na(final_data$Geneid))
sum(final_data$Geneid == "")

library(dplyr)

final_data <- final_data %>% group_by(Geneid) %>% summarise_all(sum)
write.csv(final_data, "ADmerge.csv", row.names = FALSE)


count <- read.csv("EPmerge.csv")
removelist <- c("ERR1779321","ERR1779322","ERR1779324","ERR1779325","ERR1779327","ERR1779328","ERR1779335","ERR1779337","ERR1779341","ERR1779343","ERR1779356","ERR1779358","ERR1779360","ERR1779362","ERR1779365","ERR1779367","ERR1779397","ERR1779404","ERR1779405","ERR1779436","ERR1779437","ERR1779439","ERR1779440","ERR1779447","ERR1779449","ERR1779453","ERR1779455","ERR1779459","ERR1779460","ERR1779463","ERR1779467","ERR1779472","ERR1779473","ERR1779475","ERR1779477","ERR1779479","ERR1779482","ERR1779484","ERR1779521","ERR1779522","SRR8512352","SRR8512354","SRR8512356","SRR8512358","SRR8512360","SRR8512362","SRR8512364","SRR8512366","SRR8512368","SRR8512370","SRR8512372","SRR8512374","SRR8512414","SRR8512417","SRR8512419","SRR8512423","SRR8512425","SRR8512427","SRR8512429","SRR8512431","SRR8512433")

count1 <- count[ , !(names(count) %in% removelist)]
write.csv(count1, "EPmerge_remove.csv", row.names = FALSE)
countAD <- read.csv("ADmerge.csv")
countAD1 <- countAD[ , !(names(countAD) %in% removelist)]
write.csv(countAD1, "ADmerge_remove.csv", row.names = FALSE)
