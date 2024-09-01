EPmerge_remove_afterPCA <- read.csv("EPmerge_remove_afterPCA.csv",row.names = 1)
ADmerge_remove <- read.csv("ADmerge_remove.csv",row.names = 1)
metaEP <- read.csv("metaEP.csv",row.names = 1)
metaTG <- read.csv("metaTG.csv",row.names = 1)
metaJ20<- read.csv("metaJ20.csv",row.names = 1)

EPmerge_remove_afterPCA_sum <- data.frame(colSums(EPmerge_remove_afterPCA))
ADmerge_remove_sum <- data.frame(colSums(ADmerge_remove))
EPmerge_remove_afterPCA_sum_mean <- mean(EPmerge_remove_afterPCA_sum[,1])
EPmerge_remove_afterPCA_sum_sd <- sd(EPmerge_remove_afterPCA_sum[,1])
ADmerge_remove_sum_mean <- mean(ADmerge_remove_sum[,1])
ADmerge_remove_sum_sd <- sd(ADmerge_remove_sum[,1])


EPmerge_remove_afterPCA_sum_EP <-EPmerge_remove_afterPCA_sum[metaEP$Condition == "EP",]
EPmerge_remove_afterPCA_sum_EP_WT <-EPmerge_remove_afterPCA_sum[metaEP$Condition == "WT",]
J20merge_remove_sum_J20 <-ADmerge_remove_sum[rownames(metaJ20)[metaJ20$Genotype == "J20"],]
J20merge_remove_sum_WT <-ADmerge_remove_sum[rownames(metaJ20)[metaJ20$Genotype == "WT"],]
TGmerge_remove_sum_TG <-ADmerge_remove_sum[rownames(metaTG)[metaTG$Genotype == "rTg4510"],]
TGmerge_remove_sum_WT <-ADmerge_remove_sum[rownames(metaTG)[metaTG$Genotype == "WT"],]


EPmerge_remove_afterPCA_sum_EP_mean <- mean(EPmerge_remove_afterPCA_sum_EP)
EPmerge_remove_afterPCA_sum_EP_sd <- sd(EPmerge_remove_afterPCA_sum_EP)
EPmerge_remove_afterPCA_sum_EP_WT_mean <- mean(EPmerge_remove_afterPCA_sum_EP_WT)
EPmerge_remove_afterPCA_sum_EP_WT_sd <- sd(EPmerge_remove_afterPCA_sum_EP_WT)
J20merge_remove_sum_J20_mean <- mean(J20merge_remove_sum_J20)
J20merge_remove_sum_J20_sd <- sd(J20merge_remove_sum_J20)
J20merge_remove_sum_WT_mean <- mean(J20merge_remove_sum_WT)
J20merge_remove_sum_WT_sd <- sd(J20merge_remove_sum_WT)
TGmerge_remove_sum_TG_mean <- mean(TGmerge_remove_sum_TG)
TGmerge_remove_sum_TG_sd <- sd(TGmerge_remove_sum_TG)
TGmerge_remove_sum_WT_mean <- mean(TGmerge_remove_sum_WT)
TGmerge_remove_sum_WT_sd <- sd(TGmerge_remove_sum_WT)


t_test_result_EP <- t.test(EPmerge_remove_afterPCA_sum_EP, EPmerge_remove_afterPCA_sum_EP_WT, 
                        alternative = "two.sided", var.equal = FALSE)
t_test_result_J20 <- t.test(J20merge_remove_sum_J20, J20merge_remove_sum_WT, 
                        alternative = "two.sided", var.equal = FALSE)
t_test_result_TG <- t.test(TGmerge_remove_sum_TG, TGmerge_remove_sum_WT,
                        alternative = "two.sided", var.equal = FALSE)

t_test_result_EP
t_test_result_TG
t_test_result_J20


#返回metaJ20$Genotype == "J20"的行名
