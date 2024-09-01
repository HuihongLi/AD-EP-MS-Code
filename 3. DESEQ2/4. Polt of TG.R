library(ggplot2)
library(tidyverse)
library(ggrepel)
sigresultsEP <- read.csv("sigresultsEP.csv", header = TRUE, sep = ",")

sigresultsEP <- sigresultsEP[,2]
sigresultsTG <- read.csv("Tg4510_sig_interaction-ONLY_STATS.csv", header = TRUE, sep = ",")

sigresultsTG_filter <- sigresultsTG[sigresultsTG$GenotyperTg4510.Age_months8 > 1 | sigresultsTG$GenotyperTg4510.Age_months8 < -1,]


sigresultsTG_filter_EP_DEG <- sigresultsTG_filter[sigresultsTG_filter$X %in% sigresultsEP,]
write.csv(sigresultsTG_filter_EP_DEG, "sigresultsTG_filter_EP_DEG.csv", row.names = FALSE)



df1 <- sigresultsTG_filter[,c(1,2,3)]
df1$cluster <- factor("4")
df2 <- sigresultsTG_filter[,c(1,2,4)]
df2$cluster <- factor("6")
df3 <- sigresultsTG_filter[,c(1,2,5)]
df3$cluster <- factor("8")

colnames(df1)[3] <- "log2FoldChange"
colnames(df2)[3] <- "log2FoldChange"
colnames(df3)[3] <- "log2FoldChange"
df <- rbind(df1,df2,df3)



df$Label <- ifelse(df$X %in% sigresultsEP, "EP-DEG", "NotEP-DEG")
df1$Label <- ifelse(df1$X %in% sigresultsEP, "EP-DEG", "NotEP-DEG")
df2$Label <- ifelse(df2$X %in% sigresultsEP, "EP-DEG", "NotEP-DEG")
df3$Label <- ifelse(df3$X %in% sigresultsEP, "EP-DEG", "NotEP-DEG")


df1 %>% filter(Label == "EP-DEG" & abs(log2FoldChange)>1) %>% summarise(n = n())

df1 %>% filter(abs(log2FoldChange)>1) %>% summarise(n = n())

df2 %>% filter(Label == "EP-DEG" & abs(log2FoldChange)>1) %>% summarise(n = n())
df2 %>% filter(abs(log2FoldChange)>1) %>% summarise(n = n())

df3 %>% filter(Label == "EP-DEG" & abs(log2FoldChange)>1) %>% summarise(n = n())
df3 %>% filter(abs(log2FoldChange)>1) %>% summarise(n = n())


df %>% filter(Label == "EP-DEG") %>% group_by(cluster) %>% summarise(n = n())
p <- ggplot() +
  geom_jitter(data = df,
              aes(x = cluster, y = log2FoldChange, color = Label, size = -log2(FDR_adj_LRT)),
              width = 0.4) 
p


dfbar<-data.frame(x=c('4','6','8'),
                  y=c(6.2,7,7))
dfbar1<-data.frame(x=c('4','6','8'),
                   y=c(-2.1,-3.5,-2.1))


p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1


p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = df,
              aes(x = cluster, y = log2FoldChange, color = Label,size = -log2(FDR_adj_LRT)),
              width =0.4)+scale_size_continuous(name = "-Log2(Padj)")
p2


dfcol<-data.frame(x=c('4','6','8'),
                  y=0,
                  label=c('TG_4\n(60/151)','TG_6\n(83/268)','TG_8\n(101/426)'))
mycol <- c("#4DBBD57F","#00A0877F","#3C54887F")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=1.5,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3



p6 <- p3+
  labs(x="Age (Month)",y="log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")
p6

p7 <- p6+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 15)
  )
p7
pdf("MultiDETG.pdf",width = 10,height = 8)
p7
dev.off()



