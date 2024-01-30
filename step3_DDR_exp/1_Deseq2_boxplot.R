rm(list=ls())
library(DESeq2)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(openxlsx)
library(ggsignif)
library(ggpubr)
#导入数据
load(file = "./meddata/mcounts_GSE114007.Rdata")
load(file="./meddata/step2_DE_GSE114007.Rdata")
group_info_GSE114007
#样品信息矩阵
coldata=as.data.frame(colnames(m_counts_GSE114007))
colnames(coldata)=c("tittle")

coldata$condition[str_which(coldata$tittle,'OA')]=c("OA")
coldata$condition[str_which(coldata$tittle,'normal')]=c("NORMAL")
coldata$condition[str_which(coldata$tittle,'Normal')]=c("NORMAL")
coldata$condition=factor(coldata$condition,levels=c("NORMAL","OA"))
rownames(coldata)=coldata$tittle
dim(m_counts_GSE114007)
all(colnames(m_counts_GSE114007) == rownames(coldata))

#normalized read counts
dds=DESeqDataSetFromMatrix(countData = m_counts_GSE114007,
                           colData = coldata,
                           design = ~condition)
dds=estimateSizeFactors(dds)
sizeFactors(dds)
plot(sizeFactors(dds),colSums(counts(dds)))
abline(lm(colSums(counts(dds))~sizeFactors(dds)+0))
normlzd_dds=counts(dds, normalized=T)
head(normlzd_dds)
head(m_counts_GSE114007)
normlzd_dds=as.data.frame(normlzd_dds)
normlzd_dds$Genename=rownames(normlzd_dds)
save(normlzd_dds,coldata,file="Normalized_GSE114007.R")




rm(list=ls())
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(openxlsx)
library(ggsignif)
library(ggpubr)

DDR=read.xlsx("DDR_gene.xlsx")
load("Normalized_GSE114007.R")
#normlzd_dds=read.csv("Normalized_GSE114007.csv",header = T)
#genenames=intersect(normlzd_dds$Genename,DDR$Genename) #188

#表达框筛选DDR基因
df=left_join(DDR,normlzd_dds,by="Genename") #233(有重复基因)
#type=c("BER","FA","NER","MMR","CPF","HRR","TLS","NHEJ")
#df$type=factor(df$type,levels = "BER","FA","NER","MMR","CPF","HRR","TLS","NHEJ")
df=pivot_longer(df,
                cols = Normal_Cart_10_8:OA_10,
                names_to = "sample",
                values_to = "normalied_value")
colnames(coldata)=c("sample","condition")
df=left_join(df,coldata,by="sample")


dim(df)
df$condition=factor(df$condition,levels = unique(df$condition),ordered = T)
df$type=factor(df$type,levels = unique(df$type),ordered = T)
df$normalied_value=as.numeric(df$normalied_value)

#画图
b_p=ggplot(df,aes(x=type,y=normalied_value,fill=condition))+
  geom_boxplot(width=0.6,alpha=0.8)+
  geom_signif(comparisons= list(df$condition),map_signif_level = TRUE)+
  theme_bw()

b_p=b_p+stat_compare_means(method="t.test",
                           label = "p.signif")
b_p
