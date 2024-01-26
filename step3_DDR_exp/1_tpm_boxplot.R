rm(list = ls()) 

###加载数据(begin)
load("./meddata/ID_gene_length.Rdata")
load("./meddata/mcounts_GSE114007.Rdata")
library(tidyverse)
library(stringr)
#基因表达量数据框中基因对应长度

gene_GSE114007=rownames(m_counts_GSE114007)
gene_GSE114007=as.data.frame(gene_GSE114007)
colnames(gene_GSE114007)=c("gene_name") #23710
gene_GSE114007=merge(id_GN_length,gene_GSE114007,by.y="gene_name") #20920
gene_GSE114007 #20920
gene_GSE114007$length=as.numeric(gene_GSE114007$length)
m_counts_GSE114007$gene_name=rownames(m_counts_GSE114007) #23710 gene
#数据转换
tran_df=merge(m_counts_GSE114007,id_GN_length,by.x ="gene_name" ) #20920 gene
tran_df=tran_df[!duplicated(tran_df$gene_name),] #20672
head(tran_df)
rownames(tran_df)=tran_df[,1]
tran_df=tran_df[,-1]

kb=tran_df$length /1000
head(kb)

#计算TPM
countdata= tran_df[,1:38]
head(countdata)
rpk=countdata / kb
tpm =t(t(rpk)/colSums(rpk) *100000)
save(tpm,file="m_tpm2_GSE114007.Rdata")

rm(list = ls())
load("./m_tpm2_GSE114007.Rdata")

tpm=as.data.frame(tpm)
library(readr)
library(patchwork)
library(ggpubr)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tidyr)
#对TPM值进行标化
tpm_nmlz=apply(tpm,2,function(x)return(x/mean(x)))
tpm_nmlz=log10(tpm_nmlz)
#boxplot(log10(tpm_nmlz))

m_tpm_GSE114007=tpm_nmlz
DDR=read.xlsx("DDR_gene.xlsx")
load("./meddata/mcounts_GSE114007.Rdata")
#表达框筛选DDR基因
genename=intersect(rownames(m_tpm_GSE114007),DDR$Genename) #188
DDR_tpm_GSE114007=m_tpm_GSE114007[genename, ]

#DDR_GSE114007=DDR[genename %in% DDR$Genename,]
DDR_tpm_GSE114007=DDR_tpm_GSE114007[which(rowSums(DDR_tpm_GSE114007)>0),] #186#去掉在所有样本中表达量为0的基因
genename=DDR_tpm_GSE114007
DDR_tpm_all_GSE114007=t(DDR_tpm_GSE114007)
DDR_tpm_all_GSE114007=as.data.frame(DDR_tpm_all_GSE114007)
DDR_tpm_all_GSE114007$sample=rownames(DDR_tpm_all_GSE114007)
DDR_tpm_all_GSE114007$sap_grp=group_info_GSE114007
cols=rownames(DDR_tpm_GSE114007)

df2=DDR_tpm_GSE114007
df2=as.data.frame(df2)
df2$Genename=rownames(df2)
df2=left_join(df2,DDR,by="Genename") 
sum(is.na(df2))
#df2=df2[,-39]
#all(colnames(DDR_tpm_GSE114007) == DDR_1$Genename)
df2=pivot_longer(df2,
                 cols = Normal_Cart_10_8:OA_10,
                 names_to = "sample",
                 values_to = "TPM")
load("./meddata/mcounts_GSE114007.Rdata")
#df2$sample=factor(df2$sample,levels = df2$sample,ordered = T)
df2=as.data.frame(df2)

df2= df2 %>% arrange(sample)
df2$group=c(rep("NORMAL",684),rep("OA",684)) #数字按
df2$group=factor(df2$group,levels = c("NORMAL","OA"))
#print(df2$sample)
##画图
library(ggsignif)
dim(df2)
df2$group=factor(df2$group,levels = unique(df2$group),ordered = T)
df2$type=factor(df2$type,levels = unique(df2$type),ordered = T)
df2$TPM=as.numeric(df2$TPM)

#df2$lable

b_p_1=ggplot(df2,aes(x=type,y=TPM,fill=group))+
  geom_boxplot(width=0.6,alpha=0.8)+
  #geom_signif(comparisons= list(df2$group),map_signif_level = TRUE)+
  theme_bw()

b_p_1=b_p_1+stat_compare_means(method="t.test",
                               label="p.format")
b_p_1


#b_p=ggplot(df2,aes(x=group,y=TPM,fill=type))+
  #geom_boxplot(width=0.6,alpha=0.8)+
  #geom_signif(comparisons= list(c("BER","FA","NER","MMR","CPF","HRR","TLS","NHEJ")),map_signif_level = TRUE)+
  #theme_bw()
#b_p=b_p+stat_compare_means(method="wilcox.test",
                            #   label="p.format")
#b_p
