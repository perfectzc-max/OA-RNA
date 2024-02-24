rm(list=ls())
#####https://blog.csdn.net/weixin_43843918/article/details/135899489
#if(!requireNamespace("BiocManager",quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("edgeR")
library(DESeq2)
library(edgeR)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(openxlsx)
library(ggsignif)
library(ggpubr)
#BiocManager::install("statmod")
#BiocManager::install("limma")

#install.packages("edgeR")
library(limma)

#1.导入数据
load(file = "./meddata/mcounts_GSE114007.Rdata") #23710
load(file="./meddata/mcounts_GSE51588.R")
#group_info_GSE114007

#2.使用limma对row counts进行归一化
exprSet= m_counts_GSE114007
table(group_info_GSE114007)
design=model.matrix(~0+factor(group_info_GSE114007))

colnames(design)=c("NORMAL","OA")
group_info_GSE114007=factor(group_info_GSE114007,levels=c("NORMAL","OA"))
rownames(design)=colnames(exprSet)

dge=DGEList(counts=exprSet,group=group_info_GSE114007)
keep=filterByExpr(dge)
dge=dge[keep, , keep.lib.sizes=FALSE]
dge=calcNormFactors(dge)
v=voom(dge,design,plot = TRUE,normalize="quantile")
m_counts_GSE114007_normlzd=as.data.frame(v[["E"]]) #13629


#3.合并样品信息矩阵
coldata_GSE114007=as.data.frame(colnames(m_counts_GSE114007_normlzd))
coldata_GSE114007$GSE=rep("GSE114007",38)
colnames(coldata_GSE114007)=c("tittle","GSE")
coldata_GSE51588=as.data.frame(group_info_GSE51588$title)
coldata_GSE51588$GSE=rep("GSE51588",50)
colnames(coldata_GSE51588)=c("tittle","GSE")
coldata=rbind(coldata_GSE114007,coldata_GSE51588)
#colnames(coldata)=c("tittle") 

coldata$condition[str_which(coldata$tittle,'OA')]=c("OA")
coldata$condition[str_which(coldata$tittle,'normal')]=c("NORMAL")
coldata$condition[str_which(coldata$tittle,'Normal')]=c("NORMAL")
coldata$condition=factor(coldata$condition,levels=c("NORMAL","OA"))
rownames(coldata)=coldata$tittle

library(dplyr)
#dim(m_counts_GSE114007)
#all(colnames(m_counts_GSE114007) == rownames(coldata))


#4.合并表达数据框
m_counts_GSE114007_normlzd$genename=rownames(m_counts_GSE114007_normlzd)
#diff_gene=setdiff(m_counts_GSE114007$genename,m_counts_GSE51588$genename)
m_counts=merge(m_counts_GSE114007_normlzd,m_counts_GSE51588,by="genename")#before 27328 after 18658
m_counts=m_counts[!duplicated(m_counts$genename),] #before 18658 after 12733
rownames(m_counts)=m_counts$genename
m_counts=m_counts[,-1]

m_counts_1=full_join(m_counts_GSE114007_normlzd,m_counts_GSE51588) #before 33901 after 30729
table(is.na(m_counts_1))
m_counts_1=na.omit(m_counts_1)#18658

save(coldata,m_counts,m_counts_1,file = "./meddata/step3_bind_matrix.Rdata")


rm(list = ls())
load("./meddata/step3_bind_matrix.Rdata")
#基因放到第一列
m_counts_1=select(m_counts_1,39,everything())
#样本
m_counts_2=data.frame(t(m_counts_1))
colnames(m_counts_2)=m_counts_2[1,]
m_counts_2=m_counts_2[-1,]
sample=rownames(m_counts_2)

#数据框内变量类型由字符变数字
m_counts_2=as.data.frame(lapply(m_counts_2,as.numeric))
rownames(m_counts_2)=sample
class(m_counts_2$A2M)

library(tinyarray)
coldata$GSE=factor(coldata$GSE,levels = c("GSE114007","GSE51588"))
PCA_before=prcomp(m_counts_2,
                  scale = T)
library(ggplot2)
#install_github("vqv/ggbiplot")
library(ggbiplot)
coldata_2=coldata[,1:2]
#before removing batch effect
m_pca_b=m_counts_1
m_pca_b=m_pca_b[!duplicated(m_pca_b$genename),] #去掉重复基因 12733
rownames(m_pca_b)=m_pca_b$genename
m_pca_b=m_pca_b[,-1]
m_pca_b[1:4,1:4]
boxplot(m_pca_b,las=2, cex.axis=0.6)
#normalization
m_pca_b=scale(m_pca_b)
boxplot(m_pca_b,las=2,cex.axis=0.6)
#pca visulization
draw_pca(m_pca_b,group_list = coldata_2$GSE)
#remove batch effect(sva Combat package)
##BiocManager::install("sva")
library(sva)
batch=coldata$GSE
m_pca_a=ComBat(m_pca_b,batch=batch)
#after removing batch effect
boxplot(m_pca_a, las=2,cex.axis=0.6)
draw_pca(exp=m_pca_a,group_list = coldata_2$GSE)

library("openxlsx")
DDR=read.xlsx("DDR_gene.xlsx")
m_counts=as.data.frame(data.matrix(m_pca_a)) #12733
m_counts$Genename=rownames(m_counts)

#normalized read counts
#dds=DESeqDataSetFromMatrix(countData = m_counts,
                           #colData = coldata,
                           #design = ~condition)
#dds=estimateSizeFactors(dds)
#sizeFactors(dds)
#plot(sizeFactors(dds),colSums(counts(dds)))
#abline(lm(colSums(counts(dds))~sizeFactors(dds)+0))
#normlzd_dds=counts(dds, normalized=T)
#head(normlzd_dds)
#head(m_counts_GSE114007)
#normlzd_dds=as.data.frame(normlzd_dds)
#normlzd_dds$Genename=rownames(normlzd_dds)
#save(normlzd_dds,coldata,file="Normalized_GSE114007.R")

#normlzd_dds=read.csv("Normalized_GSE114007.csv",header = T)
#genenames=intersect(normlzd_dds$Genename,DDR$Genename) #188

#5.表达框筛选DDR基因
m_counts$Genename=rownames(m_counts)
df=left_join(DDR,m_counts,by="Genename") #233(有重复基因)
type=c("BER","FA","NER","MMR","CPF","HRR","TLS","NHEJ")
df$type=factor(df$type,levels = type)
col=colnames(df)
col=gsub("-","_",col)
colnames(df)=col
colnames(coldata)=c("sample","GSE","condition")
#df=left_join(df,coldata,by="sample")
#m_counts=m_counts[,c(1:88)]
#m_counts_1=as.matrix(m_counts_1)



#boxplot(m_counts_1)
#画图
df=pivot_longer(df,
                cols = Normal_Cart_10_8 : OA_MT_9,
                names_to = "sample",
                values_to = "normalized_value")



dim(df)#20504

df$condition[str_which(df$sample,'OA')]=c("OA")
df$condition[str_which(df$sample,'ormal')]=c("NORMAL")
#df$condition[str_which(df$tittle,'Normal')]=c("NORMAL")

df$condition=factor(df$condition,levels = unique(df$condition),ordered = T)
df$type=factor(df$type,levels = unique(df$type),ordered = T)
df$normalized_value=as.numeric(df$normalized_value)

#画图
b_p=ggplot(df,aes(x=type,y=normalized_value,fill=condition))+
  geom_boxplot(width=0.6,alpha=0.8)+
  geom_signif(comparisons= list(df$condition),map_signif_level = TRUE)+
  theme_bw()

b_p=b_p+stat_compare_means(method="t.test",
                           label = "p.signif")
b_p

