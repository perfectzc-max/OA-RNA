##GSE114007 + GSE51588
rm(list = ls())
options(stringsAsFactors=F)

library(tidyverse)
library(clusterProfiler)
#install.packages("msigdbr")
library(msigdbr)
#BiocManager::install("GSVA")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
library(BiocParallel)
library(edgeR)
#1.data prepare
load("./meddata/mcounts_GSE114007.Rdata")
load("./meddata/mcounts_GSE51588.R")
#load("./meddata/step3_bind_matrix.Rdata")

##normalize GSE114007(row counts)
exprSet= m_counts_GSE114007
table(group_info_GSE114007)
design=model.matrix(~0+factor(group_info_GSE114007))

colnames(design)=c("NORMAL","OA")
group_info_GSE114007=factor(group_info_GSE114007,levels=c("NORMAL","OA"))
rownames(design)=colnames(exprSet)

dge=DGEList(counts=exprSet,group=group_info_GSE114007)
#keep=filterByExpr(dge)
#dge=dge[keep, , keep.lib.sizes=FALSE]
dge=calcNormFactors(dge)
v=voom(dge,design,plot = TRUE,normalize="quantile")
m_normal_GSE114007=as.data.frame(v[["E"]]) 
m_normal_GSE114007$genename=rownames(m_normal_GSE114007)

##bind data frame
m_counts=full_join(m_normal_GSE114007,m_counts_GSE51588,by = "genename") 
m_counts[is.na(m_counts)]=0.001
m_counts=na.omit(m_counts) #33901
m_counts=dplyr::select(m_counts,39,everything()) 
#rownames(m_counts)=m_counts$genename
##remove batch effect
#m_counts=as.data.frame(lapply(m_counts,as.numeric))
class(m_counts)
m_counts_2=data.frame(t(m_counts))
colnames(m_counts_2)=m_counts_2[1,]
m_counts_2=m_counts_2[-1,]
sample=rownames(m_counts_2)

m_counts_2=as.data.frame(lapply(m_counts_2,as.numeric))
rownames(m_counts_2)=sample
class(m_counts_2$APEX1)

library(tinyarray)
PCA_before=prcomp(m_counts_2,
                  scale = T)
PCA_before

m_pca_b=m_counts
m_pca_b=m_pca_b[!duplicated(m_pca_b$genename),]#25822
rownames(m_pca_b)=m_pca_b$genename
m_pca_b=m_pca_b[,-1]
m_pca_b[1:4,1:4]
class(m_pca_b)
#m_pca_b=as.factor(m_pca_b)
boxplot(m_pca_b,las=2, cex.axis=0.6)
#normalization
m_pca_b=scale(m_pca_b)
boxplot(m_pca_b,las=2,cex.axis=0.6)
#pca visulization
load("./meddata/step3_bindDDR_matrix.Rdata")
coldata_2=coldata[,1:2]
coldata_2$GSE=factor(coldata_2$GSE,levels = c("GSE114007","GSE51588"))
draw_pca(m_pca_b,group_list = coldata_2$GSE)
#remove batch effect(sva Combat package)
##BiocManager::install("sva")
library(sva)
batch=coldata$GSE
m_pca_a=ComBat(m_pca_b,batch=batch)
#after removing batch effect
boxplot(m_pca_a, las=2,cex.axis=0.6)
draw_pca(exp=m_pca_a,group_list = coldata_2$GSE)
write.csv(m_pca_a,"./meddata/bind_data.csv")
