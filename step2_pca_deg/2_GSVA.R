##GSVA
####https://zhuanlan.zhihu.com/p/518145829

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

#2.gene list(msigdbr package)
##KEGG
KEGG_df_all = msigdbr(species = "Homo sapiens",
                      category = "C2",
                      subcategory = "CP:KEGG")
KEGG_df = dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list = split(KEGG_df$gene_symbol,KEGG_df$gs_name)#group gene_symbol by gs_name

##GO
GO_df_all = msigdbr(species = "Homo sapiens",
                    category = "C5")
GO_df = dplyr::select(GO_df_all,gs_name,gene_symbol,gs_exact_source,gs_subcat)
GO_df = GO_df[GO_df$gs_subcat!="HPO",]
go_list = split(GO_df$gene_symbol,GO_df$gs_name)#group gene_symbol by gs_name

#3.GSVA-go
#dat=t(m_counts)
dat=as.matrix(m_pca_a)
geneset= go_list
gsva_mat = gsva(expr=dat,
                gset.idx.list = geneset,
                kcdf="Gaussian", #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                verbose=T,
                parallel.sz=parallel::detectCores())
gsva_bind_go_mat=gsva_mat
save(dat,design,gsva_bind_go_mat,coldata,file = "./meddata/gsva_114007+51588_go_matrix.R")
#################################################################################################


#4.GSVA-kegg
#dat=t(m_counts)
dat=as.matrix(m_pca_a)
#geneset= kegg_list
gsva_bind_kegg_mat = gsva(expr=dat,
                          gset.idx.list = kegg_list,
                          kcdf="Gaussian", #"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                          verbose=T,
                          parallel.sz=parallel::detectCores())

save(dat,gsva_bind_kegg_mat,coldata,file = "./meddata/gsva_114007+51588_kegg_matrix.R")
#################################################################################################