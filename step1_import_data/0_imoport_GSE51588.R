rm(list = ls())
library(GEOquery)
library(dplyr)
library(readxl)
library(tidyverse)
#这一步可以得到表达矩阵，但经常不成功，需要自己去重新下载
downGSE <- function(studyID = "GSE51588", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
  
  #保存
  #write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  #write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
}
eSet=downGSE("GSE51588")
# 表达矩阵，经常为空
exp = exprs(eSet[[1]])

#matrixfile=list.files("./GSE51588_RAW/",pattern = ".txt",full.names = T)
matrixfile

m_counts_GSE51588=as.data.frame(exp)


# 临床信息
pdata_GSE51588 = pData(eSet[[1]])
pdata_GSE51588$title = sort(pdata_GSE51588$title)
#pdata_GSE51588
#文章中的分组信息
group_info_GSE51588= pdata_GSE51588$characteristics_ch1
group_info_GSE51588=as.data.frame(group_info_GSE51588)
group_info_GSE51588$title=pdata_GSE51588$title
group_info_GSE51588$geo_accession=pdata_GSE51588$geo_accession

#save(group_info_GSE51588,pdata_GSE51588,group_info_GSE51588,file = "./meddata/mcounts_GSE51588.R")
#探针注释
GPL13497=getGEO('GPL13497',destdir = ".")
GPL13497_anno=Table(GPL13497)
ids = GPL13497_anno[,c(1,7)]
colnames(ids) = c("probe_id","symbol")
ids = ids[ids$symbol!="" & !str_detect(ids$symbol,"///"),]

m_counts_GSE51588$probe_id=rownames(m_counts_GSE51588)
m_counts_GSE51588=left_join(ids,m_counts_GSE51588,by = "probe_id")
m_counts_GSE51588=m_counts_GSE51588[,-1]
title=group_info_GSE51588$title
colnames(m_counts_GSE51588)=c("genename",title)
save(m_counts_GSE51588,pdata_GSE51588,group_info_GSE51588,file = "./meddata/mcounts_GSE51588.R")
#m_counts_GSE51588(29833*51) pdata_GSE51588(50*42)  group_info_GSE51588(50*3)
rm(list = ls())
load("./meddata/mcounts_GSE51588.Rdata")
