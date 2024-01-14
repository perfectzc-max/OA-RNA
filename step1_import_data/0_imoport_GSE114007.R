### 需要注意的点：
# 1.得到的矩阵是counts、TPM、FPKM中的哪一类，需要统一成TPM方便样本间比较！信息从文件名\文章\GEO页面去找。
# 如果没做特殊说明的话，这里的数值一般表示基因数量，即counts，那么下载到的表达量就是log2(counts+1)。
# 2.探针编号，GLP开头，探针与基因名对应，不同批次间的GENE数量可能会不一样，因为有些探针会有多对一的情况，也需要注意标准化。
# 3.single cell样本可能需要统一地重新聚类一下，然后注释，工作量很大。先暂时用文章里的分类，所以一定要找到cluster信息。
# 4.中间文件都放在"./meddata/mcounts_GSE****.Rdata"这样的一个Rdata下面，便于统一处理
rm(list = ls())
library(GEOquery)
library(dplyr)
library(readxl)
library(tidyverse)
#这一步可以得到表达矩阵，但经常不成功，需要自己去重新下载
downGSE <- function(studyID = "GSE114007", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)

  # 保存
  # write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  # write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
}
eSet=downGSE("GSE114007")
# 表达矩阵，经常为空
exp = exprs(eSet[[1]])
# 临床信息
pdata_GSE114007 = rbind(pData(eSet[[1]]),pData(eSet[[2]]))
pdata_GSE114007$title = sort(pdata_GSE114007$title)
pdata_GSE114007
#数据另下
#记录了表达矩阵的文件
# 这里文件里提到了counts，那就是原始counts矩阵
matrixfile=list.files("./",full.names = T,pattern = "GSE114007_raw")
m_counts_GSE114007_1=read_excel("GSE114007_raw_counts.xlsx") %>% as.data.frame()
#表格有两个sheet
m_counts_GSE114007_2=read_excel("GSE114007_raw_counts.xlsx", sheet = "OA") %>% as.data.frame()
m_counts_GSE114007=full_join(m_counts_GSE114007_1,m_counts_GSE114007_2,by="symbol")
row.names(m_counts_GSE114007_1)=as.character(m_counts_GSE114007_1$symbol)
row.names(m_counts_GSE114007_2)=as.character(m_counts_GSE114007_2$symbol)
m_counts_GSE114007_1=m_counts_GSE114007_1[,-1]
m_counts_GSE114007_2=m_counts_GSE114007_2[,-1]
rownames(m_counts_GSE114007)=m_counts_GSE114007$symbol
m_counts_GSE114007=m_counts_GSE114007[,-1]
#文章中的分组信息
group_info_GSE114007=str_split(pdata_GSE114007$title,"_",simplify = T)[,1]
group_info_GSE114007=str_to_upper(group_info_GSE114007)
group_info_GSE114007
#保存表达矩阵
save(m_counts_GSE114007,group_info_GSE114007,pdata_GSE114007,file = "./meddata/mcounts_GSE114007.Rdata")
