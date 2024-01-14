### 需要注意的点：
# 1.得到的矩阵是counts、TPM、FPKM中的哪一类，需要统一成TMP方便样本间比较！信息从文件名\文章\GEO页面去找。
# 如果没做特殊说明的话，这里的数值一般表示基因数量，即counts，那么下载到的表达量就是log2(counts+1)。
# 2.探针编号，GLP开头，探针与基因名对应，不同批次间的GENE数量可能会不一样，因为有些探针会有多对一的情况，也需要注意标准化。
# 3.single cell样本可能需要统一地重新聚类一下，然后注释，工作量很大。先暂时用文章里的分类，所以一定要找到cluster信息。
# 4.中间文件都放在"./meddata/mcounts_GSE****.Rdata"这样的一个Rdata下面，便于统一处理
rm(list = ls())
library(GEOquery)
library(dplyr)
#这一步可以得到表达矩阵，但经常不成功，需要自己去重新下载
downGSE <- function(studyID = "GSE168505", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)

  # 保存
  # write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  # write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
}
eSet=downGSE("GSE168505")
# 表达矩阵，经常为空
exp = exprs(eSet[[1]])
# 临床信息
pdata_GSE168505 = pData(eSet[[1]])

#数据另下
#记录了表达矩阵的文件
# 这里文件里提到了RSEM v1.3.1 (rsem-calculate-expression function w/ STAR aligning)，通常就是FPKM矩阵了（看文章以进一步确定！）
#文章里用的DESeq2，这里应该是校正后的counts数量，说的很不清楚
matrixfile=list.files("./GSE168505/",full.names = T,pattern = ".tsv")
matrixfile
m_fpkm_GSE168505=read.table(matrixfile[1],header = T)
row.names(m_fpkm_GSE168505)=m_fpkm_GSE168505$gene
m_fpkm_GSE168505=m_fpkm_GSE168505[,-1]
#文章中的分组信息
group_info_GSE168505= pdata_GSE168505$`disease state:ch1`
group_info_GSE168505
#保存表达矩阵
save(m_fpkm_GSE168505,group_info_GSE168505,pdata_GSE168505,file = "./meddata/mcounts_GSE168505.Rdata")



