### 需要注意的点：
# 1.得到的矩阵是counts、TPM、FPKM中的哪一类，需要统一成TPM方便样本间比较！信息从文件名\文章\GEO页面去找。
# 2.探针编号，GLP开头，探针与基因名对应，不同批次间的GENE数量可能会不一样，因为有些探针会有多对一的情况，也需要注意标准化。
# 3.single cell样本可能需要统一地重新聚类一下，然后注释，工作量很大。先暂时用文章里的分类，所以一定要找到cluster信息。
# 4.中间文件都放在"./meddata/mcounts_GSE****.Rdata"这样的一个Rdata下面，便于统一处理
rm(list = ls())
library(GEOquery)
library(dplyr)
downGSE <- function(studyID = "GSE104782", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
  
  # exprSet = exprs(eSet[[1]])
  # pdata = pData(eSet[[1]])
  # 
  # write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  # write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
  
}
eSet=downGSE("GSE104782")
#表达矩阵（此处为空）
exp= exprs(eSet[[1]])

#数据另下
#这里的表达矩阵是TPM矩阵
matrixfile=list.files("./GSE104782_RAW/",pattern = ".txt",full.names = T)

m_tpm=read.table(matrixfile[1],header = T)

# left_join(temp1,temp2,by="Gene")
for (i in 2:length(matrixfile)) {
  temp=read.table(matrixfile[i],header = T)
  m_tpm=left_join(temp,m_tpm,by="Gene")

}
#查看数据分布
boxplot(m_tpm[,1])

#临床测序信息
pdata_GSE104782 = pData(eSet[[1]])
#文章中的cluster信息
library(readxl)
group_cluster_info_GSE104782 <- read_excel("meddata/GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx")
# View(GSE104782_clustering_information)
# save(GSE104782_clustering_information,file = "./meddata/cluster.Rdata")


#样本信息
# GSE104782_clustering_information$Cluster
#筛掉不合格样本,1464个细胞纳入分析
row.names(m_tpm)=m_tpm$Gene
m_tpm=m_tpm[,-1]
m_tpm_GSE104782=m_tpm[,!is.na(group_cluster_info_GSE104782$Cluster)]
save(m_tpm_GSE104782,group_cluster_info_GSE104782,pdata_GSE104782,file = "./meddata/mcounts_GSE104782.Rdata")

#如果行名不是GENE
###提取芯片平台编号
gpl_number <- eSet[[1]]@annotation

#探针注释的获取
#方法1 BioconductorR包(最常用)
gpl_number 
#http://www.bio-info-trainee.com/1399.html
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
#转换为探针id与基因名称对应表格
ids <- toTable(hgu133plus2SYMBOL)
head(ids)


# 方法2 读取GPL平台的soft文件，按列取子集
##https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
if(F){
  #注：soft文件列名不统一，活学活用，有的GPL平台没有提供注释，如GPL16956
  a = getGEO(gpl_number,destdir = ".")
  b = a@dataTable@table
  colnames(b)
  ids2 = b[,c("ID","Gene Symbol")]
  colnames(ids2) = c("probe_id","symbol")
  ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]
}

# 方法3 官网下载,文件读取
##http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus

# 方法4 自主注释 
#https://mp.weixin.qq.com/s/mrtjpN8yDKUdCSvSUuUwcA
save(exp,Group,ids,gse_number,file = "step2output.Rdata")