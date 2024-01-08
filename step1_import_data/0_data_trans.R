#####参考https://blog.csdn.net/weixin_56259962/article/details/123140630

#从gtf文件调取基因
library(GenomicFeatures)
##将数据导入为txdb对象（等待时间较长）
txdb=makeTxDbFromGFF("Homo_sapiens.GRCh38.110.gtf",format = "gtf")
exons_gene=exonsBy(txdb,by="gene")
##计算总外显子长度（等待时间长）
exons_gene_lens=lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_lens[1:10]
##转换为dataframe
gene_length=sapply(exons_gene_lens,function(x){x})
id_length=as.data.frame(gene_length)
id_length$GenestableID=rownames(id_length)
#id_length=id_length[,c(1,3)]
id_length#ID+length 62754

#GN：ID+name（ensembl下载）
GN=read.table("mart_export.txt",header = T,fill=T, na.strings = c("","NA"))
GN

##查看两个数据框中基因ID是否唯一
if (length(unique(GN$GenestableID)) == nrow(GN)) {
  print("基因 ID 是唯一的。")
} else {
  print("基因 ID 不是唯一的。")
}
   
if (length(unique(id_length$GenestableID)) == nrow(id_length)) {
  print("基因 ID 是唯一的。")
} else {
  print("基因 ID 不是唯一的。")
}


rownames(GN)=GN$GenestableID
GN=GN[sort(rownames(GN)),]
#id_length=id_length[sort(rownames(id_length)),]
#id_length

#rm(exons_gene,exons_gene_lens,gene_length,txdb)

id_GN_length=merge(GN,id_length,by="GenestableID")
#id_GN_length=id_GN_length[,c(1,2,3)]

colnames(id_GN_length)=c("gene_id","gene_name","length")
id_GN_length

#save("id_GN_length",file="./meddata/ID_gene_length.Rdata")

###加载数据(begin)
load("./meddata/ID_gene_length.Rdata")
load("./meddata/mcounts_GSE114007.Rdata")

#基因表达量数据框中基因对应长度
gene_GSE114007=rownames(m_counts_GSE114007)
gene_GSE114007=as.data.frame(gene_GSE114007)
colnames(gene_GSE114007)=c("gene_name")
gene_GSE114007=merge(gene_GSE114007,id_GN_length,by="gene_name")
gene_GSE114007 #20920

tpm_factor =1e6/sum(gene_GSE114007$length/1000)
#gene_name = rownames(m_counts_GSE114007)
#m_counts_GSE114007=cbind(gene_name,m_counts_GSE114007)
as.numeric(gene_GSE114007$length)
as.numeric(tpm_factor)
#as.numeric(m_counts_GSE114007[,-1])
m_counts_GSE114007[,-1]

m_counts_GSE114007[,-1]=m_counts_GSE114007[,-1]*tpm_factor/gene_GSE114007$length
m_tpm_GSE114007=m_counts_GSE114007
save(m_tpm_GSE114007,file="./meddata/m_tpm_GSE114007.RData")


##转化公式
#countToTpm
#countToTpm =function(counts, effLen)
{
  rate = log(counts)-log(effLen)
  denom =log(sum(exp(rate)))
  exp(rate-denom +log(1e6))
}

#FpkmToTpm
#FpkmToTpm=function(fpkm){
  exp(log(fpkm)-log(sum(fpkm))+log(1e6))
}




