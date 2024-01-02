#下载人类基因组gtf数据
##https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
library(GenomicFeatures)
txdb=makeTxDbFromGFF("Homo_sapiens.GRCh38.110.gtf",format = "gtf")
exons_gene=exonsBy(txdb,by="gene")
exons_gene_lens=lapply(exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_lens[1:10]
#得到基因长度及对应ID
gene_length=sapply(exons_gene_lens,function(x){x})
id_length=as.data.frame(gene_length)
##id_length$Geneid=rownames(id_length)
id_length #id_length 62754 obs.

#ID与基因名
##下载https://www.ensembl.org/biomart/martview/dd4d77d5cfb1aea321820bb2d488d637
##选择GenestableID和Gene name导出“mart_exprot.txt”
GN=read.table("mart_export.txt",header = T,fill=T, na.strings = c("","NA"))
GN #GN 70116 obs.
#GN=na.omit(GN)
#GN ##omit 47776 obs.
rownames(GN)=GN$GenestableID
GN=GN[sort(rownames(GN)),]
#id_length=id_length[sort(rownames(id_length)),]
#id_length
#合并数据框？？行数不同
id_GN_length=cbind(GN,id_length)
#修改列名
colnames(id_GN_length)=c("gene_id","gene_name","length")
id_GN_length


#countToTpm
countToTpm =function(counts, effLen)
{
  rate = log(counts)-log(effLen)
  denom =log(sum(exp(rate)))
  exp(rate-denom +log(1e6))
}

#FpkmToTpm
FpkmToTpm=function(fpkm){
  exp(log(fpkm)-log(sum(fpkm))+log(1e6))
