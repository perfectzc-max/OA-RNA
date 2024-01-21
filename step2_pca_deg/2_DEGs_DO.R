##DO分析
rm(list = ls())
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
load("D:/LZZ/analysis/local_rna_seq/meddata/step2_DE_GSE114007.Rdata")

#dif_gene=dat[dat$change != 'stable',] #7526
#dif_gene=dif_gene[,c(1,3)]
#write.table(dif_gene,sep = "\t",row.names = F,col.names = T )

#上调基因
up_gene=dat[dat$change=="up",]
up_gene=up_gene[,c(1,3)]
#write.table(up_gene,sep = "\t",row.names = F,col.names = T )
up_genes=as.vector(up_gene[,1])
up_entrezIDs=mget(up_genes,org.Hs.egSYMBOL2EG,ifnotfound = NA)
up_entrezIDs=as.character(up_entrezIDs)
up_entrezIDs=cbind(up_gene,entreID=up_entrezIDs)#5483
#去除基因id为NA的基因
up_entrezIDs=up_entrezIDs[up_entrezIDs$entreID != 'NA',] #4801
write.csv(up_entrezIDs,file="upgene_entreID.csv")
#100187828#100505381

#DO
library(DOSE)
##手动删除csv文件中一个基因对应多个ID的行
up_gene=read.csv("upgene_entreID.csv",sep = ",",header = T)
up_gene=up_gene[,-1]
up_erich.do = DOSE::enrichDO(gene=up_gene$entreID,
                          ont = "DO",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = T)
bar_up_DO=barplot(up_erich.do)
bar_up_DO
dot_up_DO=dotplot(up_erich.do)
dot_up_DO
web_up_DO=cnetplot(up_erich.do)
web_up_DO
class(bar_up_DO)
#合并图片
up_DO=cowplot::plot_grid(bar_up_DO,dot_up_DO,web_up_DO,ncol = 3)
ggsave("up_DO.png",width = 1667 ,height =  502)

#下调基因
down_gene=dat[dat$change=="down",]
down_gene=down_gene[,c(1,3)]
#write.table(down_gene,sep = "\t",row.names = F,col.names = T )
down_genes=as.vector(down_gene[,1])
down_entrezIDs=mget(down_genes,org.Hs.egSYMBOL2EG,ifnotfound = NA)
down_entrezIDs=as.character(down_entrezIDs)
down_entrezIDs=cbind(down_gene,entreID=down_entrezIDs)#2043
#去除基因id为NA的基因
down_entrezIDs=down_entrezIDs[down_entrezIDs$entreID != 'NA',] #1799
write.csv(down_entrezIDs,file="downgene_entreID.csv")
#DO
library(DOSE)
##手动删除csv文件中一个基因对应多个ID的行
down_gene=read.csv("downgene_entreID.csv",sep = ",",header = T)
down_gene=down_gene[,-1]
down_erich.do = DOSE::enrichDO(gene=down_gene$entreID,
                             ont = "DO",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = T)
bar_down_DO=barplot(down_erich.do)
bar_down_DO
dot_down_DO=dotplot(down_erich.do)
dot_down_DO
web_down_DO=cnetplot(down_erich.do)
web_down_DO
class(bar_down_DO)

down_DO=cowplot::plot_grid(bar_down_DO,dot_down_DO,web_down_DO,ncol = 3)
down_DO
ggsave("down_DO.png",width = 1667 ,height =  502)
