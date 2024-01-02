##差异表达分析
# NOTE：
# 1.目前大家主要用三个包来分析：DESeq2、edgeR和limma。具体介绍：https://www.jianshu.com/p/0165de6fac90
# 2.暂时先用limma，比较方便，DESeq2要求counts矩阵，如果需要的话得从原始数据跑起，比较麻烦。
# 3.对于bulk，可以直接做差异表达分析，但SC样本通常需要先聚类分群注释，再比较。参考：https://zhuanlan.zhihu.com/p/614993353
rm(list = ls()) 
load(file = "./meddata/mcounts_GSE168505.Rdata")
Group_1=factor(group_info_GSE168505,levels = unique(group_info_GSE168505))
#差异分析，用limma包来做
#需要表达矩阵和Group，不需要改
library(limma)

#样本信息注释
design=model.matrix(~Group_1)
fit=lmFit(m_fpkm_GSE168505,design)

#差异分析
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
head(deg)
#为deg数据框添加几列
#1.加probe_id列，把行名变成一列
library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)
#2.加上探针注释
#一个基因对应多个探针，可以取最大值、随机去重、平均值
ids = ids[!duplicated(ids$symbol),]
#其他去重方式在zz.去重.R
deg <- inner_join(deg,ids,by="probe_id")
head(deg)
nrow(deg)

#3.加change列,标记上下调基因
logFC_t=1
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)


#4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
dim(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
dim(deg)
length(unique(deg$symbol))
save(Group,deg,logFC_t,P.Value_t,gse_number,file = "step4output.Rdata")
