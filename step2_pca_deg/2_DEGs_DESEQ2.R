##差异表达分析
# NOTE：
# 1.目前大家主要用三个包来分析：DESeq2、edgeR和limma。具体介绍：https://www.jianshu.com/p/0165de6fac90
# 2.暂时先用limma，比较方便，DESeq2要求counts矩阵，如果需要的话得从原始数据跑起，比较麻烦。
# 3.对于bulk，可以直接做差异表达分析，但SC样本通常需要先聚类分群注释，再比较。参考：https://zhuanlan.zhihu.com/p/614993353
rm(list = ls()) 
load(file = "./meddata/mcounts_GSE114007.Rdata")
#差异分析,这里有rawcounts，用DESeq2包来做，教程：https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#需要表达矩阵和Group，不需要改
# 加载包
BiocManager::install("pasilla")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
# BiocManager::install("DEGreport",force = T)
install.packages("pheatmap")
install.packages("tidyverse")
library("pasilla")
library("tidyverse")

install.packages("locfit",type="binary")
library("locfit")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")

# library("DEGreport")
# 1.count data+information table
cts =as.matrix(m_counts_GSE114007)
col_114007=  str_to_upper(colnames(m_counts_GSE114007))
colnames(cts)=col_114007
cts=as.matrix(cts)
coldata=pdata_GSE114007$characteristics_ch1.2 %>% as.data.frame()
#coldata=pdata_GSE114007$title %>% as.data.frame()
rownames(coldata)=pdata_GSE114007$title
pdata_GSE114007$title=str_to_upper(pdata_GSE114007$title)
colnames(coldata)=c("condition")
coldata$condition <- factor(coldata$condition)
rownames(coldata)
colnames(cts)
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
#构建对象
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
#加gene列
# featureData <- data.frame(gene=rownames(cts))
# mcols(dds) <- DataFrame(mcols(dds), featureData)
# 过滤掉那些 count 结果都为 0 的数据，这些没有表达的基因对结果的分析没有用
dds <- dds[ rowSums(counts(dds)) > 1, ]

# 2.差异表达分析
dds <- DESeq(dds)
res <- results(dds)
res
# 注: (1)rownames: 基因 (2)baseMean:所有样本矫正后的平均 reads 数 
# (3)log2FoldChange:取 log2 后的表达量差异 (4)pvalue:统计学差异显著性检验指标 
# (5)padj:校正后的 pvalue, padj 越小,表示基因表达差异越显著
summary(res)
plotMA(res, main="DESeq2", ylim=c(-2,2))
#计算上下调基因
sum(res$padj < 0.05, na.rm=TRUE)

# 展示pvalue最小的那个
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
#4.加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("clusterProfiler",force = TRUE)
#install.packages("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(rownames(res), 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
save(res,file = "step2DE__GSE114007.Rdata")
