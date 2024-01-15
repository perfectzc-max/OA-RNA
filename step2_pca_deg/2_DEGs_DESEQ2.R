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
# BiocManager::install("pasilla")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# # BiocManager::install("DEGreport",force = T)
# install.packages("pheatmap")
# install.packages("tidyverse")
library("pasilla")
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library(dplyr)
library(ggplot2)
library(corrplot)
library(paletteer)
library(cowplot)
library(patchwork)
library(ggplotify)
# library("DEGreport")
# 1.count data+information table
cts =as.matrix(m_counts_GSE114007)
coldata=pdata_GSE114007$characteristics_ch1.2 %>% as.data.frame()
rownames(coldata)=pdata_GSE114007$title
colnames(coldata)=c("condition")
coldata$condition <- factor(coldata$condition,levels = c("oa grade: 1","oa grade: 4"))
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
save(res,file = "step2_DE_GSE114007.Rdata")

# 简单看看情况
plotMA(res, main="DESeq2", ylim=c(-2,2))
#计算上下调基因
sum(res$padj < 0.05, na.rm=TRUE)
# 展示pvalue最小的那个
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# 3.作图
#3.1.火山图
#3.加change列,标记上下调基因
dat=as.data.frame(res) %>% na.omit()
dat=rownames_to_column(dat,var = "symbol")
dat=dat[!duplicated(dat$symbol),]
logFC_t=1
P.Value_t = 0.05
k1 = (dat$padj < P.Value_t)&(dat$log2FoldChange < -logFC_t)
k2 = (dat$padj < P.Value_t)&(dat$log2FoldChange > logFC_t)
dat <- mutate(dat,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(dat$change)

# 作图
p <- ggplot(data = dat, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p


# 可以按需求标记值得注意的gene，如：
if(T){
  #自选基因
  for_label <- dat%>% 
    filter(symbol %in% c("HADHA","LRRFIP1")) 
}
if(T){
  #p值最小的10个
  for_label <- dat %>% head(10)
}
if(F) {
  #p值最小的前3下调和前3上调
  x1 = dat %>% 
    filter(change == "up") %>% 
    head(3)
  x2 = dat %>% 
    filter(change == "down") %>% 
    head(3)
  for_label = rbind(x1,x2)
}

#加标签
#在点上加了一个空心圈
#加上了标签
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot


#3.2.差异基因热图
exp = m_counts_GSE114007[dat$symbol,]
rownames(exp) = dat$symbol
if(F){
  #全部差异基因
  cg = dat$symbol[dat$change !="stable"]
  length(cg)
}else{
  #取前10上调和前10下调
  x=dat$log2FoldChange[dat$change !="stable"] 
  names(x)=dat$symbol[dat$change !="stable"] 
  cg=names(c(head(sort(x),10),tail(sort(x),10)))
  length(cg)
}
n=exp[cg,]
dim(n)

#差异基因热图
library(pheatmap)
annotation_col=data.frame(group=coldata)
rownames(annotation_col)=colnames(n) 
heatmap_plot <- pheatmap(n,show_colnames =F,
                         scale = "row",
                         #cluster_cols = F, 
                         annotation_col=annotation_col,
                         breaks = seq(-3,3,length.out = 100)
) 
heatmap_plot

# 3.3.感兴趣基因的箱线图
g = cg
dat = t(exp[g,]) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group = coldata)

pdat = dat%>% 
  pivot_longer(cols = 2:(ncol(dat)-1),
               names_to = "gene",
               values_to = "count")

pdat$gene = factor(pdat$gene,levels = cg,ordered = T)
pdat$change = ifelse(pdat$gene %in% head(cg,10),"down","up")
box_plot = ggplot(pdat,aes(gene,count))+
  geom_boxplot(aes(fill = group$condition))+
  scale_fill_manual(values = c("blue","red"))+
  #scale_fill_paletteer_d("basetheme::minimal")+
  geom_jitter()+
  theme_bw()+
  facet_wrap(~change,scales = "free")
box_plot

# 3.4.感兴趣基因的相关性
M = cor(t(exp[g,]))
pheatmap(M)
my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10)
corrplot(M, type="upper",
         method="pie",
         order="hclust", 
         col=my_color,
         tl.col="black", 
         tl.srt=45)
cor_plot <- recordPlot() 
