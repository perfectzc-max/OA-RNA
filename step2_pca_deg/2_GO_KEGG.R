##GO analysis
rm(list = ls())  
load("./meddata/step2_DE_GSE114007.Rdata")
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(DOSE)
# 1.GO 富集分析
#加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
s2e <- bitr(dat$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
dat <- inner_join(dat,s2e,by=c("symbol"="SYMBOL"))

#输入数据：只要ENTREZID这一列就行
gene_up = dat[dat$change == 'up','ENTREZID'] 
gene_down = dat[dat$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)

#富集：go分析
#以下步骤耗时很长，设置了存在即跳过
if(!file.exists("./meddata/step2_GO_GSE114007.Rdata")){
  ego <- enrichGO(gene = gene_diff,
                  OrgDb= org.Hs.eg.db,
                  ont = "ALL",
                  readable = TRUE)
  ego_BP <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     readable = TRUE)
  #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
  save(ego,ego_BP,file = "./meddata/step2_GO_GSE114007.Rdata")
}
load("./meddata/step2_GO_GSE114007.Rdata")

#(3)可视化
#条带图
barplot(ego)
#气泡图
dotplot(ego)
dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., scale = "free") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))

#geneList 用于设置下面图的颜色
geneList = dat$log2FoldChange
names(geneList)=dat$ENTREZID

#(3)展示top通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(ego)
# 需要debug,定义两个函数
dropAsis <- function(x){
  cls <- class(x)
  structure(x, class = setdiff(cls, "AsIs"))
}
rescale.AsIs <- function(x, ...){
  # 自定义dropAsis方法
  dropAsis <- function(x){
    cls <- class(x)
    structure(x, class = setdiff(cls, "AsIs"))
  }
  # 调用本来的rescale方法
  scales:::rescale(dropAsis(x), ...)
}

cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(ego, showCategory = 3,foldChange=geneList, circular = TRUE, colorEdge = TRUE)
#Enrichment Map,这个函数最近更新过，版本不同代码会不同

Biobase::package.version("enrichplot")

if(T){
  emapplot(pairwise_termsim(ego)) #新版本
}else{
  emapplot(ego)#旧版本
}
#(4)展示通路关系 https://zhuanlan.zhihu.com/p/99789859
#goplot(ego)
goplot(ego_BP)

#(5)Heatmap-like functional classification
heatplot(ego,foldChange = geneList,showCategory = 8)


# 2.KEGG pathway analysis----
#上调、下调、差异、所有基因
#（1）输入数据
gene_up = dat[dat$change == 'up','ENTREZID'] 
gene_down = dat[dat$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)

#（2）对上调/下调/所有差异基因进行富集分析
# 会有各种网络问题：https://zhuanlan.zhihu.com/p/534214175
# 包版本问题：https://blog.csdn.net/weixin_57699783/article/details/129531280
# devtools::install_github("YuLab-SMU/clusterProfiler",force = TRUE)
# options(clusterProfiler.download.method = "wininet")
if(!file.exists(paste0("./meddata/step2_KEGG_GSE114007.Rdata"))){
  kk.up <- enrichKEGG(gene = gene_up,
                      organism ="hsa")
  kk.down <- enrichKEGG(gene = gene_down,
                        organism='hsa')
  kk.diff <- enrichKEGG(gene= gene_diff,
                        organism='hsa')
  save(kk.diff,kk.down,kk.up,file = paste0(gse_number,"_KEGG.Rdata"))
}
load("./meddata/step2_KEGG_GSE114007.Rdata")

#(3)看看富集到了吗？https://mp.weixin.qq.com/s/NglawJgVgrMJ0QfD-YRBQg
table(kk.diff@result$p.adjust<0.05)
table(kk.up@result$p.adjust<0.05)
table(kk.down@result$p.adjust<0.05)
#富集到两个

#(4)双向图
# 富集分析所有图表默认都是用p.adjust,富集不到可以退而求其次用p值，在文中说明即可
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)

source("kegg_plot_function.R")
g_kegg <- kegg_plot(up_kegg,down_kegg)
g_kegg
#g_kegg +scale_y_continuous(labels = c(2,0,2,4,6))
ggsave(g_kegg,filename = 'kegg_up_down.png')

# 3.另一个常用包：gsea作kegg和GO富集分析
#https://www.yuque.com/xiaojiewanglezenmofenshen/dbwkg1/ytawgg

#(1)查看示例数据
data(geneList, package="DOSE")
#(2)将我们的数据转换成示例数据的格式
geneList=dat$log2FoldChange
names(geneList)=dat$ENTREZID
geneList=sort(geneList,decreasing = T)
#(3)富集分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  verbose      = FALSE)
down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
#(4)可视化
g2 = kegg_plot(up_kegg,down_kegg)
g2