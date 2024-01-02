rm(list = ls())  
load(file = "./meddata/mcounts_GSE168505.Rdata")
load(file = "./meddata/mcounts_GSE114007.Rdata")
#简单的PCA 聚类，只是看一下两组间分布咋样
#如果是SC的样本这一步会更复杂，所以只看BULK 样本
# 有文章会展示这个结果：https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8806218/
#需要表达矩阵和分组信息
Group_1=factor(group_info_GSE168505,levels = unique(group_info_GSE168505))
Group_2=factor(group_info_GSE114007,levels = unique(group_info_GSE114007))
# 1.PCA 图,看一下两组间的聚类情况
dat_1=as.data.frame(t(m_fpkm_GSE168505))
dat_2=as.data.frame(t(m_counts_GSE114007))
library(FactoMineR)
library(factoextra) 
dat.pca_1 <- PCA(dat_1, graph = FALSE)
dat.pca_2 <- PCA(dat_2, graph = FALSE)
#作图
pca_plot_1 <- fviz_pca_ind(dat.pca_1,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Group_1, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)
pca_plot_1

pca_plot_2 <- fviz_pca_ind(dat.pca_2,
                           geom.ind = "point", # show points only (nbut not "text")
                           col.ind = Group_2, # color by groups
                           palette = c("#00AFBB", "#E7B800"),
                           addEllipses = TRUE, # Concentration ellipses
                           legend.title = "Groups"
)
pca_plot_2
# 2.top 10 sd 热图---- 看表达差异最大基因
cg=names(tail(sort(apply(m_fpkm_GSE168505,1,sd)),10))
n=m_fpkm_GSE168505[cg,]

# 直接画热图
library(pheatmap)
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(n) 
pheatmap(n,
         show_colnames =T,
         show_rownames = T,
         annotation_col=annotation_col
)

# 按行标准化，只按行比较
pheatmap(n,
         show_colnames =T,
         show_rownames = T,
         annotation_col=annotation_col,
         scale = "row",
         breaks = seq(-3,3,length.out = 100) #数据分配范围
)
