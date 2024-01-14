rm(list = ls())  
load(file = "./meddata/mcounts_GSE168505.Rdata")
load(file = "./meddata/mcounts_GSE114007.Rdata")
load(file = "./meddata/m_tpm_GSE114007.RData")
#简单的PCA 聚类，只是看一下两组间分布咋样
#如果是SC的样本这一步会更复杂，所以只看BULK 样本
# 有文章会展示这个结果：https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8806218/
#需要表达矩阵和分组信息
Group_1=factor(group_info_GSE168505,levels = unique(group_info_GSE168505))
Group_2=factor(group_info_GSE114007,levels = unique(group_info_GSE114007))
Group_1
Group_2
# 1.PCA 图,看一下两组间的聚类情况
dat_1=as.data.frame(t(m_fpkm_GSE168505))
dat_2=as.data.frame(t(m_counts_GSE114007))
dat_3=as.data.frame(t(m_tpm_GSE114007))
library(FactoMineR)
library(factoextra) 
dat.pca_1 <- PCA(dat_1, graph = FALSE)
dat.pca_2 <- PCA(dat_2, graph = FALSE)
dat.pca_3 <- PCA(dat_3, graph = FALSE)
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
###plot3 GSE114007counts变换tpm后数据
pca_plot_3 <- fviz_pca_ind(dat.pca_3,
                           geom.ind = "point", # show points only (nbut not "text")
                           col.ind = Group_2, # color by groups
                           palette = c("#00AFBB", "#E7B800"),
                           addEllipses = TRUE, # Concentration ellipses
                           legend.title = "Groups"
)
pca_plot_3
# 2.top 10 sd 热图---- 看表达差异最大基因
cg_1=names(tail(sort(apply(m_fpkm_GSE168505,1,sd)),10))
n_1=m_fpkm_GSE168505[cg_1,]

cg_2=names(tail(sort(apply(m_counts_GSE114007,1,sd)),10))
n_2=m_counts_GSE114007[cg_2,]

cg_3=names(tail(sort(apply(m_tpm_GSE114007,1,sd)),10))
n_3=m_tpm_GSE114007[cg_3,]

# 直接画热图
library(pheatmap)
annotation_col_1=data.frame(group=Group_1)
rownames(annotation_col_1)=colnames(n_1) 
pheatmap(n_1,
         show_colnames =T,
         show_rownames = T,
         annotation_col=annotation_col_1
)

annotation_col_2=data.frame(group=Group_2)
rownames(annotation_col_2)=colnames(n_2) 
pheatmap(n_2,
         show_colnames =T,
         show_rownames = T,
         annotation_col=annotation_col_2
)

annotation_col_3=data.frame(group=Group_2)
rownames(annotation_col_3)=colnames(n_3)
pheatmap(n_3,
         show_colnames =T,
         show_rownames = T,
         annotation_col=(annotation_col_3)
)

# 按行标准化，只按行比较
pheatmap(n_1,
         show_colnames =T,
         show_rownames = T,
         annotation_col=annotation_col_1,
         scale = "row",
         breaks = seq(-3,3,length.out = 100) #数据分配范围
)

pheatmap(n_2,
         show_colnames =T,
         show_rownames = T,
         annotation_col=annotation_col_2,
         scale = "row",
         breaks = seq(-3,3,length.out = 100) #数据分配范围
)

pheatmap(n_3,
         show_colnames =T,
         show_rownames = T,
         annotation_col=annotation_col_3,
         scale = "row",
         breaks = seq(-3,3,length.out = 100) #数据分配范围
)

