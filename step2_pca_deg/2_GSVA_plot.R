rm(list = ls())

load("./meddata/gsva_114007+51588_go_matrix.R")
load("./meddata/gsva_114007+51588_kegg_matrix.R")

##limma
design = c(rep("NORMAL",18),rep("OA",20),rep("NORMAL",10),rep("OA",40)) %>% factor(.,levels=c("NORMAL","OA"),ordered=F)
design = model.matrix(~0+factor(design))
head(design)
colnames(design)=c("NORMAL","OA")
rownames(design)=colnames(gsva_bind_go_mat)

contrast.matrix=makeContrasts("OA-NORMAL",
                              levels = design)
fit1=lmFit(gsva_bind_go_mat,design) #fitting model 
fit2=contrasts.fit(fit1,contrast.matrix) #statistical tests
efit=eBayes(fit2)

res=decideTests(efit,p.value = 0.05)
summary(res) #result
tempOutput=topTable(efit,n=Inf)  #10395
degs=na.omit(tempOutput) #10395

degs_bind_go=degs
save(design,coldata,degs_bind_go,file = "./meddata/step4_gsva_gobind_degs.result.R")

#visualization
##heatmap
padj_cutoff=0.01
log2FC_cutoff=log2(1.2)

keep <- rownames(degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC)>log2FC_cutoff, ])
length(keep) 
load("./meddata/gsva_114007+51588_go_matrix.R")
dat <- gsva_bind_go_mat[keep[1:50],] #top 50
load("./meddata/step4_gsva_gobind_degs.result.R")
gl=coldata[,3]
gl=as.data.frame(gl)
rownames(gl)=rownames(coldata)
colnames(gl)="group"
gl

pheatmap::pheatmap(dat, 
                   fontsize_row = 8,
                   height = 10,
                   width=18,
                   annotation_col = gl,
                   show_colnames = F,
                   show_rownames = T,
                   filename = paste0('GSVA_bind_go_heatmap.pdf'))

#####火山图
degs_bind_go$significance  <- as.factor(ifelse(degs_bind_go$adj.P.Val < padj_cutoff & abs(degs_bind_go$logFC) > log2FC_cutoff,
                                                 ifelse(degs_bind_go$logFC > log2FC_cutoff ,'UP','DOWN'),'NOT'))

this_title <- paste0(' Up :  ',nrow(degs_bind_go[degs_bind_go$significance =='UP',]) ,
                     '\n Down : ',nrow(degs_bind_go[degs_bind_go$significance =='DOWN',]),
                     '\n adj.P.Val <= ',padj_cutoff,
                     '\n FoldChange >= ',round(2^log2FC_cutoff,3))
#'\n是换行'
this_title

g <- ggplot(data=degs_bind_go, 
            aes(x=logFC, y=-log10(adj.P.Val),
                color=significance)) +
  #点和背景
  geom_point(alpha=0.4, size=1) +
  theme_classic()+ #无网格线
  #坐标轴
  xlab("log2 ( FoldChange )") + 
  ylab("-log10 ( adj.P.Val )") +
  #标题文本
  ggtitle( this_title ) +
  #分区颜色                   
  scale_colour_manual(values = c('blue','grey','red'))+ 
  #辅助线
  geom_vline(xintercept = c(-log2FC_cutoff,log2FC_cutoff),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(padj_cutoff),lty=4,col="grey",lwd=0.8) +
  #图例标题间距等设置
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin=unit(c(2,2,2,2),'lines'), #上右下左
        legend.title = element_blank(),
        legend.position="right")
g
ggsave(g,filename = 'GSVA_bind_go_volcano_padj.pdf',width =8,height =7.5)

#发散条图
#### 发散条形图绘制 ####
library(tidyverse)  # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(ggthemes)
library(ggprism)
p_cutoff=0.001

#载入gsva的差异分析结果
contrast.matrix=makeContrasts("OA-NORMAL",
                              levels = design)
fit1=lmFit(gsva_bind_kegg_mat,design) #fitting model 
fit2=contrasts.fit(fit1,contrast.matrix) #statistical tests
efit=eBayes(fit2)

res=decideTests(efit,p.value = 0.05)
summary(res) #result
tempOutput=topTable(efit,n=Inf)  #186
degs=na.omit(tempOutput) #186

degs_bind_kegg=degs
save(design,coldata,degs_bind_kegg,file = "./meddata/step4_gsva_keggbind_degs.result.R")


#degs <- degs_bind_kegg  
Diff <- rbind(subset(degs_bind_kegg,logFC>0)[1:20,], subset(degs_bind_kegg,logFC<0)[1:20,]) #选择上下调前20通路     
dat_plot <- data.frame(id  = row.names(Diff),
                       p   = Diff$P.Value,
                       lgfc= Diff$logFC)
dat_plot$group <- ifelse(dat_plot$lgfc>0 ,1,-1)    # 将上调设为组1，下调设为组-1
dat_plot$lg_p <- -log10(dat_plot$p)*dat_plot$group # 将上调-log10p设置为正，下调-log10p设置为负
dat_plot$id <- str_replace(dat_plot$id, "KEGG_","");dat_plot$id[1:10]
dat_plot$threshold <- factor(ifelse(abs(dat_plot$p) <= p_cutoff,
                                    ifelse(dat_plot$lgfc >0 ,'Up','Down'),'Not'),
                             levels=c('Up','Down','Not'))

dat_plot <- dat_plot %>% arrange(lg_p)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)

## 设置不同标签数量
low1 <- dat_plot %>% filter(lg_p < log10(p_cutoff)) %>% nrow()
low0 <- dat_plot %>% filter(lg_p < 0) %>% nrow()
high0 <- dat_plot %>% filter(lg_p < -log10(p_cutoff)) %>% nrow()
high1 <- nrow(dat_plot)

p <- ggplot(data = dat_plot,aes(x = id, y = lg_p, 
                                fill = threshold)) +
  geom_col()+
  coord_flip() + 
  scale_fill_manual(values = c('Up'= '#36638a','Not'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-log10(p_cutoff),log10(p_cutoff)),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('-log10(P.Value) of GSVA score') + 
  guides(fill="none")+
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'black') + #黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 黑色标签
p
ggsave("GSVA_bind_barplot_pvalue.pdf",p,width = 15,height  = 15)
