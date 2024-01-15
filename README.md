# OA-RNA
Analysis RNA-Seq data of OA samples  
2024.01.08 ZC  
1.更新GSE114007导入数据,改正了表格两个sheet导致的数据导入不全  
2.更新DESeq2代码  

2024.01.14 LZZ  
1.0_import_GSE104782.R 删掉保存Group，不存在这个对象  
2.0_import_GSE114007.R 修改表达矩阵  
3.1_PCA_heatmap.R 补充GSE114007的counts转化为tpm的热图  
4.2_DEGs_DESEQ2.R 统一GSE114007不同来源样本的title名称  
5.2_DEGs_limma.R 添加探针注释ids  

2024.01.16 ZC  
1.更新2_DEGs_DESEQ2.R，添加作图函数四个，修改分组  
2.更新2_GO_KEGG.R，富集分析，包含ENTREZ注释步骤，待完善  
