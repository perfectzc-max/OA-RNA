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

2024.01.17 ZC  
1.补齐kegg plot function文件  

2024.01.21 LZZ
1.新增DO分析代码 

2024.01.26 LZZ
1.新增step3_DDR_exp文件夹
2.新增1_tpm_boxplot.R文件（基因表达量counts以转换为TPM的方式进行标化）
3.新增DDR_gene.xlsx

2024.01.30 LZZ
1.新增1_DEseq2_boxplot.R文件（基因表达量counts以DEseq2的方式进行标化）

2024.02.01 LZZ
1. 新增1_bind_normalized_boxplot.R(对基因表达量进行quantile normalization[limma中voom方法])

2024.2.24 LZZ
1. 新增0_import_GSE51588.R 数据为标化后的表达量
2. 新增0_bind_data.R 合并GSE114007（标化后）+GSE51588表达矩阵，去批次效应
3. 新增2_GSVA.R  GSE114007（标化后）+GSE51588表达矩阵 做GSVA分析
4. 新增2_GSVA_plot.R GSVA结果作图（热图、火山图、柱状图）
