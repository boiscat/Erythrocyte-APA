
library(psych)
library(ggplot2)
library(factoextra)
library(readr)
library(tidyr)
library(ggforce)

count <- read_tsv("/home/yangyb/yangyb/Erythro_APA_Liu/linuxnd/analysis/lapa/res.count") %>% dplyr::select(1,ends_with("bam"))

# count质检：

# PCA分析 （看count数据是否能分开）----

  
  #从指定路径读取名为res.count的tsv文件，并选取第一列和以“bam”结尾的列
    count <- read_tsv("/home/yangyb/yangyb/Erythro_APA_Liu/linuxnd/analysis/lapa/res.count") %>% dplyr::select(1, ends_with("bam"))
  
  #去掉第一列，并筛选出除第一列外每行求和值不为0的数据作为df
    df <- count[-1][apply(count[-1],1,sum)!=0,]
  
  #对df进行转置，将原本每个bam文件的计数值按行变为按列
    df_t = t(df)
  
  #取df_t的列数，即表示变量的数量
    variableL <- ncol(df_t)
  
  #对df_t进行主成分分析，并将各变量标准化（scale=T）
    pca <- prcomp(df_t[, 1:variableL], scale=F)
  
  #对主成分分析结果进行可视化，主要是展示解释变量差异的程度
    fviz_eig(pca, addlabels = TRUE) + theme(axis.text = element_text(size = 12, colour = "black"))
  
  #对主成分分析结果进行可视化，将每个样本在第一主成分和第二主成分的坐标展示
    fviz_pca_ind(pca, repel = F, label = T,habillage=row.names(pca$x))
  
    row.names(pca$x)

  
# sample_list ----
  
library(DESeq2)
library(tidyverse)

count <- data.frame(count)
row.names(count) = count$Geneid...1


#也有用 rowSums(counts(dds)) >= 10, 目的：过滤低表达（似乎10是常用的标准）



count1 <- count[-1][apply(count[-1],1,sum)>=10,c(1:6)]



# 分组采样
condition <- c(rep("control",3),rep("treat",3))

colData <- data.frame(row.names = colnames(count1),condition)

dds <- DESeqDataSetFromMatrix(count1, colData, design= ~ condition)

# 使用lfcShrink,betaPrior

dds <- DESeq(dds,betaPrior=TRUE)

res <- results(dds, contrast=c("condition","treat","control"))

res1 <- data.frame(res) %>%  mutate(gene=row.names(res)) %>% filter(pvalue<0.05) 

# 加载文件

ref <- read_tsv("~/yangyb/Erythro_APA_Liu/linuxnd/analysis/ref/ensg_symbol",col_names = F)

colnames(ref) <- c("gene","symbol")

res2 <- merge(ref,res1)


# 标准化

dds.sizefactor <- estimateSizeFactors(dds)  

sizeFactors(dds.sizefactor )

count2 <- counts(dds, normalized=TRUE)

count3 <- cbind(row.names(count2),count2)

write.table(res2,"/home/yangyb/yangyb/Erythro_APA_Liu/linuxnd/analysis/expression/pabpc1_1_vs_control_deseq_raw.txt",row.names = F)
