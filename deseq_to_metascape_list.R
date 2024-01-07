test <- read_tsv('~/test/ttt')
data <- read_table("/home/yangyb/yangyb/Erythro_APA_Liu/linuxnd/analysis/expression/pabpc1_1_vs_control_deseq_raw.txt") 

colnames(data) <- gsub('\"','',colnames(data))

data$gene <- gsub('\"','',data$gene)
data$gene <- gsub('\\..*$','',data$gene)


data1 <- data 
data1$UpDown <- ifelse(data1$log2FoldChange > 0.3 & data1$pvalue < 0.05, "Up-regulated", 
                       ifelse(data1$log2FoldChange < -0.3 & data1$pvalue < 0.05, "Down-regulated", "NotSig"))

data2 <- data1 %>% filter(UpDown!="NotSig") %>% .[c('gene','UpDown')]



# 获取唯一的 UpDown 值
updown_values <- unique(data2$UpDown)

# 创建一个新的数据框存储结果
result_df <- data.frame(matrix(ncol = 1 + length(data2$gene), nrow = length(updown_values)))
colnames(result_df) <- c("UpDown", data2$gene)

# 填充数据
for (i in 1:length(updown_values)) {
  result_df[i, 1] <- updown_values[i]
  for (j in 2:(length(data2$gene) + 1)) {
    result_df[i, j] <- ifelse(data2$UpDown[j - 1] == updown_values[i], data2$gene[j - 1], "")
  }
}

# 将数据框转换为字符向量
result_vector <- apply(result_df, 1, function(x) paste(x[x != ""], collapse = ","))

# 打印结果
cat(result_vector, sep = "\n")


writeLines(result_vector, con = '/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/03-expression/deseq_go.txt')
