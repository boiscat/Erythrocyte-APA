library(tidyverse)  
 #install.packages("ggrepel")
library(ggrepel)

data <- read_table("/home/yangyb/yangyb/Erythro_APA_Liu/linuxnd/analysis/expression/pabpc1_1_vs_control_deseq_raw.txt")

colnames(data) <- gsub('\"','',colnames(data))

data$gene <- gsub('\"','',data$gene)

test <- read_table("~/test/ttt") %>% .[c(1:2)]

data1 <- merge(data,test,by="gene",all = T)




my_colors <- c("#E55858","#3C71A3",  "#66ABAF",  "#FE8E10", "#47964C", "#B47A9C","#FCC739")


# 要标注的基因
highlight_genes <- c( "PSMA5", "PTP4A2",'PABPC1','TSC22D1','IDH2','DDX17','P2RX1','RSL1D1','PARP1','ZFP36L2')

highlight_genes <- c('PABPC1','TSC22D1','IDH2','DDX17','P2RX1','RSL1D1','PARP1','ZFP36L2')


data1$UpDown <- ifelse(data1$log2FoldChange > 0.3 & data1$pvalue < 0.05, "Up-regulated", 
                      ifelse(data1$log2FoldChange < -0.3 & data1$pvalue < 0.05, "Down-regulated", "NotSig"))

# 创建火山图

# 5 3.5
ggplot(data1, aes(x = log2FoldChange, y = -log10(pvalue), color = UpDown)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_text_repel(
    data = subset(data1, symbol %in% highlight_genes),
    aes(label = symbol, color = UpDown),
    box.padding = 1,
    point.padding = 0.3,
    segment.color = "red",
    segment.size = 0.5,
    color = 'red'
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#3C71A3") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "#3C71A3") +
  scale_color_manual(values = c("NotSig" = "gray", "Up-regulated" = "#FE8E10", "Down-regulated" = "#66ABAF")) +
  theme_test() +
  labs(
    x = expression(Log[2](FoldChange)),
    y = expression(-log[10](pvalue))
  )

table(data1$UpDown
      )
