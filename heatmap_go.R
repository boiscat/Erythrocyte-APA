

# 2021.12.21


cluster_go <- as.data.frame(read_tsv("shorten_go_res.txt") )

cluster_go1 <- data.frame(cbind(cluster_go[,1], -cluster_go[c(2,3,7,4:6)]))

cluster_go1

row.names(cluster_go1) <- cluster_go1$cluster_go...1.

colnames(cluster_go1) <- gsub("\\.","-",colnames(cluster_go1))

pheatmap(cluster_go1[-1],
         cluster_cols = F,
         cluster_rows = T,
         treeheight_row = 0,
         border_color = "grey20",
         cellwidth = 17,cellheight = 12,
         color = viridis(n = 50)[50:1]
         
)


# w=6.6 h= 6.3



cluster_go <- as.data.frame(read_tsv("lengthen_go_res.txt") )

cluster_go1 <- data.frame(cbind(cluster_go[,1], -cluster_go[c(2,3,7,4:6)]))

row.names(cluster_go1) <- cluster_go1$cluster_go...1.

colnames(cluster_go1) <- gsub("\\.","-",colnames(cluster_go1))



pheatmap(cluster_go1[-1],
         cluster_cols = F,
         cluster_rows = T,
         treeheight_row = 0,
         border_color = "grey20",
         cellwidth = 17,cellheight = 12,
         color = viridis(n= 50,option = "A")[50:1]
         
)



# dotplot

cluster_go1 <- pivot_longer(cluster_go[,-c(1,2)],cols = c(1:5) ,names_to = "cluster",values_to = "LogP") %>%
  
  filter(LogP>0)


ggplot(cluster_go1, aes(x=cluster, y=name, size=LogP,fill=LogP)) +
  geom_point(shape=21, color="black") +
  scale_size(range = c(2,7.5), name="Gene Number") +
  ylab("") +
  xlab("")+
  scale_fill_continuous(low="yellow",high="#E600F5",name= expression(paste("Log"[10], "P")))+
  
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.title = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 10, colour ="black"),
        legend.title = element_text(size = 12 ),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "right")

