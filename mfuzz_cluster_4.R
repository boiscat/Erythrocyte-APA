
#BiocManager::install("Mfuzz")
#BiocManager::install("marray")
library(Mfuzz)
library(magrittr)


##Load file 1.total PPAU 2. sig PPAU


pdui <- read_tsv("~/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/04-mfuzz/02-diff_0.15_for_mfuzz_sorted_row_filled.txt") 

colnames(pdui)[-1]

##matrix file and filter

count <- data.matrix(pdui[,-1])

row.names(count) <- pdui$Gene

eset <- new("ExpressionSet",exprs = count)

eset <- filter.std(eset,min.std=0)

eset <- standardise(eset)

# choose optimal cluster number by Gap statistic method

library(factoextra)
library(NbClust)

#fviz_nbclust(count, kmeans, method = "gap_stat",nboot = 500) + labs(subtitle = "Gap statistic method")

#elbow method (alternative)

fviz_nbclust(count, kmeans, method = "wss") +labs(subtitle = "Elbow method")

c <- 4

m <- mestimate(eset)

cl <- mfuzz(eset, c = c, m = m)

# 查看每个cluster中的基因个数

cl$size

# 提取某个cluster下的基因

cl$cluster[cl$cluster == 8]




write.table(data.frame(cl$cluster),"~/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/04-mfuzz/mfuzz_apa_events_filtered_na_gene_list.txt",quote = F,sep = '\t')

# 查看基因和cluster之间的membership

a<- cl$membership

####plot

mfuzz.plot2(eset, cl, mfrow = c(2,2),
            min.mem     =0.8,
            colo = "fancy",
            col.sub = "blue",
            col = "blue",
            centre= TRUE,
            centre.lwd=3,
            cex.main = 1.5,
            cex.lab=1.5,
            ylab = "Relative PDUI level",
            #            single = 4,
            x11 = FALSE,
            time.labels = c('BFU-E','CFU-E','Pro-E','Early-E','Late-E','Poly','Ortho')
)




#add colorbar
mfuzzColorBar(main="Membership", cex.main=1, col = "fancy" )


#save file

setwd("/home/yangyb/cardiomyo_diff_APA/04-pdui/02-mfuzz")

for(i in 1:5){
  
  write.table(str_split(names(cl$cluster[cl$cluster==i]),"\\|",simplify = T)[,2],paste0("cluster_",i,"_gene.txt"),sep = "\t",row.nam
              es = F,col.names = F,quote = F)
  
}

write.table(cl$membership,"cluster_membership",sep = "\t",row.names = T,quote = F)


test <- as.data.frame(cl$membership)

test1 <- test[rev(order(test$`3`)),]
