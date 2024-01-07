## ---------------------------
##
## Purpose of script: 表达分析（TPM）
##
## R version: 3.6.1
##
## Author: Dr. Yanbo Yang
##
## Date Created: 2022-01-09
##
## Copyright (c) Jimmy Yang, 2022
##
## Email: 404495360@qq.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


setwd("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/01-data/hsc1_proe1_APA/03-expression/")

raw_count <- read_tsv("raw_count.txt")

raw_count1 <- raw_count[,-grep("Geneid",colnames(raw_count))]

colnames(raw_count1)[1:3] <- "Pro-E"
colnames(raw_count1)[4:6] <- "Early-E"
colnames(raw_count1)[7:9] <- "Late-E"
colnames(raw_count1)[10:12] <- "Poly"
colnames(raw_count1)[13:15] <- "Ortho"
colnames(raw_count1)[16:18] <- "BFU-E"
colnames(raw_count1)[19:21] <- "CFU-E"

raw_count2 <- cbind(gene=raw_count$Geneid...1,raw_count1)

library(GenomicFeatures)

txdb <- makeTxDbFromGFF("~/yangyb/Erythro_APA_Liu/02-reference/Homo_sapiens.GRCh37.87.gtf",format="gtf")

exons_gene <- exonsBy(txdb, by = "gene")

gene_length <- lapply(exons_gene,function(x){sum(width(reduce(x)))})

gene_length1 <- as.data.frame(cbind(Geneid=names(unlist(gene_length)),Length=as.numeric(unlist(gene_length))))

gene_length1 <- gene_length1[match(raw_count2$gene,gene_length1$Geneid),]

gene_length1$Geneid <- as.character(gene_length1$Geneid)

gene_length1$Length <- as.numeric(as.character(gene_length1$Length))


#03 calculate TPM
library(snowfall) 

library(doSNOW)


cl <- makeCluster(40)
registerDoSNOW(cl)
iterations <- dim(raw_count1)[1]
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


tpm <- raw_count1


# 第一步，用长度矫正count

tpm1<-foreach(i = 1:iterations, .options.snow = opts) %dopar%
  
  {
    tpm[i,] <- tpm[i,]/gene_length1$Length[i]
  }

tpm1<- do.call(rbind.data.frame, tpm1) 

close(pb)

stopCluster(cl) 


# 第二步，计算矫正后count总值sum(x)

# 第三步x/sum(x)

tpm2 <- apply(tpm1,2,function(x){
  
  x/sum(x)*10^6
  
})

tpm_res <- cbind(gene=raw_count$Geneid...1,tpm2[,c(16:21,1:15)])

write.table(tpm_res,"/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/01-data/hsc1_proe1_APA/03-expression/tpm_matrix.txt",quote = F,sep = '\t',row.names = F)



