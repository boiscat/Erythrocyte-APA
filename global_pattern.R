pdui <- read_tsv("~/yangyb/Erythro_APA_Liu/merged_APA_20211208/01-data/hsc1_proe1_APA/pdui_res/pdui_result.txt")

pdui1 <- pdui[,-c(1:4)]


median_val <- apply(pdui1, 2, function(x){median(x,na.rm = T)})

data <- data.frame(names(median_val),median_val)


data$names.median_val.[1:3] <- "Pro-E"
data$names.median_val.[4:6] <- "Early-E"
data$names.median_val.[7:9] <- "Late-E"
data$names.median_val.[10:12] <- "Poly"
data$names.median_val.[13:15] <- "Ortho"
data$names.median_val.[16:18] <- "BFU-E"
data$names.median_val.[19:21] <- "CFU-E"

colnames(data)[1] <- "sample"

data1 <- data[c(16:21,1:15),] 

data1$sample <- factor(data1$sample,level=unique(data1$sample))

# 3.6 4.4

library(ggsci)

ggplot(data1,aes(x=sample,y=median_val,fill=sample))+geom_boxplot(width=0.3)+
  
theme_classic()+labs(x="",y="Relative 3'UTR Length (%)")+ylim(c(0.2,0.8))+scale_fill_nejm()+
  
theme(axis.text.x = element_text( size=10,colour = "black"),
      axis.text.y = element_text(size=10,colour = "black"),
      title = element_text(size=12,colour = "black"))+ guides(fill=FALSE)



#7 5

tmp_short <- table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[13:15], 1, function(x){median(x,na.rm=T)})> 0.2)
tmp_long <- table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[13:15], 1, function(x){median(x,na.rm=T)})< -0.2)


df <- data.frame(group=c("Shorter","Longer"),value=c(as.numeric(tmp_short[2]),as.numeric(tmp_long[2])))

p6 <- ggplot(df,aes(x="",y=value,fill=group))+
  
  geom_bar(stat = 'identity',width=1, position = 'stack',fill=c("#BFDD85","#7CA1D6"))+
  
  coord_polar(theta = 'y')+theme_void()+
  
  geom_text(label=df$value)


(p1+p2+p3)/(p4+p5+p6)



table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[1:3], 1, function(x){median(x,na.rm=T)})> 0.2)
table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[1:3], 1, function(x){median(x,na.rm=T)})< -0.2)

table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[4:6], 1, function(x){median(x,na.rm=T)})> 0.2)
table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[4:6], 1, function(x){median(x,na.rm=T)})< -0.2)

table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[7:9], 1, function(x){median(x,na.rm=T)})> 0.2)
table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[7:9], 1, function(x){median(x,na.rm=T)})< -0.2)



table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[13:15], 1, function(x){median(x,na.rm=T)})> 0.2)
table(apply(pdui1[16:18], 1, function(x){median(x,na.rm=T)})- apply(pdui1[13:15], 1, function(x){median(x,na.rm=T)})< -0.2)



