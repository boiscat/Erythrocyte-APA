setwd("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/01-data/hsc1_proe1_APA/01-pdui_res")

library(reshape2)

pdui <- read_tsv("pdui_result.txt") %>% .[-c(2:4)]

pdui1 <- pdui  %>% melt()

anno <- data.frame(variable= colnames(pdui)[-1])

anno$state[1:3] <- "Pro-E"
anno$state[4:6] <- "Early-E"
anno$state[7:9] <- "Late-E"
anno$state[10:12] <- "Poly"
anno$state[13:15] <- "Ortho"
anno$state[16:18] <- "BFU-E"
anno$state[19:21] <- "CFU-E"

pdui2 <- merge(pdui1,anno) %>% .[-1] %>% na.omit()

pdui3 <- pivot_wider(pdui2,names_from = state,values_from = value,values_fn = median) 

stage_list <- unique(anno$state)[c(6,7,1:5)]

pdui4 <- pdui3[,c(1,7,8,2:6)]  %>% .[-1]

pdui3$Gene

lengthen_matrix <- data.frame(vector())

shorten_matrix <- data.frame(vector())

for (i in 1:(length(stage_list)-1)) {
  
  
      tmp_pdui <- pdui4[c(i,i+1)]
      
      colnames(tmp_pdui) <- c("ahead","fore")
      
      tmp_pdui1 <-  mutate(tmp_pdui,event_id = pdui3$Gene,
                                        diff = fore-ahead,
                                      result = case_when(diff >= 0.15 ~"Lengthen",
                                                         diff <= -0.15 ~ "Shorten",
                                                         TRUE ~ "Normal"
                                                         )) %>% 
        
                    separate(event_id,c("a","gene","c","strand"),sep = "\\|",remove=F) 
 
  # 切换event_id和gene
      
      lengthen <- data.frame( filter(tmp_pdui1,result=="Lengthen") %>% .$event_id )
      
      colnames(lengthen) <- paste0(stage_list[i],"_",stage_list[i+1])
      
      lengthen_matrix <- rbind.fill(lengthen_matrix,lengthen)
      
      
      
      shorten <- data.frame(filter(tmp_pdui1,result=="Shorten") %>% .$event_id)
      
      colnames(shorten) <- paste0(stage_list[i],"_",stage_list[i+1])
       
      shorten_matrix <-  rbind.fill(shorten_matrix,shorten)
    
    
  }




write.table(data.frame(t(lengthen_matrix[-1])),"lengthen_transcript_list.txt",quote = F,sep = ",",row.names = T,col.names = F)

write.table(data.frame(t(shorten_matrix[-1])),"shorten_transcript_list.txt",quote = F,sep = ",",row.names = T,col.names = F)

