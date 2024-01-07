setwd("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/01-data/hsc1_proe1_APA/pdui_res")

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

pdui4 <- pdui3[,c(1,7,8,2:6)] %>% .[-1] 


      
      pdui_matrix <- data.frame(matrix(ncol = 7,nrow = 7))
      
      colnames(pdui_matrix) <- stage_list
      
      row.names(pdui_matrix) <- stage_list

      



for (i in 1:(length(stage_list)-1)) {
  
  for (j in 1:(length(stage_list)-1)) {
  
  if (i>j) {
    
  }else{ 

  tmp_pdui <- pdui4[c(i,j+1)]
  
  colnames(tmp_pdui) <- c("ahead","fore")
  
  tmp_pdui1 <-  mutate(tmp_pdui,diff= fore-ahead,
                       result = case_when(diff >= 0.15 ~"Lengthen",
                                          diff <= -0.15 ~ "Shorten",
                                          TRUE ~ "Normal"
                                          )
                       )
  pdui_matrix[i,j+1]=table(tmp_pdui1$result)[3]
  
  pdui_matrix[j+1,i]= -table(tmp_pdui1$result)[1]
  #table(tmp_pdui1$result)[3]
  }
  }
}
      
      
      pdui_matrix1   <- data.frame( melt(data.frame(row = row.names(pdui_matrix),pdui_matrix)) )
      
      pdui_matrix1$variable  <- sub("\\.","-",pdui_matrix1$variable)
      
      pdui_matrix1$variable <- factor(pdui_matrix1$variable,levels = stage_list)
      
      pdui_matrix1$row <- factor(pdui_matrix1$row,levels = stage_list) 

      null_matrix <- data.frame(matrix(data= NA, nrow = 7,ncol = 7) %>% melt())
      

  ggplot() + 
  
    geom_point(data = pdui_matrix1, aes(x=variable, y=row,size=abs(value),fill=value),shape=21)+ 
  
    geom_tile(data = null_matrix, aes(x=Var1, y=Var2,size=abs(value)),color="grey",fill=NA,size=0.5) +
    
    
    scale_size(range = c(2,10), name="Gene Number") + theme_void()+ theme(
      
      axis.text  = element_text(size=14),  axis.text.x  = element_text(angle = 40)
      
    )+scale_fill_gradient2(low="blue", mid="white", high="red")

   

# 6 7

