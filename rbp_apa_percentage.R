

pdui <- read_tsv("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/01-pdui_and_motif_res/pdui_result.txt") %>% .[c(1,20:25,5:19)]

tt <- read_tsv("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/script/02-diff_0.15_for_mfuzz_sorted.txt")

tt <-  apply(pdui[-1],1,function(x){max(x,na.rm = T)-min(x,na.rm = T)>=0.15})

pdui_diff <- pdui[pdui$Gene %in% tt$Gene,] 

tpm <- read_tsv("~/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/03-expression/tpm_matrix.txt")

ref <- read_tsv("~/cardiomyo_diff_APA/02-reference/trend_factor_id.txt") %>% 
  
       filter(`APA factor`=="key factor") %>% .[-2]

colnames(ref) <- c("gene_name","gene")

tpm1 <- merge(ref,tpm)



rbp_name <- vector()

apa_name <- vector()

p_v <- vector()

r <- vector()

 
    for (rbp_num in 1:dim(tpm1)) {
      
                                          print(rbp_name)
        
        for (apa_num in 1:dim(pdui_diff)) {
                      
                                          rbp_name <- append(rbp_name,tpm1[rbp_num,2])
                                          
                                          apa_name <- append(apa_name,as.character(pdui_diff[apa_num,1]))
              
                                          tmp  <- cor.test( as.numeric(tpm1[rbp_num,-c(1:2)]),as.numeric(pdui_diff[apa_num,-1]),method = "pearson")
                                        
                                          p_v <- append(p_v,tmp[["p.value"]])
                                          
                                          r <- append(r,tmp[["estimate"]][["cor"]])
                
                                }
    }


res <- data.frame(rbp_name=rbp_name,apa_name=apa_name, r=r,p_v=p_v,fdr=p.adjust(p_v,method = "fdr")) %>% 
  
                  filter(fdr<0.05,r^2>0.3) %>% 
  
                  group_by(rbp_name) %>% 
  
                  dplyr::summarise(percentage=n()/dim(pdui_diff)[1]) %>% 
  
                  arrange(percentage)

res$rbp_name <- factor(res$rbp_name,levels = rev(res$rbp_name))

library(ggsci)
ggplot(res,aes(x=rbp_name,y=percentage,fill=rbp_name))+geom_bar(stat= "identity",width = 0.7) +
  
        theme_classic()+ xlab("")+ ylab("Related APA events (%)")+
  
        scale_y_continuous(expand=c(0,0))+
  
        scale_fill_igv()+ 
  
        theme(
          axis.text.x = element_text(size=10,colour = "black",angle =40,vjust = 0.6),
          axis.text  = element_text(size=10,colour = "black"),
          axis.title   = element_text(size=12,colour = "black")
            
              )+guides(fill=F)
        
# 3 w=6

res1 <- data.frame(rbp_name=rbp_name,apa_name=apa_name, r=r,p_v=p_v,fdr=p.adjust(p_v,method = "fdr")) %>% 
  
        filter(fdr<0.05,r^2>0.3) %>% 
        
        filter(rbp_name=="PABPC1")


lengthen <- tt %>% mutate(diff=CFU-BFU) %>% filter(diff>=0.15 ) %>% filter(Gene %in% res1$apa_name)

shorten <- tt %>% mutate(diff=CFU-BFU) %>% filter(diff<= -0.15) %>% filter(Gene %in% res1$apa_name)

