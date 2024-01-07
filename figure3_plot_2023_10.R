
# PA活性分数----

tpm <- read_tsv("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/03-expression/tpm_matrix.txt") 


ensembl_ref_procoding <- read_tsv("~/cardiomyo_diff_APA/02-reference/ensembl.identifiers_2019") %>% 
  
  filter(`Gene type`=="protein_coding" )

tpm_filter <- tpm %>% 
  
  filter(gene %in% ensembl_ref_procoding$`Gene stable ID`)


cpsf_list <- read_tsv("~/cardiomyo_diff_APA/04-pdui/03-gene_list/GSEA_cpsf_pathway_M22041",col_names = F) %>% .$X1

ref_cpsf <- read_tsv("~/cardiomyo_diff_APA/02-reference/ensembl.identifiers_2019") %>% 
  
  filter(`Gene type`=="protein_coding" & `Transcript type`=="protein_coding" & `Gene name` %in% cpsf_list ) %>% 
  
  dplyr::select(1,4,5) %>% 
  
  unique() %>% .[-2,]


ttt <- tpm_filter %>% filter(gene %in% ref_cpsf$`Gene stable ID`)

tpm_scale <- data.frame(Geneid = tpm_filter$gene,scale(tpm_filter[,-1])) 


tpm_scale_cpsf_median <- filter(tpm_scale,Geneid %in% ref_cpsf$`Gene stable ID`) %>% .[-1] %>% 
  
  apply(2,median) 

tpm_scale_cpsf_median1 <- data.frame(sample=names(tpm_scale_cpsf_median),zscore_med=tpm_scale_cpsf_median) %>% mutate(time1=sub("\\.\\.\\..*","",sample),
                                                                                                                      time=sub("\\.","-",time1)
                                                                                                                      )

tpm_scale_cpsf_median1$time <- factor(tpm_scale_cpsf_median1$time,levels = unique(tpm_scale_cpsf_median1$time))


ggplot(tpm_scale_cpsf_median1,aes(x=time,y=zscore_med))+geom_boxplot(fill="#5CB4E5")+labs(y="APA machinery activity",x="Day")+theme_classic()+
  
  theme(axis.text.x = element_text( size=10,colour = "black"),
        axis.text.y = element_text(size=10,colour = "black"),
        title = element_text(size=11))+guides(fill=FALSE)


# RBP expression----



ref <- read_tsv("~/cardiomyo_diff_APA/02-reference/trend_factor_id.txt")

tpm <- read_tsv("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/03-expression/tpm_matrix.txt") %>% 
  
  filter(gene %in%  ref$`Gene stable ID`) %>% melt() 


  tpm$variable <- sub("\\.\\.\\..*","",tpm$variable)

tpm1 <- dcast(tpm, gene ~ variable,median) 

colnames(ref)[3] <- "gene"

gene_list <- ref[ref$`APA factor`== "key factor",] 

tpm2 <- as.data.frame(merge(ref,tpm1)) %>% filter(`Gene name` %in%  gene_list$`Gene name` )

tpm2 <- tpm2[,c(1:5,10,6:7,9,8)]

"BFU" = "BFU-E",
"CFU" = "CFU-E",
"proerythroblast" = "Pro-E",
"early_basophilic" = "Early baso-E",
"late_basophilic" = "Late baso-E",
"polychromatic" = "Poly-E",
"orthochromatic" = "Ortho-E"))


row.names(tpm2) <- tpm2$`Gene name`




tt <- t(scale(t(tpm2[-c(1:3)]))) 






pheatmap(tt,cluster_rows = T, cluster_cols = F,  number_color = "black", display_numbers = F,
         
         show_colnames = T,show_rownames = T, border_color = "white",
         
         col = colorRampPalette(c('blue', "white ",'red'))(200),
         #color = magma(100)[5*(1:20)],
         
         fontsize_number = 6,heatmap_legend_param = list(title = ""))





