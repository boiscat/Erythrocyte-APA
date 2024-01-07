

pdui <- read_tsv("~/yangyb/Erythro_APA_Liu/merged_APA_20211208/01-data/hsc1_proe1_APA/pdui_res/pdui_result.txt")

bed_file <- read_tsv("/home/yangyb/cardiomyo_diff_APA/02-reference/hg19_refseq_extracted_3UTR_2019.bed",col_names = F) %>%
  
  mutate(utr_l= X3-X2) %>%
  
  dplyr::filter(X1!="chrX" & X1!="chrY" & X1!="chrM")

multi_APA <- bed_file[bed_file$X4 %in% pdui$Gene,] %>%
  
  mutate(type="multi_APA")

single_APA <- bed_file[!(bed_file$X4 %in% pdui$Gene),] %>%
  
  mutate(type="single_APA")

bed_res <- full_join(multi_APA,single_APA)

ggplot(bed_res,aes(x=type,y=utr_l))+geom_boxplot(width=0.4,fill="grey90")+ylim(c(1,3000))+ labs(x="",y="3'UTR Length (nt)")+
  
  theme(  axis.text.x = element_text(size=11,colour = "black"),
          axis.text.y = element_text(hjust = 1,size=11,colour = "black"),
          title = element_text(size=12),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.8))+
  
  scale_fill_tron()+ stat_compare_means()

compare_means(utr_l~type,bed_res)


