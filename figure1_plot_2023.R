

# aUTR计算，UTR中位数计算 ----

pdui <- read_tsv("/home/yangyb/yangyb/Erythro_APA_Liu/merged_APA_20211208/02-result/01-pdui_and_motif_res/pdui_result.txt") %>% 
  
  separate(Gene,c("Transcript_ID","Gene_Name","3","strand"),sep = "\\|") %>% .[-c(2:5)]

three_utr_annotation <- read_tsv("/home/yangyb/software/DaPars2/3utr_anotation.txt") 

three_utr_annotation$Transcript_ID[duplicated(three_utr_annotation$Transcript_ID)]

pdui1 <- merge(three_utr_annotation,pdui)   %>% 

  separate(`3'UTR_Position`,c("chr","utr"),sep = ":") %>% 
  
  separate(utr,c("site1","site2"),sep = "-") %>% 
  
  separate(Loci,c("chr_loci","utr_loci"),sep = ":") %>% 
  
  separate(utr_loci,c("site1_loci","site2_loci"),sep = "-") %>% 
  
  filter(Predicted_Proximal_APA>=site1 & Predicted_Proximal_APA<=site2)


dup_name <- pdui1$Transcript_ID[duplicated(pdui1$Transcript_ID)]

pdui1 <- pdui1[!(pdui1$Transcript_ID %in%  dup_name),]


pdui1$site1 <- as.numeric(as.character(pdui1$site1))

pdui1$site2 <- as.numeric(as.character(pdui1$site2))

pdui1$site1_loci <- as.numeric(as.character(pdui1$site1_loci))

pdui1$site2_loci <- as.numeric(as.character(pdui1$site2_loci))




pdui2 <- mutate(pdui1,cutr=case_when(
  Strand=="+" ~ (Predicted_Proximal_APA - site1),
  Strand=="-" ~ (site2 - Predicted_Proximal_APA)
  
)) %>% 
  mutate(autr=case_when(
    Strand=="+" ~ (site2 - Predicted_Proximal_APA),
    Strand=="-" ~ (Predicted_Proximal_APA - site1)
    
  )) %>% 
  
  mutate(utr=abs(site1-site2)) %>% 
  
  mutate(utr1=abs(site1_loci-site2_loci)) %>% 
  
  mutate(utr_state=case_when(
    utr<500 ~ "<500 bp",
    utr>=500 & utr< 1000 ~ "500-1000 bp",
    utr>=1000 & utr< 1500 ~ "1000-1500 bp",
    utr>=1500 & utr< 2000 ~ "1500-2000 bp",
    utr>=2000 ~ ">=2000 bp"
  )) %>% filter(utr<=utr1) %>% filter(cutr>10) 



# 3'UTR和aUTR的相关性----

# w= 3.72 h=3.4

pdui2$utr_state=factor(pdui2$utr_state,levels = unique(pdui2$utr_state)[c(2,1,5,4,3)])


# 图1：aUTR不同长度的数量 ----

# 5 3.5

library(ggsci)

my_colors <- c("#E55858","#3C71A3",  "#66ABAF",  "#FE8E10", "#47964C", "#B47A9C","#FCC739")


 ggplot(pdui2) +
  geom_density(aes(x = autr, color = utr_state)) +
  xlim(0, 3000) +
  labs(x = "aUTR length (nt)", y = "Density",color = "3'UTR length") +
  theme_classic() +
  theme(axis.text.x = element_text(size=11, colour = "black"),
        axis.text.y = element_text(size=11, colour = "black"),
        title = element_text(size=12)) +
  scale_color_manual(values = my_colors)

# 图2：cUTR不同长度的数量 ----

 # 5 3.5
 
ggplot(pdui2) +
  geom_density(aes(x = cutr, color = utr_state)) +
  xlim(0, 3000) +
  labs(x = "cUTR length (nt)", y = "Density",color = "3'UTR length") +
  theme_classic() +
  theme(axis.text.x = element_text(size=11, colour = "black"),
        axis.text.y = element_text(size=11, colour = "black"),
        title = element_text(size=12)) +
   scale_color_manual(values = my_colors)


# 图3：3'UTR和aUTR的相关性 ----

# 4.5 4
 
ggplot(pdui2, aes(x = utr, y = autr)) + 
  geom_point(color = "#3C71A3") +  # change point color to blue
  xlab("3'UTR size of APA events (nt)") + 
  ylab("aUTR size of APA events (nt)") + 
  geom_smooth(method = "lm", color = "#FE8E10", alpha = 0.4) +  # change line color to red
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, colour = "black"),
    axis.text  = element_text(size = 10, colour = "black"), 
    strip.text = element_text(size = 12)
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~~")), size = 5)


# 图4：3'UTR和cUTR的相关性 ----
 
 # 4.5 4

ggplot(pdui2, aes(x = utr, y = cutr)) + 
  geom_point(color = "#3C71A3") +  # change point color to blue
  xlab("3'UTR size of APA events (nt)") + 
  ylab("cUTR size of APA events (nt)") + 
  geom_smooth(method = "lm",color = "#FE8E10", alpha = 0.4) +  # change line color to red
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, colour = "black"),
    axis.text  = element_text(size = 10, colour = "black"), 
    strip.text = element_text(size = 12)
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~~")), size = 5)



# aUTR和分析 ----


# 使用正则表达式获取与阶段相关的列名
stage_related_columns <- grep("^(SRR\\d+_GSM\\d+)_hs_(.+?)_\\d+_Homo_sapiens_RNA-Seq_PDUI$", colnames(pdui2), value = TRUE)

# 从这些列名中提取唯一的阶段标识
stages_extracted <- unique(sub("^(SRR\\d+_GSM\\d+)_hs_(.+?)_\\d+_Homo_sapiens_RNA-Seq_PDUI$", "\\2", stage_related_columns))

# 按照生物学意义对阶段进行排序
stages_ordered <- c("BFU", "CFU", "proerythroblast", "early_basophilic", "late_basophilic", "polychromatic", "orthochromatic")

stages_extracted <- intersect(stages_ordered, stages_extracted)


stage_columns_list <- lapply(stages_extracted, function(stage) {
  grep(paste0("_hs_", stage, "_\\d+_Homo_sapiens_RNA-Seq_PDUI$"), colnames(pdui2), value = TRUE)
})

all_stage_columns <- unlist(stage_columns_list)

new_pdui <- pdui2[, c("utr_state", all_stage_columns)]

long_pdui <- new_pdui %>%
  gather(key = "Stage", value = "APA_Value", all_stage_columns)

long_pdui$Simplified_Stage <- sub("^.*_hs_(.+?)_\\d+_Homo_sapiens_RNA-Seq_PDUI$", "\\1", long_pdui$Stage)


# 分析: Interaction between Stage and UTR State on APA Value ----


anova_result <- aov(APA_Value ~ Simplified_Stage * utr_state, data = long_pdui, na.action = na.omit)
summary(anova_result)

# Stage：该因子的影响是高度显著的（p < 2e-16），这意味着不同的Stage与APA值有显著的关系。

# utr_state：这个因子的影响也是高度显著的（p < 2e-16），这意味着不同的utr_state与APA值有显著的关系。

# Stage:utr_state（交互作用）：交互作用也是高度显著的（p < 2e-16）。这表示Stage和utr_state之间有显著的交互作用，即一个因子的影响取决于另一个因子的水平。
#
# 这表明mRNA的不同3'UTR长度类别可能在红系分化中有不同的作用。
#

# 图3: Interaction between Stage and UTR State on APA Value ----

unique(long_pdui$Simplified_Stage)

long_pdui <- long_pdui %>%
  mutate(Simplified_Stage = recode(Simplified_Stage,
                                   "BFU" = "BFU-E",
                                   "CFU" = "CFU-E",
                                   "proerythroblast" = "Pro-E",
                                   "early_basophilic" = "Early baso-E",
                                   "late_basophilic" = "Late baso-E",
                                   "polychromatic" = "Poly-E",
                                   "orthochromatic" = "Ortho-E"))

long_pdui$Simplified_Stage <- factor(long_pdui$Simplified_Stage, levels = unique(long_pdui$Simplified_Stage))

# 5 4

ggplot(long_pdui, aes(x = Simplified_Stage, y = APA_Value, color = utr_state, group = utr_state)) + ylim(0,1)+
  stat_summary(fun.data = mean_se, geom = "line", position = position_dodge(0.2)) + 
  stat_summary(fun.data = mean_se, geom = "point", position = position_dodge(0.2)) +
  theme_bw() +
  labs(
        title = "",x="",
       y = "Mean PDUI value",color = "3'UTR length") +
  theme_test() +
  theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(size=11, angle = 45, hjust = 1,colour = "black"),
        axis.text.y = element_text(size=11, colour = "black"),
        title = element_text(size=12)) +
  scale_color_manual(values = my_colors)


