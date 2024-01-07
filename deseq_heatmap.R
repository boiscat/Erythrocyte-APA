# Loading required libraries
library(ggplot2)
library(reshape2)

# Data Frame
data <- data.frame(
  Description = c("VEGFA-VEGFR2 signaling", "Neutrophil degranulation", "regulation of transferase activity",
                  "intracellular protein transport", "Cytokine Signaling in Immune system", "Signaling by Rho GTPases",
                  "Transcriptional Regulation by TP53", "cell activation", "Hemostasis", "Translation",
                  "ribonucleoprotein complex biogenesis", "S Phase",
                  "mitochondrion organization", "Cell Cycle", "Processing of Capped Intron-Containing Pre-mRNA",
                  "DNA metabolic process", "protein-RNA complex assembly", "mitochondrial gene expression"),
  Down_regulated = c(-9.591063819, -9.764965132, -5.808135128, -11.24622934, -3.9361338, -3.297540971,
                     -17.95409526, 0, 0, -100, -66.78753282, -23.7441405,
                     -27.80356744, -34.95988576, -29.93063, -21.12419469, -28.51186463, -27.65781544),
  Up_regulated = c(-12.19276297, -14.71078104, -11.87286064, -7.396480892, -26.91253019, -18.34996512,
                   -3.74234747, -28.02756519, -25.80036744, 0, 0, 0,
                   0, 0, 0, 0, 0, 0)
)

# Cap values at 35 and replace 0 with NA
data$Down_regulated[data$Down_regulated > 35] <- 35
data$Up_regulated[data$Up_regulated > 35] <- 35
data$Down_regulated[data$Down_regulated == 0] <- NA
data$Up_regulated[data$Up_regulated == 0] <- NA

# Melting the data for use with ggplot
melted_data <- melt(data, id.vars = 'Description')

# Plotting the heatmap
ggplot(melted_data, aes(variable, Description, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("#FFFFFF", "#FF0000"), na.value = "#FFFFFF", limits = c(NA, 35)) +
  theme_minimal() +
  theme(axis.text.x = element_text( hjust = 1),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(fill = "-LogP")
