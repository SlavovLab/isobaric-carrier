# Does another search engine, besides MaxQuant, also show the same trend in increasing
# quality of PSM as the amount of cells in the carrier increases? 

# Looking at SpectroMine (Spectronaut / Biognosys) search results on the same files

# Useful functions
source("code/functions_parameters.R")

# Load the output of the SpectroMine search
df<-read.delim("ic16.xls")

# Remove duplicate PSMs per experiment
df$mr<-paste0(df$PEP.StrippedSequence,df$R.FileName)
df<-remove.duplicates(df, "mr")

# Only keep peptides in common to all 10 runs for a controlled comparison
kp<-names(table(df$PEP.StrippedSequence))[table(df$PEP.StrippedSequence)==10]

# Fix the autosampler mixup (200 and 300 cell carrier sampels swapped)
df$rf<-df$R.FileName
df$rf[df$R.FileName=="hs0042.raw"]<-"hs0043.raw"
df$rf[df$R.FileName=="hs0043.raw"]<-"hs0042.raw"
df$rf[df$R.FileName=="hs0101.raw"]<-"hs0102.raw"
df$rf[df$R.FileName=="hs0102.raw"]<-"hs0101.raw"

# Bin all high scores for easier visual comparison
df$PSM.Score[df$PSM.Score>350]<-350

# Keep only peptides observed in all runs
df2<-df[df$PEP.StrippedSequence%in%kp, ]

# Labels
carrier_labs<-c("100","200","300","400","600","800", "100","200","300","400","600","800")

# Plot the scores
ggplot(df2, aes(y=rf, x=PSM.Score)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "green") + 
  coord_flip() +
  theme_bw() + 
  scale_y_discrete(labels= carrier_labs) + 
  ylab("\nNumber of cells in isobaric carrier") + 
  xlab("Spectronaut Score\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))



