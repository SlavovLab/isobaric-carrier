# Do we see increased signal and less missing data 
# in the reporter ion channels when we accumulate 
# more ions (using increased MS2 fill times) ?? 

# Source useful functions
source("code/functions_parameters.R")

# Load files: 
ev1<-read.delim("dat/hs237-240/evidence.txt")

# Filter peptides to <1% FDR and remove reverse and contaminant hits
ev1<-ev1[ev1$PEP<0.02, ]
ev1<-ev1[-which(ev1$Potential.contaminant=="+"), ]
ev1<-ev1[-which(ev1$Reverse=="+"), ]

# Take the mean reporter ion intensity in the single cell channels
ev1$scRI<-log2( rowMeans(ev1[,68:73]) )

# Calculate the T-cell to Monocyte peptide quantitation ratios from the two carrier channels
ev1$carJU<-log2(ev1$Reporter.intensity.1/ev1$Reporter.intensity.2)
ev1$carJU[abs(ev1$carJU)==Inf]<-NA
ev1$carJU[abs(ev1$carJU)==-Inf]<-NA

# Calculate the mean T-cell to Monocyte peptide quantitation ratios from the single cell channels
ev1$scJU<-log2( rowMeans(ev1[,seq(68,73,2)], na.rm=T) / rowMeans(ev1[,seq(69,73,2)], na.rm=T))
ev1$scJU[abs(ev1$scJU)==Inf]<-NA
ev1$carJU[abs(ev1$carJU)==-Inf]<-NA

# Remove duplicate PSMs per experiment
ev1$ms<-paste0(ev1$Modified.sequence,ev1$Charge)
ev1$mr<-paste0(ev1$ms,ev1$Raw.file)
ev1<-remove.duplicates(ev1, "mr")

# Only keep peptides in common to all 3 runs for a controlled comparison
kp<-names(table(ev1$ms))[table(ev1$ms)==length(unique(ev1$Raw.file))]

# Keep only peptides observed in all runs
ev1<-ev1[ev1$ms%in%kp, ]

# Fill times, recycled variable name, ignore that
carrier_labs<-c("100 ms","200 ms","300 ms","600 ms")

# Plot the scores
p2<-ggplot(ev1, aes(y=Raw.file, x=scRI, fill="red")) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_manual(values=rep(rgb(0.8,0,0), length(unique(ev1$Raw.file)))) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "green") + 
  coord_flip() +
  xlim(c(5,16)) +
  theme_bw() + 
  scale_y_discrete(labels= carrier_labs) + 
  ylab("\n Max fill time (ms)") + 
  ggtitle("") +
  rremove("legend") +
  xlab("scRI, log2\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

# Bin the highest intensities for easier visualization
ev1$Score[ev1$Score>250]<-250

p1<-ggplot(ev1, aes(y=Raw.file, x=Score, fill="red")) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_manual(values=rgb(0.8,0,0)) + 
  stat_summary(
    geom = "point",
    fun.y = "mean",
    col = "black",
    size = 3,
    shape = 24,
    fill = "green") + 
  coord_flip() +
  #xlim(c(5,16)) +
  theme_bw() + 
  scale_y_discrete(labels= carrier_labs) + 
  ylab("\n Max fill time (ms)") + 
  ggtitle("") +
  rremove("legend") +
  xlab("Andromeda Score\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))



# Look at the accuracy of quantitation

# Initialize variables
pear_cor<-c()
fts_x<-c()
pear_cor_CIup<-c()
pear_cor_CIdown<-c()
fts<-c(100,200,300, 600)

# Compute the correlation of the mean single cell J/U ratio to the bulk J/U ratio
# Compute 80% confidence intervals with resampling
set.seed(42)
for(X in unique(ev1$Raw.file)){
  
  # For each experiment
  evt<-ev1[ev1$Raw.file==X, ]
  
  pear_cor_x<-c()
  # Subsample and compute correlation 1000 times
  for(i in 1:1000){
    
    # Sample
    row_x<-sample(1:length(evt$carJU), size=length(evt$carJU)/10, replace = T)
    
    # Compute
    pear_cor_x<-c(pear_cor_x, cor(evt$carJU[row_x], evt$scJU[row_x], method = "spearman", use="complete.obs"))
    
  }
  
  # Take the mean correlation of the resampling
  pear_cor<-c(pear_cor, mean(pear_cor_x, na.rm = T))
  # Use the resampling to compute the confidence intervals
  pear_cor_CIup<-c(pear_cor_CIup, quantile(pear_cor_x, probs=0.1, na.rm = T))
  pear_cor_CIdown<-c(pear_cor_CIdown, quantile(pear_cor_x, probs=0.90, na.rm = T))
  
  # Record which experiment these values are computed from
  fts_x<-c(fts_x, fts[which(unique(ev1$Raw.file)==X)])
  
}

# Combine into data frame format
cors_df_ev2<-data.frame(fts_x, pear_cor, pear_cor_CIup, pear_cor_CIdown)
cors_df_ev2

p3<-ggplot(cors_df_ev2, aes(y=pear_cor, x=fts_x)) + 
  geom_point(size=4)+
  geom_errorbar(aes(ymin=pear_cor_CIup, ymax=pear_cor_CIdown), width=20) +
  theme_bw() + 
  #xlim(c(0,700))+
  ylim(c(0.5,0.9))+
  scale_x_continuous(limits=c(0, 700), breaks=c(0,100,200,300,400,500, 600,700), labels= c("","100 ms","200 ms","300 ms","","", "600 ms","")) +
  xlab("\n Max fill time (ms)") + 
  ylab("Correlation\n") + 
  #scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

p3b<-ggplot(cors_df_ev2, aes(y=pear_cor, x=as.character(fts_x))) + 
  geom_dotplot(y=pear_cor, binwidth = 1, dotsize = 0.1)+
  geom_errorbar(aes(ymin=pear_cor_CIup, ymax=pear_cor_CIdown), width=0.1) +
  theme_bw() + 
  #xlim(c(0,700))+
  ylim(c(0.5,0.9))+
  #scale_x_continuous(limits=c(0, 700), breaks=c(0,100,200,300,400,500, 600,700), labels= c("","100 ms","200 ms","300 ms","","", "600 ms","")) +
  xlab("\n Max fill time (ms)") + 
  ylab("Correlation\n") + 
  #scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))



# Show missing data improvement

pctNA<-function(x){ 
  
  return( ( length(which(is.na(x))) + length(which(x==0)) ) / length(x))
  
}

evm<-melt(ev1[,c(which(colnames(ev1)%in%c("Raw.file")), 68:73)  ]) %>% group_by(Raw.file, variable) %>% mutate(pct = pctNA(value))

evm<-remove.duplicates(evm, c("Raw.file","variable","pct"))

evm$pct<-evm$pct*100

evm$ft<-c(100,200,300,600)[as.numeric(factor(evm$Raw.file))]

p4<-ggplot(evm, aes(x=ft, y=pct, color=variable)) +
geom_jitter(width = 30, size=2) +
  #geom_dotplot(binaxis = "y", stackdir = "center", dotsize=2, position = position_dodge(width=0.5)) +
  theme_bw() + 
  scale_y_continuous(labels = comma) + 
  #scale_x_discrete(breaks=1:6,labels=c("100 ms","200 ms","300 ms","","", "600 ms")) +
  scale_x_continuous(limits=c(0, 700), breaks=c(0,100,200,300,400,500, 600,700), labels= c("","100 ms","200 ms","300 ms","","", "600 ms","")) +
  ylab("% missing data\n") + 
  #rremove("legend") +
  xlab("\n Max fill time (ms)") + 
  ggtitle("")+
  theme(legend.position="right") + 
  theme(legend.title = element_blank()) + 
  scale_color_discrete(name = "Samples", labels = c("1 HEK-293", "1 U-937", "1 HEK-293", "1 U-937","1 HEK-293", "1 U-937")) + 
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.x = element_text(color = "grey20", size = 18, angle = 00, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.ticks.x = element_blank() )

p1 + p3 + p2 + p4 +plot_layout(nrow=2)
