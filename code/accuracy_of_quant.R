# How does the number of cells in the isobaric carrier effect
# the single cell quantitation ?

# Does filling the orbitrap with more ions change the effect ?

# Is the trend dependant on the single cell reporter ion intensity ? 


# Source useful functions
source("code/functions_parameters.R")

########################################################################
##### Low AGC series
####
##
#

# Load files: 
exps<-c(paste0("dat/carrier_202004xx/hs00",41:46,"/evidence.txt") )
ev1<-read.delim(exps[1])

for(X in exps[-c(1)]){
  
  ev1<-rbind(ev1, read.delim(X))
  
}

# Correct autosampler mixup of 200 and 300-cell carrier samples
ev1$rf<-ev1$Raw.file
ev1$rf[ev1$Raw.file=="hs0042"]<-"hs0043"
ev1$rf[ev1$Raw.file=="hs0043"]<-"hs0042"
ev1$Raw.file<-as.factor(ev1$rf)

# Filter peptides to <1% FDR and remove reverse and contaminant hits
ev1<-ev1[ev1$PEP<0.02, ]
ev1<-ev1[-which(ev1$Potential.contaminant=="+"), ]
ev1<-ev1[-which(ev1$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev1$ms<-paste0( ev1$Modified.sequence, ev1$Charge )

# Remove duplicate peptides
ev1<-remove.duplicates(ev1, c("ms", "Raw.file"))

# Take only peptides observed in all 6 experiments
table_peps<-table(ev1$ms)
intersect_peps<-names(table_peps)[table_peps==length(unique(ev1$Raw.file))]

# Calculate the T-cell to Monocyte peptide quantitation ratios from the two carrier channels
ev1$carJU<-log2(ev1$Reporter.intensity.1/ev1$Reporter.intensity.2)
ev1$carJU[abs(ev1$carJU)==Inf]<-NA

# Calculate the mean T-cell to Monocyte peptide quantitation ratios from the single cell channels
ev1$scJU<-log2(rowMeans(ev1[,seq(72,79,2)], na.rm=T) / rowMeans(ev1[,seq(73,79,2)], na.rm=T))
ev1$scJU[abs(ev1$scJU)==Inf]<-NA

# For visual purposes only: Bin the low intensity single cell reporter ion values
ev1[, 72:79][ev1[, 72:79]<=100]<-100

# Calculate the mean single cell signal across both cell types
ev1$meanSCRI<-log10( rowMeans( ev1[, 72:79], na.rm=T) )
ev1$meanSCRI[abs(ev1$meanSCRI)==Inf]<-NA

# Labels for later: 
carrier_amounts<-c(100, 200, 300, 400, 600, 800)

# For a controlled comparison, consider only peptides observed in every experiment:
ev2<-ev1[ev1$ms%in%intersect_peps, ]

# Initialize variables
pear_cor<-c()
carrier_amount_x<-c()
pear_cor_CIup<-c()
pear_cor_CIdown<-c()

# Compute the correlation of the mean single cell J/U ratio to the bulk J/U ratio
# Compute 80% confidence intervals with resampling
for(X in levels(ev2$Raw.file)){
  
  # For each experiment
  evt<-ev2[ev2$Raw.file==X, ]
  
  set.seed(42)
  pear_cor_x<-c()
  
  # Subsample and compute correlation 1000 times
  for(i in 1:1000){
    
    # Sample
    row_x<-sample(1:length(evt$carJU), size=length(evt$carJU)/10, replace = T)
    
    # Compute
    pear_cor_x<-c(pear_cor_x, cor(evt$carJU[row_x], evt$scJU[row_x], method = "pearson", use="complete.obs"))
    
  }
  
  # Take the mean correlation of the resampling
  pear_cor<-c(pear_cor, mean(pear_cor_x, na.rm = T))
  # Use the resampling to compute the confidence intervals
  pear_cor_CIup<-c(pear_cor_CIup, quantile(pear_cor_x, probs=0.1, na.rm = T))
  pear_cor_CIdown<-c(pear_cor_CIdown, quantile(pear_cor_x, probs=0.90, na.rm = T))
  
  # Record which experiment these values are computed from
  carrier_amount_x<-c(carrier_amount_x, carrier_amounts[which(levels(ev2$Raw.file)==X)])
  
}

# Combine into data frame format
cors_df_ev2<-data.frame(carrier_amount_x, pear_cor, pear_cor_CIup, pear_cor_CIdown)
cors_df_ev2$carrier_amount_x<-factor(as.character(cors_df_ev2$carrier_amount_x), levels = unique(as.character(cors_df_ev2$carrier_amount_x) ) )

# Record trimmed data frame as another variable for later
interp.df<-data.frame(1:length(cors_df_ev2$pear_cor), cors_df_ev2$pear_cor, cors_df_ev2$pear_cor_CIup,cors_df_ev2$pear_cor_CIdown); colnames(interp.df)<-c("x","y","ci_u","ci_d")
interp.df1<-interp.df

# ???
ev2_1<-ev2

# Plot the mean single cell reporter ion intensity as a function of number of cells in the 
# isobaric carrier
p1<-ggplot(ev2, aes(y=Raw.file, x=meanSCRI)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  coord_flip() +
  theme_bw() + 
  scale_y_discrete(labels= carrier_amounts) + 
  ylab("") + 
  ggtitle("Reporter ion intensities in single cell channels") +
  xlab(expression("RI,"~log["10"])) + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

# Plot the Andromeda score as a function of number of cells in the 
# isobaric carrier
p2<-ggplot(ev2, aes(y=Raw.file, x=Score)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  coord_flip() +
  theme_bw() + 
  scale_y_discrete(labels= carrier_amounts) + 
  ylab("") + 
  xlab("Andromeda score\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

# Plot the accuracy of quantitation as a function of number of cells in the 
# isobaric carrier
p3<-ggplot(cors_df_ev2, aes(y=pear_cor, x=carrier_amount_x)) +  
  geom_point(size=4)+
  theme_bw() + 
  scale_x_discrete(labels= carrier_amounts) +
  ggtitle("Pearson correlations of single cell ratios\n to bulk ratio, (Jurkat / monocyte)") + 
  xlab("\nCarrier cell amount") + 
  ylab("Correlation\n") + 
  scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))


# Now, recompute the correlation values with resample, just as before, but for only a subset of the 
# peptides: those with single cell reporter ion intensity >2000. Since this results in skewed distributions
# for some of the samples, try to control for the distribution by binning across the range of reporter 
# ions >2000 and only take one peptide per bin, then compute correlation: 

# For curiousity, take a look at a range of single cell reporter ions, not just the high bin (>2000)
scRIq<-quantile(ev2$meanSCRI, na.rm=T)
scRIq[1]<-log10(800)
scRIq[2]<-log10(850)
scRIq[3]<-log10(950)
scRIq[4]<-log10(2000)
scRIq[5]<-5.7

# Number of bins
bx<-10

# Initiate variables
idf.c<-data.frame(x<-c(),  y<-c(), ci_u<-c(),ci_d<-c(), quant<-c() )

# Loop over the ranges of scRI above
for(j in 1:4){
  
  # Grab the peptides that fall into the given range of scRI
  ev2t<-ev2[ (ev2$meanSCRI > scRIq[j] )&(ev2$meanSCRI <= scRIq[j+1]),  ]
  
  # Compute the bin intervals
  bins20<-seq(scRIq[j], scRIq[j+1], length.out = bx+1)
  
  # Initiate variables
  spear_cor<-c()
  pear_cor<-c()
  pear_cor_sd<-c()
  pear_cor_CIup<-c()
  pear_cor_CIdown<-c()
  carrier_amount_x<-c()
  
  # For each experiment, compute the accuracy of correlation. This code is the same as above, for the 
  # set of peptides in common between all experiments
  for(X in levels(ev2t$Raw.file)){
    
    evt<-ev2t[ev2t$Raw.file==X, ]

    spear_cor<-c(spear_cor, cor(evt$carJU, evt$scJU, method = "spearman", use="complete.obs"))

    set.seed(42)
    pear_cor_x<-c()
    for(i in 1:1000){
      
      rk<-c()
      for(k in 1:bx){
        if(length(which(evt$meanSCRI > bins20[k] & evt$meanSCRI <= bins20[k+1])) > 0){
          rkt<-sample(which(evt$meanSCRI > bins20[k] & evt$meanSCRI <= bins20[k+1]),size = 1, replace = F)
          rk<-c(rk,rkt)
        }
      }
      
      pear_cor_x<-c(pear_cor_x, cor(evt$carJU[rk], evt$scJU[rk], method = "pearson", use="complete.obs"))
      
    }
    
    pear_cor<-c(pear_cor, median(pear_cor_x, na.rm = T))
    pear_cor_CIup<-c(pear_cor_CIup, quantile(pear_cor_x, probs=0.1, na.rm = T))
    pear_cor_CIdown<-c(pear_cor_CIdown, quantile(pear_cor_x, probs=0.90, na.rm = T))
    carrier_amount_x<-c(carrier_amount_x, carrier_amounts[which(levels(ev2$Raw.file)==X)])
    
  }
  
  # Record results
  cors_df_ev2<-data.frame(carrier_amount_x, pear_cor, pear_cor_CIup, pear_cor_CIdown)
  cors_df_ev2$carrier_amount_x<-factor(as.character(cors_df_ev2$carrier_amount_x), levels = unique(as.character(cors_df_ev2$carrier_amount_x) ) )
  cors_df_ev2$quant<-names(scRIq)[j+1]
  interp.dfx<-data.frame(1:length(cors_df_ev2$pear_cor), cors_df_ev2$pear_cor, cors_df_ev2$pear_cor_CIup,cors_df_ev2$pear_cor_CIdown,  cors_df_ev2$quant); colnames(interp.dfx)<-c("x","y","ci_u","ci_d", "quant")
  idf.c<-rbind( idf.c, interp.dfx)
  
}

# Rename variable holding results, will be written over... for convenience.
idf.c1<-idf.c

# Alternative plot of the correlations 
p3<-ggplot(interp.df, aes(y=y, x=x)) + 
  geom_point(size=4)+
  geom_errorbar(aes(ymin=ci_d, ymax=ci_u), width=.2) +
  theme_bw() + 
  scale_x_continuous(breaks=interp.df$x, labels= carrier_amounts) +
  ggtitle("Pearson correlations of single cell ratios\n to bulk ratio, (Jurkat / monocyte)") + 
  xlab("\nCarrier cell amount") + 
  ylab("Correlation\n") + 
  scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  geom_smooth(method="loess", n=10, span=10, se=F)

# Plot of the correlation using peptides with scRI >2000
p4<-ggplot(idf.c[idf.c$quant=="100%", ], aes(y=y, x=x)) + 
  geom_point(size=4)+
  geom_errorbar(aes(ymin=ci_d, ymax=ci_u), width=.2) +
  theme_bw() + 
  scale_x_continuous(breaks=interp.df$x, labels= carrier_amounts) +
  ggtitle("Correlations calculated using only high-abundance peptides") + 
  xlab("\nCarrier cell amount") + 
  ylab("Correlation\n") + 
  scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  geom_smooth(method="loess", n=10, span=10, se=F)

# Visualize
#p1 + p2 + p3 + p4 +plot_layout(ncol=1, nrow=4)



########################################################################
##### High AGC series
####
##
#

# Note: the following code is the same as the first section, just on different
# experimental data. Variables are overwritten from the first section

# Load files: 
exps<-c(paste0("dat/carrier_202004xx/hs0",100:105,"/evidence.txt") )
ev1<-read.delim(exps[1])

for(X in exps[-c(1)]){
  
  ev1<-rbind(ev1, read.delim(X))
  
}

# Correct autosampler mixup of 200 and 300-cell carrier samples
ev1$rf<-ev1$Raw.file
ev1$rf[ev1$Raw.file=="hs0042"]<-"hs0043"
ev1$rf[ev1$Raw.file=="hs0043"]<-"hs0042"
ev1$Raw.file<-as.factor(ev1$rf)

# Filter peptides to <1% FDR and remove reverse and contaminant hits
ev1<-ev1[ev1$PEP<0.02, ]
ev1<-ev1[-which(ev1$Potential.contaminant=="+"), ]
ev1<-ev1[-which(ev1$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev1$ms<-paste0( ev1$Modified.sequence, ev1$Charge )

# Remove duplicate peptides
ev1<-remove.duplicates(ev1, c("ms", "Raw.file"))

# Take only peptides observed in all 6 experiments
table_peps<-table(ev1$ms)
intersect_peps<-names(table_peps)[table_peps==length(unique(ev1$Raw.file))]

# Calculate the T-cell to Monocyte peptide quantitation ratios from the two carrier channels
ev1$carJU<-log2(ev1$Reporter.intensity.1/ev1$Reporter.intensity.2)
ev1$carJU[abs(ev1$carJU)==Inf]<-NA

# Calculate the mean T-cell to Monocyte peptide quantitation ratios from the single cell channels
ev1$scJU<-log2(rowMeans(ev1[,seq(72,79,2)], na.rm=T) / rowMeans(ev1[,seq(73,79,2)], na.rm=T))
ev1$scJU[abs(ev1$scJU)==Inf]<-NA

# For visual purposes only: Bin the low intensity single cell reporter ion values
ev1[, 72:79][ev1[, 72:79]<=100]<-100

# Calculate the mean single cell signal across both cell types
ev1$meanSCRI<-log10( rowMeans( ev1[, 72:79], na.rm=T) )
ev1$meanSCRI[abs(ev1$meanSCRI)==Inf]<-NA

# Labels for later: 
carrier_amounts<-c(100, 200, 300, 400, 600, 800)

# For a controlled comparison, consider only peptides observed in every experiment:
ev2<-ev1[ev1$ms%in%intersect_peps, ]

# Initialize variables
pear_cor<-c()
carrier_amount_x<-c()
pear_cor_CIup<-c()
pear_cor_CIdown<-c()

# Compute the correlation of the mean single cell J/U ratio to the bulk J/U ratio
# Compute 80% confidence intervals with resampling
for(X in levels(ev2$Raw.file)){
  
  # For each experiment
  evt<-ev2[ev2$Raw.file==X, ]
  
  set.seed(42)
  pear_cor_x<-c()
  
  # Subsample and compute correlation 1000 times
  for(i in 1:1000){
    
    # Sample
    row_x<-sample(1:length(evt$carJU), size=length(evt$carJU)/10, replace = T)
    
    # Compute
    pear_cor_x<-c(pear_cor_x, cor(evt$carJU[row_x], evt$scJU[row_x], method = "pearson", use="complete.obs"))
    
  }
  
  # Take the mean correlation of the resampling
  pear_cor<-c(pear_cor, mean(pear_cor_x, na.rm = T))
  # Use the resampling to compute the confidence intervals
  pear_cor_CIup<-c(pear_cor_CIup, quantile(pear_cor_x, probs=0.1, na.rm = T))
  pear_cor_CIdown<-c(pear_cor_CIdown, quantile(pear_cor_x, probs=0.90, na.rm = T))
  
  # Record which experiment these values are computed from
  carrier_amount_x<-c(carrier_amount_x, carrier_amounts[which(levels(ev2$Raw.file)==X)])
  
}

# Combine into data frame format
cors_df_ev2<-data.frame(carrier_amount_x, pear_cor, pear_cor_CIup, pear_cor_CIdown)
cors_df_ev2$carrier_amount_x<-factor(as.character(cors_df_ev2$carrier_amount_x), levels = unique(as.character(cors_df_ev2$carrier_amount_x) ) )

# Record trimmed data frame as another variable for later
interp.df<-data.frame(1:length(cors_df_ev2$pear_cor), cors_df_ev2$pear_cor, cors_df_ev2$pear_cor_CIup,cors_df_ev2$pear_cor_CIdown); colnames(interp.df)<-c("x","y","ci_u","ci_d")
interp.df1<-interp.df

# ???
ev2_1<-ev2

# Plot the mean single cell reporter ion intensity as a function of number of cells in the 
# isobaric carrier
p1<-ggplot(ev2, aes(y=Raw.file, x=meanSCRI)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  coord_flip() +
  theme_bw() + 
  scale_y_discrete(labels= carrier_amounts) + 
  ylab("") + 
  ggtitle("Reporter ion intensities in single cell channels") +
  xlab(expression("RI,"~log["10"])) + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

# Plot the Andromeda score as a function of number of cells in the 
# isobaric carrier
p2<-ggplot(ev2, aes(y=Raw.file, x=Score)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  coord_flip() +
  theme_bw() + 
  scale_y_discrete(labels= carrier_amounts) + 
  ylab("") + 
  xlab("Andromeda score\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

# Plot the accuracy of quantitation as a function of number of cells in the 
# isobaric carrier
p3<-ggplot(cors_df_ev2, aes(y=pear_cor, x=carrier_amount_x)) +  
  geom_point(size=4)+
  theme_bw() + 
  scale_x_discrete(labels= carrier_amounts) +
  ggtitle("Pearson correlations of single cell ratios\n to bulk ratio, (Jurkat / monocyte)") + 
  xlab("\nCarrier cell amount") + 
  ylab("Correlation\n") + 
  scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))


# Now, recompute the correlation values with resample, just as before, but for only a subset of the 
# peptides: those with single cell reporter ion intensity >2000. Since this results in skewed distributions
# for some of the samples, try to control for the distribution by binning across the range of reporter 
# ions >2000 and only take one peptide per bin, then compute correlation: 

# For curiousity, take a look at a range of single cell reporter ions, not just the high bin (>2000)
scRIq<-quantile(ev2$meanSCRI, na.rm=T)
scRIq[1]<-log10(800)
scRIq[2]<-log10(850)
scRIq[3]<-log10(950)
scRIq[4]<-log10(2000)
scRIq[5]<-5.7

# Number of bins
bx<-10

# Initiate variables
idf.c<-data.frame(x<-c(),  y<-c(), ci_u<-c(),ci_d<-c(), quant<-c() )

# Loop over the ranges of scRI above
for(j in 1:4){
  
  # Grab the peptides that fall into the given range of scRI
  ev2t<-ev2[ (ev2$meanSCRI > scRIq[j] )&(ev2$meanSCRI <= scRIq[j+1]),  ]
  
  # Compute the bin intervals
  bins20<-seq(scRIq[j], scRIq[j+1], length.out = bx+1)
  
  # Initiate variables
  spear_cor<-c()
  pear_cor<-c()
  pear_cor_sd<-c()
  pear_cor_CIup<-c()
  pear_cor_CIdown<-c()
  carrier_amount_x<-c()
  
  # For each experiment, compute the accuracy of correlation. This code is the same as above, for the 
  # set of peptides in common between all experiments
  for(X in levels(ev2t$Raw.file)){
    
    evt<-ev2t[ev2t$Raw.file==X, ]
    
    spear_cor<-c(spear_cor, cor(evt$carJU, evt$scJU, method = "spearman", use="complete.obs"))
    
    set.seed(42)
    pear_cor_x<-c()
    for(i in 1:1000){
      
      rk<-c()
      for(k in 1:bx){
        if(length(which(evt$meanSCRI > bins20[k] & evt$meanSCRI <= bins20[k+1])) > 0){
          rkt<-sample(which(evt$meanSCRI > bins20[k] & evt$meanSCRI <= bins20[k+1]),size = 1, replace = F)
          rk<-c(rk,rkt)
        }
      }
      
      pear_cor_x<-c(pear_cor_x, cor(evt$carJU[rk], evt$scJU[rk], method = "pearson", use="complete.obs"))
      
    }
    
    pear_cor<-c(pear_cor, median(pear_cor_x, na.rm = T))
    pear_cor_CIup<-c(pear_cor_CIup, quantile(pear_cor_x, probs=0.1, na.rm = T))
    pear_cor_CIdown<-c(pear_cor_CIdown, quantile(pear_cor_x, probs=0.90, na.rm = T))
    carrier_amount_x<-c(carrier_amount_x, carrier_amounts[which(levels(ev2$Raw.file)==X)])
    
  }
  
  # Record results
  cors_df_ev2<-data.frame(carrier_amount_x, pear_cor, pear_cor_CIup, pear_cor_CIdown)
  cors_df_ev2$carrier_amount_x<-factor(as.character(cors_df_ev2$carrier_amount_x), levels = unique(as.character(cors_df_ev2$carrier_amount_x) ) )
  cors_df_ev2$quant<-names(scRIq)[j+1]
  interp.dfx<-data.frame(1:length(cors_df_ev2$pear_cor), cors_df_ev2$pear_cor, cors_df_ev2$pear_cor_CIup,cors_df_ev2$pear_cor_CIdown,  cors_df_ev2$quant); colnames(interp.dfx)<-c("x","y","ci_u","ci_d", "quant")
  idf.c<-rbind( idf.c, interp.dfx)
  
}

# Alternative plot of the correlations 
p3<-ggplot(interp.df, aes(y=y, x=x)) + 
  geom_point(size=4)+
  geom_errorbar(aes(ymin=ci_d, ymax=ci_u), width=.2) +
  theme_bw() + 
  scale_x_continuous(breaks=interp.df$x, labels= carrier_amounts) +
  ggtitle("Pearson correlations of single cell ratios\n to bulk ratio, (Jurkat / monocyte)") + 
  xlab("\nCarrier cell amount") + 
  ylab("Correlation\n") + 
  scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  geom_smooth(method="loess", n=10, span=10, se=F)

# Plot of the correlation using peptides with scRI >2000
p4<-ggplot(idf.c[idf.c$quant=="100%", ], aes(y=y, x=x)) + 
  geom_point(size=4)+
  geom_errorbar(aes(ymin=ci_d, ymax=ci_u), width=.2) +
  theme_bw() + 
  scale_x_continuous(breaks=interp.df$x, labels= carrier_amounts) +
  ggtitle("Correlations calculated using only high-abundance peptides") + 
  xlab("\nCarrier cell amount") + 
  ylab("Correlation\n") + 
  scale_y_continuous(limits=c(0,1))+ 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  geom_smooth(method="loess", n=10, span=10, se=F)

# Visualize
p1 + p2 + p3 + p4 +plot_layout(ncol=1, nrow=4)


