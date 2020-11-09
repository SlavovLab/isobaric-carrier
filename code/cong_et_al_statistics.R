# All MaxQuant output referenced below was downloaded from the Cong et al data repository
# at ProteomeXchagne: PXD016921 

# Source functions
source("code/functions_parameters.R")

# Read in MaxQuant output
msms<-read.delim("dat/txt/msmsScans.txt") 

# Count the number of MS2 scans taken
df_ms<-aggregate( Scan.index ~ Raw.file, data = msms, FUN=max)

# Read in MaxQuant output
ev<-read.delim("dat/txt/evidence.txt") 

# Calculate FDR on a per experiment basis
for(X in unique(ev$Raw.file)){
  ev$fdr[ev$Raw.file==X]<-calc_fdr(ev$PEP[ev$Raw.file==X])
}

# Filter peptides to 1% FDR
ev<-ev[ev$fdr<0.01, ]
ev$index<-1

# Count the number of PSMs at 1% FDR 
df_psm<-aggregate( index ~ Raw.file, data = ev, FUN=sum)

# Read in MaxQuant output
ap<-read.delim("dat/txt/allPeptides.txt") 

# Count the number of MS1 features identified by MaxQuant
df_ap<-data.frame(df_ms$Raw.file, 
                  c(length(which(ap$Charge[ap$Raw.file==unique(ap$Raw.file)[2]]>=2)), 
                    length(which(ap$Charge[ap$Raw.file==unique(ap$Raw.file)[1]]>=2)),
                    length(which(ap$Charge[ap$Raw.file==unique(ap$Raw.file)[10]]>=2)),
                    length(which(ap$Charge[ap$Raw.file==unique(ap$Raw.file)[4]]>=2)),
                    length(which(ap$Charge[ap$Raw.file==unique(ap$Raw.file)[5]]>=2)),
                    length(which(ap$Charge[ap$Raw.file==unique(ap$Raw.file)[3]]>=2)),
                    length(which(ap$Charge[ap$Raw.file==unique(ap$Raw.file)[6]]>=2)))
                  )
colnames(df_ap)<-c("Raw.file","feat")

# Trim out unused experiments 
df_ap2<-df_ap[c(1,3:7),]
df_ms2<-df_ms[c(1,3:7),]
df_psm2<-df_psm[c(1,3:7),-3]

# Take only peptide-like features (those with charge > +1)
ap2<-ap[ap$Charge>1, ]

# Trim out unused experiments and rename all the single cell experiments to "1 cell" 
ap2<-ap2[ap2$Raw.file%in%c("Single cell 4", "100 HeLa cells" ,"Blank","Single cell 1","Single cell 2","Single cell 3" ), ]
ap2$type<-as.character(ap2$Raw.file)
ap2$type[ap2$type%in%c("Single cell 4","Single cell 1","Single cell 2","Single cell 3" )]<-"1 cell"

# Convert data from factor to character
df_ap2$Raw.file <- as.character(df_ap2$Raw.file )
df_psm2$Raw.file<- as.character(df_psm2$Raw.file )

# Rename all the single cell experiments to "1 cell" 
df_ap2$Raw.file[3:6]<-"1 cell"
df_psm2$Raw.file[3:6]<-"1 cell"

# Combine results into one data frame
df_cb<-cbind(df_ap2, df_ms2)[,-3]
df_cb<-cbind(df_cb, df_psm2[,2])

p1<-ggplot(melt(df_cb[c(3,4,5,6,1),-5 ]), aes(x = Raw.file, fill = variable, y = value)) +
  geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize=2) +
  theme_pubr() +
  scale_x_discrete(labels=c("1 cell","100 cells")) +
  scale_fill_manual(values=c("black","grey50","lightgreen"), labels=c("Features","MS2 Scans","PSMs, 1%FDR")) +
  ylab("Count\n") +
  scale_y_log10() +
  xlab("") +
  rremove("legend") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.ticks = element_blank() )

# Calculate the percent increase observed in each statistic when going from a single cell to 100 cells
type<-c("features","features","features","features",
        "ms2","ms2","ms2","ms2",
        "psms","psms","psms","psms",
        "pm","pm","pm","pm")
benefit<-c( df_cb$feat[1]/(df_cb$feat[3:6]),
            df_cb$Scan.index[1]/(df_cb$Scan.index[3:6]),
            df_cb$`df_psm2[, 2]`[1]/(df_cb$`df_psm2[, 2]`[3:6]),
            (df_cb$`df_psm2[, 2]`[1]/ df_cb$Scan.index[1])/
              (df_cb$`df_psm2[, 2]`[3:6]/ (df_cb$Scan.index[3:6])) )

benefit<-benefit*100-100
df_x<-data.frame(type, benefit)

p2<- ggplot(df_x, aes(x=type, y=benefit, fill=type)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", position="dodge", dotsize =2) +
  scale_fill_manual(values = c("grey40","grey60","lightgreen", "springgreen3")) + 
  theme_pubr() + 
  scale_y_continuous(labels = comma, limits=c(0,1000)) + 
  scale_x_discrete(labels=c("Features","MS2 scans",  "PSMs / MS2","PSMs, 1% FDR")) +
  ylab("% Increase\n") + 
  rremove("legend") +
  xlab("") + 
  ggtitle("")+
  theme(axis.text.x = element_text(color = "grey20", size = 18, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 18, angle = 0, hjust = 1, vjust = 0, face = "plain"), 
        axis.title.x = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 18, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.ticks.x = element_blank() )

p1 + p2 + plot_layout(nrow=1, widths=c(1.1,1))
