
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

xx1<-read.delim("dat/evidence_whs271.txt")

kp<-intersect(colnames(ev1), colnames(xx1))

ev1<-rbind(xx1[xx1$Raw.file%in%"whs271", kp], ev1[,kp])


# Correct autosampler mixup of 200 and 300-cell carrier samples
ev1$rf<-ev1$Raw.file
ev1$rf[ev1$Raw.file=="hs0042"]<-"hs0043"
ev1$rf[ev1$Raw.file=="hs0043"]<-"hs0042"
ev1$Raw.file<-as.factor(ev1$rf)

# Filter peptides to <1% FDR and remove reverse and contaminant hits
ev1<-ev1[ev1$PEP<0.1, ]
#ev1<-ev1[-which(ev1$Potential.contaminant=="+"), ]
#ev1<-ev1[-which(ev1$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev1$ms<-paste0( ev1$Modified.sequence, ev1$Charge )

# Remove duplicate peptides
ev1<-remove.duplicates(ev1, c("ms", "Raw.file"))

# Code for plot adopted from DO-MS publication, Huffman et al. 
plotdata <- ev1[,c('Raw.file', 'PEP','Type')]
plotdata <- plotdata %>% dplyr::filter(Type != "MULTI-MATCH")
plotdata <- plotdata %>% dplyr::select('Raw.file', 'PEP')

peps <- seq(log10(max(c(min(plotdata$PEP)), 1e-5)), log10(max(plotdata$PEP)), length.out=500)
peps <- c(log10(.Machine$double.xmin), peps)

plotdata <- plotdata %>%
  dplyr::mutate(bin=findInterval(PEP, 10**peps)) %>%
  dplyr::group_by(Raw.file, bin) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::mutate(cy=cumsum(n),
                pep=10**peps[bin])


maxnum <- c()
rawnames <- c()

for(X in unique(plotdata$Raw.file)){
  maxnum <- c(maxnum, max(plotdata$cy[plotdata$Raw.file %in% X]) )
  rawnames <- c(rawnames, X)
}

names(maxnum) <- rawnames
rank_exp <- maxnum[order(maxnum, decreasing = T)]
rank_exp_ord <- seq(1, length(rank_exp),1)
names(rank_exp_ord) <- names(rank_exp)
plotdata$rank_ord <- NA

for(X in levels(plotdata$Raw.file)) {
  plotdata$rank_ord[plotdata$Raw.file %in% X] <- rank_exp_ord[X]
}

cc <- scales::seq_gradient_pal('darkgreen', 'orange', 'Lab')(seq(0, 1, length.out=(length(rank_exp_ord)-1)))

cc<-c(cc, cc[length(cc)])

plotdata$Resolution<-"70k"
plotdata$Resolution[plotdata$Raw.file%in%c("whs271")]<-"35k"

p1<-ggplot(plotdata, aes(x=pep, color=Raw.file, y=cy, group=Raw.file, lty=Resolution)) + 
  geom_line(size = 1.4) +
  scale_colour_manual(name='Experiment', values=cc, labels=names(rank_exp_ord)) +
  #scale_colour_manual(name='Experiment', values=c("#F88000","#008000","#F88000","#008000")) +
  #scale_colour_manual(name='Experiment', values=c(rgb(0.3,0.8,0,3/4),rgb(0.3,0.9,0,2/4),rgb(0,0,0.8,3/4),rgb(0,0,0.9,2/4))) +
  coord_flip() + 
  theme_bw()+
  ylim(c(0,6000))+
  scale_x_log10(limits=c(.00009,.1), breaks=c(.0001,.001,.01,.1), 
                labels=scales::trans_format('log10', scales::math_format(10^.x))) + 
  theme(legend.key=element_rect(fill='white')) +
  xlab('Posterior error probability\n') + ylab('\nPeptides ordered by PEP') + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))


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
ev1$rf[ev1$Raw.file=="hs0101"]<-"hs0102"
ev1$rf[ev1$Raw.file=="hs0102"]<-"hs0101"
ev1$Raw.file<-as.factor(ev1$rf)

# Filter peptides to <1% FDR and remove reverse and contaminant hits
ev1<-ev1[ev1$PEP<0.1, ]
#ev1<-ev1[-which(ev1$Potential.contaminant=="+"), ]
#ev1<-ev1[-which(ev1$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev1$ms<-paste0( ev1$Modified.sequence, ev1$Charge )

# Remove duplicate peptides
ev1<-remove.duplicates(ev1, c("ms", "Raw.file"))



plotdata <- ev1[,c('Raw.file', 'PEP','Type')]
plotdata <- plotdata %>% dplyr::filter(Type != "MULTI-MATCH")
plotdata <- plotdata %>% dplyr::select('Raw.file', 'PEP')

peps <- seq(log10(max(c(min(plotdata$PEP)), 1e-5)), log10(max(plotdata$PEP)), length.out=500)
peps <- c(log10(.Machine$double.xmin), peps)

plotdata <- plotdata %>%
  dplyr::mutate(bin=findInterval(PEP, 10**peps)) %>%
  dplyr::group_by(Raw.file, bin) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::mutate(cy=cumsum(n),
                pep=10**peps[bin])


maxnum <- c()
rawnames <- c()

for(X in unique(plotdata$Raw.file)){
  maxnum <- c(maxnum, max(plotdata$cy[plotdata$Raw.file %in% X]) )
  rawnames <- c(rawnames, X)
}

names(maxnum) <- rawnames
rank_exp <- maxnum[order(maxnum, decreasing = T)]
rank_exp_ord <- seq(1, length(rank_exp),1)
names(rank_exp_ord) <- names(rank_exp)
plotdata$rank_ord <- NA

for(X in levels(plotdata$Raw.file)) {
  plotdata$rank_ord[plotdata$Raw.file %in% X] <- rank_exp_ord[X]
}

cc <- scales::seq_gradient_pal('darkgreen', 'orange', 'Lab')(seq(0, 1, length.out=(length(rank_exp_ord))))

plotdata$Resolution<-"70k"

p2<-ggplot(plotdata, aes(x=pep, color=Raw.file, y=cy, group=Raw.file, lty=Resolution)) + 
  geom_line(size = 1.4) +
  scale_colour_manual(name='Experiment', values=cc, labels=names(rank_exp_ord)) +
  #scale_colour_manual(name='Experiment', values=c("#F88000","#008000","#F88000","#008000")) +
  #scale_colour_manual(name='Experiment', values=c(rgb(0.3,0.8,0,3/4),rgb(0.3,0.9,0,2/4),rgb(0,0,0.8,3/4),rgb(0,0,0.9,2/4))) +
  coord_flip() + 
  theme_bw()+
  ylim(c(0,6000))+
  scale_x_log10(limits=c(.00009,.1), breaks=c(.0001,.001,.01,.1), 
                labels=scales::trans_format('log10', scales::math_format(10^.x))) + 
  theme(legend.key=element_rect(fill='white')) +
  xlab('Posterior error probability\n') + ylab('\nPeptides ordered by PEP') + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        title = element_text(color = "grey20", size = 15, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))




# Visualize
p1 + p2 +plot_layout(ncol=2, nrow=1)


