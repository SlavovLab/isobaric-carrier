##### High AGC series
####
##
#


# Load files: 
exps<-c(paste0("dat/carrier_202004xx/hs0",100:105,"/msms.txt") )
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
ev1<-ev1[ev1$PEP<0.02, ]
ev1<-ev1[-which(ev1$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev1$ms<-paste0( ev1$Modified.sequence, ev1$Charge )

# Remove duplicate peptides
ev1<-remove.duplicates(ev1, c("ms", "Raw.file"))

# Take only peptides observed in all 6 experiments
table_peps<-table(ev1$ms)
intersect_peps1<-names(table_peps)[table_peps==length(unique(ev1$Raw.file))]


########################################################################
##### Low AGC series
####
##
#

# Load files: 
exps<-c(paste0("dat/carrier_202004xx/hs00",41:46,"/msms.txt") )
ev2<-read.delim(exps[1])

for(X in exps[-c(1)]){
  
  ev2<-rbind(ev2, read.delim(X))
  
}

# Correct autosampler mixup of 200 and 300-cell carrier samples
ev2$rf<-ev2$Raw.file
ev2$rf[ev2$Raw.file=="hs0042"]<-"hs0043"
ev2$rf[ev2$Raw.file=="hs0043"]<-"hs0042"
ev2$Raw.file<-as.factor(ev2$rf)

# Filter peptides to <1% FDR and remove reverse and contaminant hits
ev2<-ev2[ev2$PEP<0.02, ]
ev2<-ev2[-which(ev2$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev2$ms<-paste0( ev2$Modified.sequence, ev2$Charge )

# Remove duplicate peptides
ev2<-remove.duplicates(ev2, c("ms", "Raw.file"))

# Take only peptides observed in all 6 experiments
table_peps<-table(ev2$ms)
intersect_peps2<-names(table_peps)[table_peps==length(unique(ev2$Raw.file))]

#############################

ev1$`MS2 AGC`<-"1e6"
ev2$`MS2 AGC`<-"5e4"

ev1$obs_every_exp<-"no"
ev1$obs_every_exp[ev1$ms%in%intersect_peps1]<-"yes"

ev2$obs_every_exp<-"no"
ev2$obs_every_exp[ev2$ms%in%intersect_peps1]<-"yes"

ev1$cells_in_carrier <- c(100,200,300,400,600,800) [as.numeric( factor( ev1$Raw.file ) ) ]
ev2$cells_in_carrier <- c(100,200,300,400,600,800) [as.numeric( factor( ev2$Raw.file ) ) ]

ev3<-rbind(ev1,ev2)

carrier_labs<-c( c(100,200,300,400,600,800), c(100,200,300,400,600,800) )
#############################

# Score 

ev3$Score[ev3$Score>200]<-200
ggplot(ev3[ev3$obs_every_exp=="yes",], aes(y=rf, x=Score, fill=`MS2 AGC`)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_manual(values=c(rgb(0.8,0,0),rgb(0,0,0.8))) + 
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
  ylab("\n Number of cells in isobaric carrier") + 
  ggtitle("") +
  xlab("Andromeda Score\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

#############################

# Fragments 

ev3$Number.of.matches[ev3$Number.of.matches>40]<-40
ggplot(ev3[ev3$obs_every_exp=="yes",], aes(y=rf, x=Number.of.matches, fill=`MS2 AGC`)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_manual(values=c(rgb(0.8,0,0),rgb(0,0,0.8))) + 
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
  ylab("\n Number of cells in isobaric carrier") + 
  ggtitle("") +
  xlab("# Fragments\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

#############################

# small sample RIs - every experiment 

ev3$scRI<-log10( rowMeans(ev3[,76:83], na.rm=T) )
ggplot(ev3[ev3$obs_every_exp=="yes",], aes(y=rf, x=scRI, fill=`MS2 AGC`)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_manual(values=c(rgb(0.8,0,0),rgb(0,0,0.8))) + 
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
  ylab("\n Number of cells in isobaric carrier") + 
  ggtitle("Peptides observed in every experiment") +
  xlab("Small sample\n RI intensity, log10\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))

#############################

# small sample RIs - not every experiment 

ggplot(ev3[ev3$obs_every_exp=="no",], aes(y=rf, x=scRI, fill=`MS2 AGC`)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_manual(values=c(rgb(0.8,0,0),rgb(0,0,0.8))) + 
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
  ylab("\n Number of cells in isobaric carrier") + 
  ggtitle("Peptides not observed in every experiment") +
  xlab("Small sample\n RI intensity, log10\n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))


##### High AGC series
####
##
#


# Load files: 
exps<-c(paste0("dat/carrier_202004xx/hs0",100:105,"/msmsScans.txt") )
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
# ev1<-ev1[ev1$PEP<0.02, ]
# ev1<-ev1[-which(ev1$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev1$ms<-paste0( ev1$Modified.sequence, ev1$Charge )

# Remove duplicate peptides
ev1<-remove.duplicates(ev1, c("ms", "Raw.file"))

# Take only peptides observed in all 6 experiments
table_peps<-table(ev1$ms)
intersect_peps1<-names(table_peps)[table_peps==length(unique(ev1$Raw.file))]


########################################################################
##### Low AGC series
####
##
#

# Load files: 
exps<-c(paste0("dat/carrier_202004xx/hs00",41:46,"/msmsScans.txt") )
ev2<-read.delim(exps[1])

for(X in exps[-c(1)]){
  
  ev2<-rbind(ev2, read.delim(X))
  
}

# Correct autosampler mixup of 200 and 300-cell carrier samples
ev2$rf<-ev2$Raw.file
ev2$rf[ev2$Raw.file=="hs0042"]<-"hs0043"
ev2$rf[ev2$Raw.file=="hs0043"]<-"hs0042"
ev2$Raw.file<-as.factor(ev2$rf)

# Filter peptides to <1% FDR and remove reverse and contaminant hits
# ev2<-ev2[ev2$PEP<0.02, ]
# ev2<-ev2[-which(ev2$Reverse=="+"), ]

# Create a variable for ions (Sequence+charge)
ev2$ms<-paste0( ev2$Modified.sequence, ev2$Charge )

# Remove duplicate peptides
ev2<-remove.duplicates(ev2, c("ms", "Raw.file"))

# Take only peptides observed in all 6 experiments
table_peps<-table(ev2$ms)
intersect_peps2<-names(table_peps)[table_peps==length(unique(ev2$Raw.file))]

#############################

ev1$`MS2 AGC`<-"1e6"
ev2$`MS2 AGC`<-"5e4"

ev1$obs_every_exp<-"no"
ev1$obs_every_exp[ev1$ms%in%intersect_peps1]<-"yes"

ev2$obs_every_exp<-"no"
ev2$obs_every_exp[ev2$ms%in%intersect_peps1]<-"yes"

ev1$cells_in_carrier <- c(100,200,300,400,600,800) [as.numeric( factor( ev1$Raw.file ) ) ]
ev2$cells_in_carrier <- c(100,200,300,400,600,800) [as.numeric( factor( ev2$Raw.file ) ) ]

ev3<-rbind(ev1,ev2)

carrier_labs<-c( c(100,200,300,400,600,800), c(100,200,300,400,600,800) )

#############################
# Fill time 

ggplot(ev3[ev3$obs_every_exp=="yes",], aes(y=rf, x=Ion.injection.time, fill=`MS2 AGC`)) +
  geom_density_ridges(stat = "binline", bins = 30, scale = 0.95, draw_baseline = FALSE) +
  scale_fill_manual(values=c(rgb(0.8,0,0),rgb(0,0,0.8))) + 
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
  ylab("\n Number of cells in isobaric carrier") + 
  ggtitle("") +
  xlab("Fill time (ms) \n") + 
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 30, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))


