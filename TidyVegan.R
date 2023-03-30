#Author: David Mark
#Purpose: Calculate diversity of BGCs in Micromonospora ori-proximal and ori-distal regions (Fig. 5. )
#Ecological comparisons of Ori and Ter regions in Micromonospora chromosomes

#Import packages----
library(tidyverse)
library(ggpubr)
library(vegan)
setwd("D:/Actinos/Micromonosporaceae/Micromonospora/Diversity/")

#Load presence-absence tables----
oriPAT <- read.csv(
  "./BGCpresence_absence_tables/oriclustersfreqtable.csv", 
  header = TRUE,
  row.names = 1
  )

terPAT <- read.csv(
  "./BGCpresence_absence_tables/terbgctable.csv", 
  header = TRUE,
  row.names = 1
  )

#Convert to matrix----
MoriPAT <- as.matrix(oriPAT)
MterPAT <- as.matrix(terPAT)

#Calculate richness and diversity of loci----
oriDsh <- diversity(MoriPAT, index = "shannon")
oriDsi <- diversity(MoriPAT, index = "simpson")
oriDinv <- diversity(MoriPAT, index = "invsimpson")
oriRichness <- specnumber(MoriPAT)

terDsh <- diversity(MterPAT, index = "shannon")
terDsi <- diversity(MterPAT, index = "simpson")
terDinv <- diversity(MterPAT, index = "invsimpson")
terRichness <- specnumber(MterPAT)

#Convert diversity statistics into data frames to be plotted----
oriDshdf <- as.data.frame(oriDsh)
oriDsidf <- as.data.frame(oriDsi)
oriDinvdf <- as.data.frame(oriDinv)
terDshdf <- as.data.frame(terDsh)
terDsidf <- as.data.frame(terDsi)
terDinvdf <- as.data.frame(terDinv)

#Merge ori and ter dfs----
ShanDF <- bind_cols(oriDshdf, terDshdf)
SimpsonDF <- bind_cols(oriDsidf, terDsidf)
InvDF <- bind_cols(oriDinvdf, terDinvdf)

#Plot density estimates of diversity----
##Plot shannon diversity----
Shanplot <- ggplot()+
  geom_density(data=ShanDF,
               aes(x=terDsh,
                   fill=factor("Terminus"),
                   alpha=0.5))+
  geom_density(data = ShanDF,
               aes(x = oriDsh,
                   fill = factor("Origin"),
                   alpha=0.5))+
  xlab("Shannon Diversity")+
  guides(alpha="none")+
  theme(legend.title = element_blank())
Shanplot

##Plot Simpson diversity----
Simpsonplot <- ggplot()+
  geom_density(data=SimpsonDF,
               aes(x=terDsi,
                   fill=factor("Terminus"),
                   alpha=0.5))+
  geom_density(data = SimpsonDF,
               aes(x = oriDsi,
                   fill = factor("Origin"),
                   alpha=0.5))+
  xlab("Simpson Diversity")+
  guides(alpha="none")+
  theme(legend.title = element_blank())
Simpsonplot

##Plot Inverse Simpson diversity----
Invplot <- ggplot()+
  geom_density(data=InvDF,
               aes(x=terDinv,
                   fill=factor("Terminus"),
                   alpha=0.5))+
  geom_density(data = InvDF,
               aes(x = oriDinv,
                   fill = factor("Origin"),
                   alpha=0.5))+
  xlab("Inverse Simpson Diversity")+
  guides(alpha="none")+
  theme(legend.title = element_blank())
Invplot

#Plot Violin Plots of Diversities----
##Make appropriate data frames----
###Shannon Diversity Frame----
ShanlineDF <- data.frame("Values"=c(ShanDF$oriDsh, ShanDF$terDsh), "Host"=c(row.names(ShanDF), row.names(ShanDF)), "Locus"=c("oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh","oriDsh", "terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh","terDsh"))

###Simpson Diversity Frame----
SimpLineDF <- data.frame("Values"=c(SimpsonDF$oriDsi, SimpsonDF$terDsi), "Host"=ShanlineDF$Host, "Locus"=ShanlineDF$Locus)

###Inverse Simpson Diversity Frame----
InvLineDF <- data.frame("Values"=c(InvDF$oriDinv, InvDF$terDinv), "Host"=ShanlineDF$Host, "Locus"=ShanlineDF$Locus)

##Wilcox Tests ----

ShanTest <- wilcox.test(ShanDF$oriDsh, ShanDF$terDsh, 
                        paired = TRUE,
                        alternative = "less",
                        exact = F)
InvSimpTest <- wilcox.test(InvDF$oriDinv, InvDF$terDinv,
                           paired = TRUE,
                           alternative = "less",
                           exact = F)
SimpTest <- wilcox.test(SimpsonDF$oriDsi, SimpsonDF$terDsi, 
                        paired = TRUE, 
                        alternative = "less",
                        exact = F)

##ggPlot the Violins----
Shanbox <- 
  ggplot(ShanlineDF)+
  geom_violin(aes(x=factor(Locus), y=Values, fill=factor(Locus)))+
  geom_point(aes(x=factor(Locus), y=Values))+
  geom_line(aes(y=Values, x=factor(Locus), group=Host))+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("ori-proximal", "ori-distal"))+
  xlab("Locus")+
  ylab("Shannon Diversity")+
  ggtitle("A")+
  geom_text(aes(y=2.5, x=1), label="p=2.48e-06")
Shanbox


Invbox <-
  ggplot(InvLineDF)+
  geom_violin(aes(x=factor(Locus), y=Values, fill=factor(Locus)))+
  geom_point(aes(x=factor(Locus), y=Values))+
  geom_line(aes(y=Values, x=factor(Locus), group=Host))+
  scale_x_discrete(labels=c("ori-proximal", "ori-distal"))+
  xlab("Locus")+
  ylab("Inverse Simpson Diversity")+
  theme_classic()+
  theme(legend.position = "none")+
  geom_text(aes(label="p=2.23e-06", x="oriDsh", y=10.5))+
  ggtitle("B")
Invbox

ggsave(plot=Invbox, filename = "Invbox.png", dpi=300)

Simpbox <-
  ggplot(SimpLineDF)+
  geom_violin(aes(x=factor(Locus), y=Values, fill=factor(Locus)))+
  geom_point(aes(x=factor(Locus), y=Values))+
  geom_line(aes(y=Values, x=factor(Locus), group=Host))+
  theme_classic()+
  theme(legend.position = "none")+
  scale_x_discrete(labels=c("ori-proximal", "ori-distal"))+
  xlab("Locus")+
  ylab("Simpson Diversity")+
  geom_text(aes(label="p=2.23e-06", x="oriDsh", y=0.93))+
  ggtitle("C")
Simpbox

A1 <- ggarrange(Shanbox, Invbox, Simpbox, ncol = 1)
ggsave(file="D:/Figures/Divplots.png", dpi=300, height = 15, width = 5)
