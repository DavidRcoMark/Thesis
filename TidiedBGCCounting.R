#Intro----
#Author: David Mark
#Purpose: Count BGCs based on chromosomal locus (ori-proximal or ori-distal)
#Contributes to: Figure 5.19

#Import Data----
library(tidyverse)
library(ggpubr)

setwd("F:/Actinos/Micromonosporaceae/Micromonospora/BGCs/")
circBGCs <- read.csv("MicrosNormBGCs.csv",
                     header = TRUE)

#Plot Density of hybrid BGCs

HybridDensityPlot <- ggplot(data = circBGCs[grep(x = circBGCs$BGC.Type, pattern = ","),],
       aes(x=NMT))+
  geom_density(fill="white", show.legend = NA)+
  geom_point(aes(y=0, col=NMT<25|NMT>75))+
  geom_vline(xintercept = 50, col="red")+
  theme_classic()+
  theme(legend.position = "bottom")+
  ggtitle("A")+
  labs(col="Locus", x="Normalised Position", y="BGC Density")+
  scale_color_discrete(label=c("ori-proximal", "ori-distal"))

#Prep Data Frame of Hybrid BGCs ----
HybridBGCS <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = ","),]
HybridBGCS
HybridRowN <- row.names(HybridBGCS)

OriHybrids <- 
  HybridBGCS %>% filter(NMT>25 & NMT<75)
OHClustercount <- count(OriHybrids,Host.Organism)
OHClustercount <- data.frame("Organism"=OHClustercount$Host.Organism, 
                             "OriginHybridClusters"=OHClustercount$n)

TerHybrids <-
  HybridBGCS %>% filter(NMT<25|NMT>75)
THClustercount <- count(TerHybrids,Host.Organism)
THClustercount <- data.frame("Organism" = THClustercount$Host.Organism,
                             "TerHybridClusters" =THClustercount$n)

#Merge Cluster counts and count organisms without clusters ----
HybridClusters <- full_join(THClustercount, 
                            OHClustercount,
                            by="Organism")
HybridClusters[is.na(HybridClusters)] <- 0

wilcox.test(HybridClusters$TerHybridClusters,
            HybridClusters$OriginHybridClusters,
            paired = T,
            alternative = "greater")

#Prep data frame for plotting ----
HybridPaird <- data.frame("Organism"=c(HybridClusters$Organism, 
                                       HybridClusters$Organism),
                          "Values"=c(HybridClusters$OriginHybridClusters,
                                     HybridClusters$TerHybridClusters),
                          "Locus"=c(rep("1Origin", 
                                        each=28), 
                                    rep("2Mid-Chromosome", 
                                        each=28)))

#Plot the frame ----
HybridClusterPlot <- ggplot(HybridPaird)+
  geom_violin(aes(x=Locus, y=Values, fill=factor(Locus)))+
  geom_point(aes(x=Locus, y=Values))+
  geom_line(aes(x=c(Locus),
                y=Values,
                group=Organism))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("Number of Clusters")+
  scale_y_continuous(breaks = 0:15)+
  scale_x_discrete(labels = c("ori-proximal", "ori-distal"))+
  geom_label(aes(x="1Origin", y=6, label="P=2.756e-06"))+
  ggtitle("B")
HybridClusterPlot


#Repeat for Non-Hybrids----
#Prep Data Frame of Non-Hybrid BGCs----
NonhybridBGCs <- circBGCs[-as.numeric(HybridRowN), ] 
NonhybridBGCs

OriNH <- 
  NonhybridBGCs %>% filter(NMT>25 & NMT<75)
ONHClustercount <- count(OriNH,Host.Organism)
ONHClustercount <- data.frame("Organism"=ONHClustercount$Host.Organism, 
                              "OriginNonHybridClusters"=ONHClustercount$n)

TerNH <-
  NonhybridBGCs %>% filter(NMT<25|NMT>75)
TerNHClustercount <- count(TerNH, Host.Organism)
TerNHClustercount <- data.frame("Organism"=TerNHClustercount$Host.Organism, 
                                "TerNonHybridClusters"=TerNHClustercount$n)

#Merge Cluster counts and count organisms without clusters----
NonHybridClusters <- full_join(TerNHClustercount, 
                            ONHClustercount,
                            by="Organism")
NonHybridClusters[is.na(NonHybridClusters)] <- 0

wilcox.test(NonHybridClusters$TerNonHybridClusters,
            NonHybridClusters$OriginNonHybridClusters,
            paired = T)

#Prep data frame for plotting----
NonHybridPaird <- data.frame("Organism"=c(NonHybridClusters$Organism, 
                                          NonHybridClusters$Organism),
                            "Values"=c(NonHybridClusters$OriginNonHybridClusters,
                                      NonHybridClusters$TerNonHybridClusters),
                            "Locus"=c(rep("1Origin", 
                                        each=28), 
                                    rep("2Mid-Chromosome", 
                                        each=28)))
NonHybridPaird

NonHybridClusterPlot <- ggplot(NonHybridPaird)+
  geom_violin(aes(x=Locus, y=Values, fill=factor(Locus)))+
  geom_point(aes(x=Locus, y=Values))+
  geom_line(aes(x=c(Locus),
                y=Values,
                group=Organism))+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("Number of Clusters")+
  ggtitle("C")+
  scale_y_continuous(breaks = 0:15)+
  scale_x_discrete(labels = c("ori-proximal", "ori-distal"))+
  geom_label(aes(x="1Origin", y=14, label="P=0.00424"))
NonHybridClusterPlot

#Draw Figure ----
ArrangedFig <- ggarrange(HybridDensityPlot, HybridClusterPlot, NonHybridClusterPlot, ncol = 1)
ggsave(ArrangedFig, filename =  "Filename", dpi=300, width = 5, height = 15)
