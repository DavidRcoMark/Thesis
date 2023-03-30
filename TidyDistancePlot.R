#Author: David Mark
#Purpose: Plot BGCs by distance from origin of replication (Fig. 5.10)

library(tidyverse)
library(ggpubr)
library(ggnewscale)
setwd("F:/Actinos/")
circBGCs <- read.csv("D:/Actinos/Micromonosporaceae/Micromonospora/BGCs/MicrosNormBGCs.csv", header = TRUE)

#Generate DF of biosynthetic features of interest ----
Nurps <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "NRPS"),]
NAGGN <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "NAGGN"),]
T3PKS <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "T3PKS"),]
T2PKS <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "T2PKS"),]
Turps <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "terpene"),]
Sides <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "siderophore"),]
T1PKS <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "T1PKS"),]

#Plot distance from origin ----
DistancePlot <- ggplot(data=Nurps, aes(x="NRPS",
                                       y=sqrt((NMT-50)^2)))+
  geom_violin()+
  geom_violin(data = Turps, aes(x="Terpene"))+
  geom_violin(data = Sides, aes(x="Siderophore"))+
  geom_violin(data = T2PKS, aes(x="T2PKS"))+
  geom_violin(data = T3PKS, aes(x="T3PKS"))+
  geom_violin(data = T1PKS, aes(x="T1PKS"))+
  geom_violin(data = NAGGN, aes(x="NAGGN"))+
  new_scale_fill()+
  geom_jitter(aes(fill=sqrt((Norm.Mid-50)^2)<25), shape=21, col="black", size=2)+
  geom_jitter(data = Turps, aes(x="Terpene", fill=sqrt((Norm.Mid-50)^2)<25), shape=21, col="black", size=2)+
  geom_jitter(data = Sides, aes(x="Siderophore", fill=sqrt((Norm.Mid-50)^2)<25), shape=21, col="black", size=2)+
  geom_jitter(data = T2PKS, aes(x="T2PKS", fill=sqrt((Norm.Mid-50)^2)<25), shape=21, col="black", size=2)+
  geom_jitter(data = T3PKS, aes(x="T3PKS", fill=sqrt((Norm.Mid-50)^2)<25), shape=21, col="black", size=2)+
  geom_jitter(data = T1PKS, aes(x="T1PKS", fill=sqrt((Norm.Mid-50)^2)<25), shape=21, col="black", size=2)+
  geom_jitter(data = NAGGN, aes(x="NAGGN", fill=sqrt((Norm.Mid-50)^2)<25), shape=21, col="black", size=2)+
  theme_classic()+
  labs(x = "Cluster Type", y = "Distance from Origin", fill="ori-distal")+
  theme(legend.position = "bottom")

ggsave("D:/Figures/BGCDistancePlot.png", DistancePlot, width = 5, height = 5, dpi=300)

DistancePlot
