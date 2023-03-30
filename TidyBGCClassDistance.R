#Intro ----
#Author: David Mark
#Purpose: Plot distance of BGCs from the origin of replication
#Contributes to: Fig 5.10

#Import data ----
library(tidyverse)
setwd("F:/Actinos/Micromonosporaceae/Micromonospora/BGCs/")
circBGCs <- read.csv("MicrosNormBGCs.csv", header = TRUE)

#Organise BGCs by function ----

Nurps <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "NRPS"),]
Naggn <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "NAGGN"),]
T3PKS <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "T3PKS"),]
T2PKS <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "T2PKS"),]
Turps <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "terpene"),]
Sides <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "siderophore"),]
T1PKS <- circBGCs[grep(x = circBGCs$BGC.Type, pattern = "T1PKS"),]

#Plot BGC Distances ----
DistancePlot <- ggplot(data=Nurps, aes(x="NRPS", y=sqrt((NMT-50)^2)))+
  geom_violin(fill = "darkgreen")+
  geom_violin(data = Naggn, aes(x="NAGGN"), fill = "green")+
  geom_violin(data = Turps, aes(x="Terpene"), fill = "purple")+
  geom_violin(data = Sides, aes(x="Siderophore"), fill = "darkred")+
  geom_violin(data = T2PKS, aes(x="T2PKS"), fill = "orange")+
  geom_violin(data = T3PKS, aes(x="T3PKS"), fill = "darkorange")+
  geom_violin(data = T1PKS, aes(x="T1PKS"), fill = "pink")+
  geom_jitter()+
  geom_jitter(data = Naggn, aes(x="NAGGN"))+
  geom_jitter(data = Turps, aes(x="Terpene"))+
  geom_jitter(data = Sides, aes(x="Siderophore"))+
  geom_jitter(data = T2PKS, aes(x="T2PKS"))+
  geom_jitter(data = T3PKS, aes(x="T3PKS"))+
  geom_jitter(data = T1PKS, aes(x="T1PKS"))+
  theme_classic()+
  labs(x = "Cluster Type", y = "Distance from Origin")+
  geom_abline(intercept = 25, 
              slope = 0,
              colour = "red")
DistancePlot

#ggsave(plot = DistancePlot, filename = "F:/Figures/Distanceplot.png")

#Print positions of BGCs ----
list("NRPS"=summary(sqrt((Nurps$NMT-50)^2)),
"NAGGN"=summary(sqrt((Naggn$NMT-50)^2)),
"Terpene"=summary(sqrt((Turps$NMT-50)^2)),
"Siderophore"=summary(sqrt((Sides$NMT-50)^2)),
"T2PKS"=summary(sqrt((T2PKS$NMT-50)^2)),
"T3PKS"=summary(sqrt((T3PKS$NMT-50)^2)),
"T1PKS"=summary(sqrt((T1PKS$NMT-50)^2)))
