#Author: David Mark
#Plots chord diagram of BGCs present in organisms (Fig. 4.8.)

setwd("D:/Actinos/Micromonosporaceae/Micromonospora/BGCs/AtacamaStrains/")
Atacamas <- read.csv("BGCsforcircos.csv", header = F)
Atacamas <- Atacamas[,1:2]
colnames(Atacamas)=c("Organism", "Class")

library(tidyverse)
require(circlize)
view(Atacamas)

ggplot(Atacamas)+
  geom_bar(aes(x=Organism, fill=Class), colour="black")+
  theme_classic()

AtacamasForChordDiagram <- Atacamas %>% select(Class, Organism)

Chordcolours <- c("M. sp. O5"="#ff0010",
                  "M. sp. O3"="#0075dc",
                  "M. sp. PH63"="#2bce48",
                  "betalactone"="#FFFFFF",
                  "lanthipeptide-class-i"="#f0a3ff",
                  "lanthipeptide-class-ii"="#993f00",
                  "NAGGN"="#4c005f",
                  "NRPS"="#191919",
                  "NRPS,betalactone,transAT-PKS,lanthipeptide-class-ii"="#005c31",
                  "NRPS,lanthipeptide-class-i"="#ffcc99",
                  "NRPS,lanthipeptide-class-ii"="#808080",
                  "NRPS,LAP,T1PKS"="#94ffb5",
                  "NRPS,other,betalactone"="#8f7c00",
                  "NRPS,siderophore,T1PKS"="#9dcc00",
                  "NRPS.T1PKS"="#c20088",
                  "NRPS,T1PKS,LAP"="#003380",
                  "NRPS,thioamitides,lanthipeptide-class-ii"="#ffa405",
                  "other"="#ffa8bb",
                  "other,NRPS"="#426600",
                  "RiPP-like"="#5ef1f2",
                  "siderophore"="#00998f",
                  "T1PKS"="#e0ff66",
                  "T1PKS,NRPS"="#740aff",
                  "T1PKS,NRPS,betalactone,lanthipeptide-class-ii"="#ffff80",
                  "T1PKS,NRPS,siderophore"="#ffff00",
                  "T2PKS"="#ff5005",
                  "T3PKS"="Cyan",
                  "terpene"="darkorchid4",
                  "thiopeptide,LAP"="Magenta")
Chordcolours

BGCchords <- chordDiagram(AtacamasForChordDiagram,
             grid.col = Chordcolours,
             grid.border = "Black",
             link.border = "Black",
             transparency = 0.1, 
            )

svg("./ChordDiagramTest")
