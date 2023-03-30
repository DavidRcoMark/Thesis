#Author: David Mark
#Purpose: Outputs plots of GCSkew as calculated by skewpy (Fig. 5.5.)

library(tidyverse)
library(ggpubr)
GC <- read.table("./GCskew.tsv", header = T)
GC <- GC %>% group_by(Sequence) %>%
  mutate("Nindex"=(100*Index/max(Index))) %>% 
  filter(Sequence!="Micromonospora_L5"&Sequence!="Micromonospora_B006") %>%
  mutate("CumulIndex"=cumsum(GC_Skew_20kb))

CumGCPlot <- ggplot(GC, aes(x=Nindex, y=CumulIndex, col=Sequence, fill=Sequence))+
  geom_area()+
  theme_classic()+
  theme(legend.position = "none")+
  facet_wrap(~Sequence, ncol=1)
CumGCPlot

GCPlot <- ggplot(GC, aes(x=Nindex, y=GC_Skew_20kb, col=Sequence, fill=Sequence))+
  geom_area()+
  theme_classic()+
  theme(legend.position = "none")+
  facet_wrap(~Sequence, ncol=1)
GCPlot

Both <- ggarrange(GCPlot, CumGCPlot)

All <- ggplot(GC, aes(x=Nindex, y=GC_Skew_20kb, col=Sequence, fill=Sequence))+
  geom_area()+
  theme_classic()+
  theme(legend.position = "none")
All

AllCumul <- ggplot(GC, aes(x=Nindex, y=CumulIndex))+
  geom_area()+
  theme_classic()+
  theme(legend.position = "none")
AllCumul
  
ggsave("D:/Figures/GCSkew.png", dpi=300, height = 30, width = 10)
