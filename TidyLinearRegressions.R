#Author: David Mark
#Purpose: Linear regression of chromosome size v. no. of BGCs and chromosome size v. percentage of chromosome occupied by BGCs (Fig. 5.8)
#Import data----
library(tidyverse)
setwd("D:/Actinos/Micromonosporaceae/Micromonospora/BGCs/")
circBGCs <- read.csv("MicrosNormBGCs.csv", header = TRUE)
circBGCs

#Linear regresion of No. BGCs vs. chromosome size ----
LinearReg <-lm(data = (circBGCs %>% count(Chromosome.Size)),
               formula= n~Chromosome.Size)
summary(LinearReg)

LinearRegplot <- circBGCs %>% 
  count(Chromosome.Size) %>% 
  ggplot(aes(x=Chromosome.Size, y=n))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()+
  ylab("Number of Clusters")+
  xlab("Chromosome Size (bp)")+
  ggtitle("A")+
  geom_label(label="Adjusted R^2=0.32, P=0.0009844",
             x=median(circBGCs$Chromosome.Size),
             y=35)
LinearRegplot

#Lijnear regression of Chromosome commitment vs BGC size ----
BGCCommitment <- circBGCs %>% 
  group_by(Chromosome.Size) %>% 
  summarise(sum(Size)) %>%
  mutate(Commitment = `sum(Size)`/Chromosome.Size*100)

CommitmentPlot <-
  ggplot(BGCCommitment,
         aes(x = Chromosome.Size, y = Commitment)) +
  geom_point() +
  geom_smooth(method = "lm")+
  theme_classic()+
  xlab("Chromosome Size (bp)")+
  ylab("Chromosomal Commitment (%)")+
  ggtitle("B")+
  geom_label(label="Adjusted R^2=0.20, P=0.00949",
             x=median(circBGCs$Chromosome.Size),
             y=21)
CommitmentPlot

summary(lm(data = BGCCommitment, formula=Commitment~Chromosome.Size))

ggsave(LinearRegplot, filename = "F:/Figures/CountvSize.png", dpi=300)
ggsave(CommitmentPlot, filename = "F:/Figures/CommitmentPlot.png", dpi=300)