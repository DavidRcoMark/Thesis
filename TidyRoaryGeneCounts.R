library(tidyverse)
library(ggpubr)
#Tidy Unique Genes ----
setwd("D:/Treework/20221025AutoMLSTRun/Roary/Roarywopataloonsplitparalogs/")
UniqueGenes <- read.table("./Roary on data 181, data 179, and others Number of Unique Genes.txt", sep="\t")
UniqueGenes <- data.frame(t(UniqueGenes)) %>%
  mutate("Number of Chromosomes"=c(1:ncol(UniqueGenes))) %>%
  pivot_longer(cols = c(1:10))

UniqueGenes <- UniqueGenes[,-2]
colnames(UniqueGenes)=c("Number of Chromosomes", "Unique Genes")

UniqueGenes %>% ggplot(aes(x = `Number of Chromosomes`, y=`Unique Genes`))+
                         geom_smooth()
#Tidy Conserved Genes ----
ConservedGenes <- read.table("./Roary on data 181, data 179, and others Number of Conserved Genes.txt", sep="\t")
ConservedGenes <- data.frame(t(ConservedGenes)) %>%
  mutate("Number of Chromosomes"=c(1:ncol(ConservedGenes))) %>%
  pivot_longer(cols = c(1:10))

ConservedGenes <- ConservedGenes[,-2]
colnames(ConservedGenes)=c("Number of Chromosomes", "Conserved Genes")
ConservedGenes %>% ggplot(aes(x=`Number of Chromosomes`, y=`Conserved Genes`))+
  geom_smooth()

#Tidy Total Genes ----
TotalGenes <- read.table("./Roary on data 181, data 179, and others Number of Genes in Pan Geneome.txt", sep="\t")
TotalGenes <- data.frame(t(TotalGenes)) %>%
  mutate("Number of Chromosomes"=c(1:ncol(TotalGenes))) %>%
  pivot_longer(cols = c(1:10))

TotalGenes <- TotalGenes[,-2]
colnames(TotalGenes)=c("Number of Chromosomes", "Total Genes")

#Tidy New Genes ----
NewGenes <- read.table("./Roary on data 181, data 179, and others Number of New Genes.txt", sep="\t")

NewGenes <- data.frame(t(NewGenes)) %>%
  mutate("Number of Chromosomes"=c(1:ncol(NewGenes))) %>%
  pivot_longer(cols = c(1:10))

NewGenes <- NewGenes[,-2]
colnames(NewGenes)=c("Number of Chromosomes", "New Genes")

ggplot(NewGenes, aes(x=`Number of Chromosomes`, y=`New Genes`))+
  geom_smooth()

##Plot ----

RoaryTotalAndUnique <- TotalGenes %>% ggplot(aes(x=`Number of Chromosomes`, y=`Total Genes`, fill=factor("Total Genes"), col=factor("Total Genes")))+
  geom_smooth()+
  geom_point(alpha=0.5)+
  geom_smooth(data = UniqueGenes, aes(x=`Number of Chromosomes`, y=`Unique Genes`, fill=factor("Unique Genes"), col=factor("Unique Genes")))+
  geom_point(data = UniqueGenes, alpha=0.5, aes(x=`Number of Chromosomes`, y=`Unique Genes`, fill=factor("Unique Genes"), col=factor("Unique Genes")))+
  theme_classic()+
  labs(x="Number of Genomes", y="Number of Genes")+
  theme(legend.title = element_blank())

RoaryConservedAndNew <- ConservedGenes %>% ggplot(aes(x=`Number of Chromosomes`, y=`Conserved Genes`, fill=factor("Conserved Genes"), col=factor("Conserved Genes")))+
  geom_smooth()+
  geom_point(alpha=0.5)+
  geom_smooth(data = NewGenes, aes(x=`Number of Chromosomes`, y=`New Genes`, fill=factor("New Genes"), col=factor("New Genes")))+
  geom_point(data = NewGenes, alpha=0.5, aes(x=`Number of Chromosomes`, y=`New Genes`, fill=factor("New Genes"), col=factor("New Genes")))+
  theme_classic()+
  labs(x="Number of Genomes", y="Number of Genes")+
  theme(legend.title = element_blank())

RoaryEstimates <- ggarrange(RoaryConservedAndNew, RoaryTotalAndUnique, common.legend = F, legend = "bottom")

#Prep Plot of Core, Shell, Cloud, and Accessory genes

RoarySummary <- read.table("./Roary on data 181, data 179, and others Summary statistics.tabular", sep = "\t")
RoarySummary <- RoarySummary[1:4,]
RoarySummary$V1=paste(RoarySummary$V1, RoarySummary$V2, "n =", RoarySummary$V3)
RoaryPlot <- RoarySummary %>% ggplot(aes(x="", fill=V1, y=(V3)))+
  geom_bar(stat="identity", col="black")+
  coord_polar("y")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank())
RoaryPlot

RoaryForThesis <- ggarrange(RoaryEstimates, RoaryPlot, ncol = 1)
ggsave("D:/Treework/20221025AutoMLSTRun/Roary/Figures/Withoutpataloonsplit.png", plot = RoaryForThesis, dpi = 300, width = 10)
