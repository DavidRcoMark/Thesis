#Intro ----
#Author: David Mark
#Purpose: Draws Heatmap of BLASTp Results
#Contributes to: Figure 4.3.

library(tidyverse)
library(ggnewscale)
library(msa)

#Import BLASTp Results ----
BLAST <- read.table("F:/BLAST/BlastResults/BlastP.tabular")
colnames(BLAST)=c("qseqid",
                  "sseqid",
                  "pident",
                  "length",
                  "mismatch",
                  "gapopen",
                  "qstart",
                  "qend",
                  "sstart",
                  "send",
                  "evalue",
                  "bitscore")

#Tidy BLASTp Results ----
BLAST <- BLAST %>% separate(col = sseqid, into =c("Org", "CDS"), sep = "_", remove = F)
BLAST$Org <- gsub(pattern = "PNLHFBPO", replacement = "PH63", x=BLAST$Org)       
BLAST$Org <- gsub(pattern = "NIMBDPOJ", replacement = "O3", x=BLAST$Org)       
BLAST$Org <- gsub(pattern = "BFMPOHDH", replacement = "O5", x=BLAST$Org)       

BestHits <- BLAST %>% group_by(qseqid, Org) %>% summarise(sseqid, Org, evalue, pident) %>% filter(evalue==min(evalue))
      
BestHits$sseqid <- gsub(pattern = "PNLHFBPO", replacement = "PH63", x=BestHits$sseqid)
BestHits$sseqid <- gsub(pattern = "NIMBDPOJ", replacement = "O3", x=BestHits$sseqid)       
BestHits$sseqid<- gsub(pattern = "BFMPOHDH", replacement = "O5", x=BestHits$sseqid) 

#Draw heatmap----
BlastMap <- BestHits %>% ggplot(aes(x=qseqid, y=Org, fill=pident))+
  geom_tile()+
  theme_classic()+
  scale_fill_viridis_c(option = "D")+
  labs(x="Developmental Regulator", y="Organism", fill="Identity (%)")+
  scale_y_discrete(labels=c("M. sp. O3", "M. sp. O5", "M. sp. PH63"))+
  new_scale_fill()+
  geom_label(aes(label=pident, fill=pident>70), alpha=1, show.legend = F)
BlastMap

ggsave(BlastMap, filename = "D:/Figures/Blastmap.png", dpi=300)
