#Author: David Mark
#Purpose: Calculate codon usage and extract TTA dependent genes (Fig. 4.5. and 4.6.)

library(coRdon)
library(tidyverse)
library(ggfortify)
setwd("D:/CordonForAtacamaMicros/")

#Prep O3 Codon Count----
O3 <- readSet(file = "./FFNs/O3CDS.fasta")
O3Codons <- codonTable(O3)
O3Counts <- codonCounts(O3Codons)
O3SumOfCounts <- data.frame(as.list(colSums(O3Counts)))
O3SumOfCounts <-
  data.frame("Count" = t(O3SumOfCounts), "Codon" = row.names(t(O3SumOfCounts)))
O3SumOfCounts <- O3SumOfCounts %>% mutate("Strain" = "M. sp. O3")

#Repeat for PH63----
PH63 <- readSet(file = "./FFNs/PH63CDS.fasta")
PH63Codons <- codonTable(PH63)
PH63Counts <- codonCounts(PH63Codons)
PH63SumOfCounts <- data.frame(as.list(colSums(PH63Counts)))
PH63SumOfCounts <-
  data.frame("Count" = t(PH63SumOfCounts), "Codon" = row.names(t(PH63SumOfCounts)))
PH63SumOfCounts <-
  PH63SumOfCounts %>% mutate("Strain" = "M. sp. PH63")

#Repeat for O5 ----
O5 <- readSet(file = "./FFNs/O5CDS.fasta")
O5Codons <- codonTable(O5)
O5Counts <- codonCounts(O5Codons)
O5SumOfCounts <- data.frame(as.list(colSums(O5Counts)))
O5SumOfCounts <-
  data.frame("Count" = t(O5SumOfCounts), "Codon" = row.names(t(O5SumOfCounts)))
O5SumOfCounts <- O5SumOfCounts %>% mutate("Strain" = "M. sp. O5")

#Prep Data frame for plotting ----
SumOfCounts <- O3SumOfCounts %>%
  bind_rows(O5SumOfCounts) %>%
  bind_rows(PH63SumOfCounts)

#Plot Codon Counts ----
CountPlot <-
  SumOfCounts %>% arrange(Codon) %>% ggplot(aes(x = Count, y = Codon, fill =
                                                  Strain)) +
  geom_col(col = "black", position = "dodge") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 45),
        legend.position = "top")
CountPlot

ggsave(CountPlot,
       filename = "./plots/CodonCounts.png",
       dpi = 300,
       height = 15)


#Prep Tables of TTA Containing Genes----
##O3----
O3TTAGenes = NULL
O3TTAGenes <-
  as.data.frame(O3Counts, row.names = O3Codons@ID) %>% filter(TTA > 0)

O3TTAGenes <-
  O3TTAGenes %>% slice(-c(
    grep(x = row.names(O3TTAGenes), pattern = "*tRNA-"),
    grep(x = row.names(O3TTAGenes), pattern = "*S ribosomal RNA*"))) %>%
  select(TTA)

##O5----
O5TTAGenes = NULL
O5TTAGenes <- as.data.frame(O5Counts, row.names = O5Codons@ID) %>%
  filter(TTA > 0)

O5TTAGenes <- O5TTAGenes %>%
  slice(-c(
    grep(x = row.names(O5TTAGenes), pattern = "*tRNA-"),
    grep(x = row.names(O5TTAGenes), pattern = "*S ribosomal RNA*"))) %>%
  select(TTA)


##PH63----
PH63TTAGenes = NULL

PH63TTAGenes <- as.data.frame(PH63Counts, row.names = PH63Codons@ID) %>%
  filter(TTA > 0)

PH63TTAGenes <- PH63TTAGenes %>%
  slice(-c(
    grep(x = row.names(PH63TTAGenes), pattern = "*tRNA-"),
    grep(x = row.names(PH63TTAGenes), pattern = "*S ribosomal RNA*"))) %>%
  select(TTA)

##Write tables of TTA containing genes ----
write.csv(O5TTAGenes, file = "./O5TTAs.csv")
write.csv(O3TTAGenes, file = "./O3TTAs.csv")
write.csv(PH63TTAGenes, file = "./PH63TTAs.csv")

#Query to confirm presence of bldA ----
O3Codons@ID[grep(x=O3Codons@ID, pattern = "tRNA-Leu")]
O5Codons@ID[grep(x=O5Codons@ID, pattern = "tRNA-Leu")]
PH63Codons@ID[grep(x=PH63Codons@ID, pattern = "tRNA-Leu")]

ngenes <- data.frame("Organism"=c("M. sp. O3",
                                  "M. sp. O5",
                                  "M. sp. PH63"),
                     "Number of TTA Genes"=c(length(O3TTAGenes$TTA),
                                             length(O5TTAGenes$TTA),
                                             length(PH63TTAGenes$TTA)))
nplot <- ngenes %>% ggplot(aes(x=Organism, fill=Organism, y=Number.of.TTA.Genes))+
                    geom_col(position="dodge", col="black")+
  theme_classic()+
  labs(y="Number of TTA Containing Genes")+
  theme(legend.position = "none")
nplot

ggsave(nplot, filename = "./plots/genecount.png", dpi=300)
