#Generates multiple sequence alignments of conserved regulators ----

library(Biostrings)
library(msa)
library(ggmsa)
library(ggplot2)

#Import Data ----
setwd("D:/BLAST/AASeqs/")
O3 <- readAAStringSet("./O3Proteome.fasta")
PH63 <- readAAStringSet("./PH63Proteome.fasta")
O5 <- readAAStringSet("./O5Proteome.fasta")
ScoRegs <- readAAStringSet("./CoelicolorRegulators.fasta")

PH63@ranges@NAMES <- gsub(pattern = "PNLHFBPO", replacement = "PH63", x=PH63@ranges@NAMES)       
O3@ranges@NAMES <- gsub(pattern = "NIMBDPOJ", replacement = "O3", x=O3@ranges@NAMES)       
O5@ranges@NAMES <- gsub(pattern = "BFMPOHDH", replacement = "O5", x=O5@ranges@NAMES) 

#Generate MSAs and Plots ----
##AdpA ----
AdpASet <- c(ScoRegs[3], 
             O5[grep(pattern = "O5_06151", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_01650", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_00771", x=PH63@ranges@NAMES)])
AdpAAligned <- msa(AdpASet)
UnmaskedAdpAAlignment <- unmasked(AdpAAligned)
AdpAPlot <- ggmsa(UnmaskedAdpAAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##Amfr ----
AmfRSet <- c(ScoRegs[4], 
             O5[grep(pattern = "O5_05566", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_01047", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_05028", x=PH63@ranges@NAMES)])
AmfRAligned <- msa(AmfRSet)
UnmaskedAmfRAlignment <- unmasked(AmfRAligned)
AmfRPlot <- ggmsa(UnmaskedAmfRAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##BldD ----
BldDSet <- c(ScoRegs[2], 
                        O5[grep(pattern = "O5_02268", x=O5@ranges@NAMES)],
                        O3[grep(pattern = "O3_04205", x=O3@ranges@NAMES)],
                        PH63[grep(pattern = "PH63_02127", x=PH63@ranges@NAMES)])
BldDAligned <- msa(BldDSet)
UnmaskedBldDAlignment <- unmasked(BldDAligned)
BldDPlot <- ggmsa(UnmaskedBldDAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##BldM ----
BldMSet <- c(ScoRegs[6], 
             O5[grep(pattern = "O5_01751", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_00252", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_04197", x=PH63@ranges@NAMES)])
BldMAligned <- msa(BldMSet)
UnmaskedBldMAlignment <- unmasked(BldMAligned)
BldMPlot <- ggmsa(UnmaskedBldMAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##BldN ----
BldNSet <- c(ScoRegs[1], 
             O5[grep(pattern = "O5_05918", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_01384", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_05370", x=PH63@ranges@NAMES)])
BldNAligned <- msa(BldNSet)
UnmaskedBldNAlignment <- unmasked(BldNAligned)
BldNPlot <- ggmsa(UnmaskedBldNAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()
BldNPlot

##RsbN ----
RsbNSet <- c(ScoRegs[5], 
             O5[grep(pattern = "O5_05918", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_01384", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_05370", x=PH63@ranges@NAMES)])
RsbNAligned <- msa(RsbNSet)
UnmaskedRsbNAlignment <- unmasked(RsbNAligned)
RsbNPlot <- ggmsa(UnmaskedRsbNAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##WhiA ----
WhiASet <- c(ScoRegs[7], 
             O5[grep(pattern = "O5_04547", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_00099", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_04045", x=PH63@ranges@NAMES)])
WhiAAligned <- msa(WhiASet)
UnmaskedWhiAAlignment <- unmasked(WhiAAligned)
WhiAPlot <- ggmsa(UnmaskedWhiAAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##WhiB ----
WhiBSet <- c(ScoRegs[8], 
             O5[grep(pattern = "O5_01069", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_00968", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_01189", x=PH63@ranges@NAMES)])
WhiBAligned <- msa(WhiBSet)
UnmaskedWhiBAlignment <- unmasked(WhiBAligned)
WhiBPlot <- ggmsa(UnmaskedWhiBAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##WhiG ----
WhiGSet <- c(ScoRegs[9], 
             O5[grep(pattern = "O5_03987", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_04462", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_03593", x=PH63@ranges@NAMES)])
WhiGAligned <- msa(WhiGSet)
UnmaskedWhiGAlignment <- unmasked(WhiGAligned)
WhiGPlot <- ggmsa(UnmaskedWhiGAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

##WhiI ----
WhiISet <- c(ScoRegs[10], 
             O5[grep(pattern = "O5_03577", x=O5@ranges@NAMES)],
             O3[grep(pattern = "O3_01138", x=O3@ranges@NAMES)],
             PH63[grep(pattern = "PH63_04260", x=PH63@ranges@NAMES)])
WhiIAligned <- msa(WhiISet)
UnmaskedWhiIAlignment <- unmasked(WhiIAligned)
WhiIPlot <- ggmsa(UnmaskedWhiIAlignment, 
                  seq_name = T,
                  char_width=0,
                  show.legend = T,
                  color = "Chemistry_AA")+
  geom_msaBar()

#Save Plots ----
ggsave("D:/Figures/AdpAAlignment.png",
       plot = AdpAPlot,
       width = 40,
       height = 2,
       dpi=300)

ggsave("D:/Figures/AmfRAlignment.png",
       plot = AmfRPlot,
       width = 42,
       height = 2,
       dpi=300)

ggsave("D:/Figures/BldDAlignment.png",
       plot = BldDPlot,
       width = 18,
       height = 2,
       dpi=300)

ggsave("D:/Figures/BldMAlignment.png",
       plot = BldMPlot,
       width = 24,
       height = 2,
       dpi=300)

ggsave("D:/Figures/BldNAlignment.png",
       plot = BldNPlot,
       width = 33,
       height = 2,
       dpi=300)

ggsave("D:/Figures/RsbNAlignment.png",
       plot = RsbNPlot,
       width = 43,
       height = 2,
       dpi=300)

ggsave("D:/Figures/WhiAAlignment.png",
       plot = WhiAPlot,
       width = 33,
       height = 2,
       dpi=300)

ggsave("D:/Figures/WhiBAlignment.png",
       plot = WhiBPlot,
       width = 12,
       height = 2,
       dpi=300)

ggsave("D:/Figures/WhiGAlignment.png",
       plot = WhiGPlot,
       width = 29,
       height = 2,
       dpi=300)

ggsave("D:/Figures/WhiIAlignment.png",
       plot = WhiIPlot,
       width = 29,
       height = 2,
       dpi=300)