#Author: David Mark
#Purpose: Plot strand bias of complete Micromonospora chromosomes (Fig. 5.4.)

library(tidyverse)
library(ape)
library(ggpubr)
setwd(
  "D:/Actinos/Micromonosporaceae/Micromonospora/Micromonospora Assemblies/gffs/Prokka on collection 142_ gff/"
)
GFFs <- lapply(FUN = read.gff, X = list.files(pattern = ".gff"))

##Process GFFs ----

M28 <- GFFs[[1]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6644909)) #NMT for each gene midpoint = 100*gene midpoint/Chromosome size

Maur110 <- GFFs[[2]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7255264))

Maur27029 <- GFFs[[3]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7025559))

Maurni <- GFFs[[4]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6758600))

Mcarb <- GFFs[[5]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7528335))

Mchok <- GFFs[[6]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6897819))

Mcor <- GFFs[[7]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6929787))

Mcox <- GFFs[[8]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6769793))

Mcran <- GFFs[[9]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6839926))

Mechaur <- GFFs[[10]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7203233))

Mechfus <- GFFs[[11]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7002627))

Mechspo <- GFFs[[12]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7775586))

Mend <- GFFs[[13]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6586761))

MHM <- GFFs[[14]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7565212))

Mins <- GFFs[[15]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6708094))

Mkrab <- GFFs[[16]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7074838))

Mmar <- GFFs[[17]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6673976))

Mnara <- GFFs[[18]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6612057))

Mpur <- GFFs[[19]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6674234))

Mrif <- GFFs[[20]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7011369))

Msag <- GFFs[[21]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6930736))

Msiam <- GFFs[[22]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6250219))

Mterm <- GFFs[[23]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6717200))

Mtul <- GFFs[[24]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7228472))

Mvir <- GFFs[[25]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7073985))

MWMMA <- GFFs[[26]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6483419))

MWMMC <- GFFs[[27]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 6281187))

Mzam <- GFFs[[28]] %>% filter(type == "gene") %>%
  mutate("NMT" = (100 * (((
    start + end
  ) / 2)) / 7096032))

GFFlist <- list(
  M28,
  Maur110,
  Maur27029,
  Maurni,
  Mcarb,
  Mchok,
  Mcor,
  Mcox,
  Mcran,
  Mechaur,
  Mechfus,
  Mechspo,
  Mend,
  MHM,
  Mins,
  Mkrab,
  Mmar,
  Mnara,
  Mpur,
  Mrif,
  Msag,
  Msiam,
  Mterm,
  Mtul,
  Mvir,
  MWMMA,
  MWMMC,
  Mzam
)

BoundGFFs <- bind_rows(GFFlist)
TestGFF <- BoundGFFs

TestGFF <-
  TestGFF %>%
  mutate("Midpoint" = NMT) %>%
  mutate("NMT" = ifelse(NMT > 50, NMT - 50, NMT + 50))

write.table(file = "./boundgff.csv", x=TestGFF, sep = ",", )
#Assign leading and lagging strands----

RightReplichore <- TestGFF %>%
  filter(NMT > 50)
LeftReplichore <- TestGFF %>%
  filter(NMT < 50)

LeftReplichore$strand <- gsub("-", "Leading", LeftReplichore$strand)
LeftReplichore$strand <-gsub("+", "Lagging", LeftReplichore$strand, fixed = T)

LeftReplichore %>%
  group_by(strand) %>%
  count()

RightReplichore$strand <-gsub("+", "Leading", RightReplichore$strand, fixed = T)
RightReplichore$strand <-gsub("-", "Lagging", RightReplichore$strand, fixed = T)
RightReplichore %>%
  group_by(strand) %>%
  count()

StrandCorrected <- bind_rows(LeftReplichore, RightReplichore)
StrandCorrected %>%
  group_by(strand) %>%
  count()

#Kernel Density Estimates ----

LeadLagPlot <-
  ggplot(StrandCorrected, aes(x = Midpoint, fill = strand)) +
  geom_density(adjust = 0.25, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 50, col = "red") +
  labs(x = "Position", y = "Gene Density")

LeadLagPlotByOrg <-
  ggplot(StrandCorrected, aes(x = Midpoint, fill = strand)) +
  geom_density(adjust = 0.25, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap( ~ seqid, ncol = 7) +
  geom_vline(xintercept = 50, col = "red") +
  labs(x = "Position", y = "Gene Density")

LeadLagBoth <-
  ggarrange(LeadLagPlotByOrg,
            LeadLagPlot,
            common.legend = T,
            ncol = 1)

PlusMinusPlot <- ggplot(TestGFF, aes(x = NMT, fill = strand)) +
  geom_density(adjust = 0.25, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 50, col = "red") +
  labs(x = "Position", y = "Gene Density")

PlusMinusPlotByOrg <- ggplot(TestGFF, aes(x = NMT, fill = strand)) +
  geom_density(adjust = 0.25, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 50, col = "red") +
  facet_wrap( ~ seqid, ncol = 7) +
  labs(x = "Position", y = "Gene Density")+
  ggtitle("A")

StrandBoth <-
  ggarrange(PlusMinusPlotByOrg,
            PlusMinusPlot,
            common.legend = T,
            ncol = 1)

ggsave(
  LeadLagBoth,
  filename = "D:/Figures/LeadLagPlotMidPoint.png",
  dpi = 300,
  width = 20,
  height = 15
)

ggsave(
  StrandBoth,
  filename = "D:/Figures/PlusMinusPlotMidPoint.png",
  dpi = 300,
  width = 20,
  height = 15
)

#T test of strands

Leading <- StrandCorrected %>%
  group_by(seqid) %>%
  count(strand) %>% filter(strand=="Leading")
Lagging <- StrandCorrected %>%
  group_by(seqid) %>%
  count(strand) %>% filter(strand=="Lagging")

LeadLagT <- t.test(Lagging$n, Leading$n)

Pluscount <- TestGFF %>%
  group_by(seqid) %>%
  count(strand) %>% filter(strand=="+")

Minuscount <- TestGFF %>%
  group_by(seqid) %>%
  count(strand) %>% filter(strand=="-")

PlusMinusT <- t.test(Pluscount$n, Minuscount$n)
#Draw Boxplots of Gene Count ----

PlusMinusBoxplot <- TestGFF %>%
  group_by(seqid) %>%
  count(strand) %>%
  ggplot(aes(x = strand, y = n, fill = factor(strand))) +
  geom_boxplot() +
  geom_line(aes(x = strand, group = seqid))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="Strand", y="Number of Genes")+
  ylim(c(2000,4000))+
  geom_text(aes(x="+", y=3500), label="NS")

LeadLagBoxplot <- StrandCorrected %>%
  group_by(seqid) %>%
  count(strand) %>%
  ggplot(aes(x = strand, y = n, fill = factor(strand))) +
  geom_boxplot() +
  geom_line(aes(x = strand, group = seqid))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="Strand", y="Number of Genes")+
  ylim(c(2000,4000))+
  geom_text(aes(x="Leading", y=3900), label="P=1.7003e-30")

BothBoxes <- ggarrange(PlusMinusBoxplot, LeadLagBoxplot)
ggsave(plot = BothBoxes, filename = "D:/Figures/Strandbias.png", dpi=300)
