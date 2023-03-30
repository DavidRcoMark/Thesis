#Intro----
#Author: David Mark
#Purpose: Generates a csv for clinker to group genes by antiSMASH predicted role

#Import data----
library(genbankr) #For manipulating gbks
library(tidyverse)
setwd("Dir of antiSMASH annotated GenBank files")
Genbanks <- list.files(pattern = ".gbk", all.files = T)
ParsedGenbanks <- lapply(Genbanks, parseGenBank)

#Process imported Genbanks ----
SplitGenbanks <- unlist(ParsedGenbanks, recursive = F)
ParsedGenbanks <- NULL
FeatureIndices <- grep(x=names(SplitGenbanks), pattern = "FEATURES")
justfeatures <-SplitGenbanks[FeatureIndices]
justfeatures <- unlist(justfeatures, recursive = F)

#Coerce genes of interest into data frame ----
GenesToBeColoured <- grep(x=justfeatures, pattern = "gene_kind")
BiosyntheticFeatures <- justfeatures[GenesToBeColoured]
names(BiosyntheticFeatures) = seq(length(names(BiosyntheticFeatures)))
gene_kind <- (BiosyntheticFeatures %>% imap_dfr("gene_kind"))
locus_tag <- (BiosyntheticFeatures %>% imap_dfr("locus_tag"))

#Spits out data frame into csv for clinker ----
forclinker <- data.frame("locus_tag"=unlist(locus_tag),
           "gene_kind"=unlist(gene_kind))
write.table(x=forclinker,
            file = "./forclinker.csv",
            col.names = FALSE, row.names =  F,
            sep=",")
