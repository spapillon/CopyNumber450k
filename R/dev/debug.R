library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)

# load("~/git/CopyNumber450k/data/control_RGset.RData")
# pData(control_RGset)$Sample_Group[1] <- "case"

source("~/git/CopyNumber450k/R/dev/MethylCNVDataSet.R")
source("~/git/CopyNumber450k/R/dev/FromRGChannelSet.R")

mcds <- MethylCNVDataSetFromRGChannelSet(control_RGset)
