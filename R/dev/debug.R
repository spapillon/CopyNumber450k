library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)

# load('~/git/CopyNumber450k/data/control_RGset.RData')
# pData(control_RGset)$Sample_Group[1] <- 'case'

source("~/git/CopyNumber450k/R/dev/extractFromRGChannelSet.R")
source("~/git/CopyNumber450k/R/dev/normalization.functional.R")
source("~/git/CopyNumber450k/R/dev/normalization.quantile.R")

source("~/git/CopyNumber450k/R/dev/MethylCNVDataSet.R")
source("~/git/CopyNumber450k/R/dev/MethylCNVDataSet.transformation.R")


mcds <- MethylCNVDataSetFromRGChannelSet(control_RGset)
pData(mcds)$sex <- predictSampleSexes(mcds)

mcds.f <- normalize(mcds, "functional")
mcds.q <- normalize(mcds, "quantile") 
