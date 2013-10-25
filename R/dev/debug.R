################################################################################ 

# library(minfi)

# library(IlluminaHumanMethylation450kmanifest)

# library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)

################################################################################ 

# load('~/git/CopyNumber450k/data/control_RGset.RData')

# pData(control_RGset)$Sample_Group[1] <- 'case'

# updateObject(control_RGset)

################################################################################ 

source("~/git/CopyNumber450k/R/dev/extract.R")
source("~/git/CopyNumber450k/R/dev/normalization.functional.R")
source("~/git/CopyNumber450k/R/dev/normalization.quantile.R")

source("~/git/CopyNumber450k/R/dev/class.R")
source("~/git/CopyNumber450k/R/dev/methods.R")
source("~/git/CopyNumber450k/R/dev/comparison.R")
source("~/git/CopyNumber450k/R/dev/plots.R")

################################################################################ 

mcds <- MethylCNVDataSetFromRGChannelSet(control_RGset)
pData(mcds)$sex <- predictSampleSexes(mcds)

mcds.f <- normalize(mcds, "functional")
mcds.q <- normalize(mcds, "quantile")

################################################################################  
