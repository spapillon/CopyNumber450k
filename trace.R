# TODO: Add comment
# 
# Author: spapillo
###############################################################################


# load everything (RGSet)
library(minfi)
source('~/git/CopyNumber450k/R/generics.R')
source('~/git/CopyNumber450k/R/extractFromRGSet450k.R')
source('~/git/CopyNumber450k/R/CNVObject.R')

source('~/git/CopyNumber450k/FunNormCN450k.R')


path <- '~/Documents/iChange/data_ETMR'
RGset <- read.450k.exp(base = path, targets = read.450k.sheet(path))

CNVobj <- CNVObject(RGset)

CNVobj <- filterSNPProbes(CNVobj)

CNVobj <- normalize(CNVobj)

CNVobj <- buildSegments(CNVobj)

CNVobj <- createFilters(CNVobj)