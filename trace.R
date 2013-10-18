# TODO: Add comment
# 
# Author: spapillo
###############################################################################


# load everything (RGSet)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
library(matrixStats)
library(gmodels)
library(preprocessCore) 

source('~/git/CopyNumber450k/R/generics.R')
source('~/git/CopyNumber450k/R/extractFromRGSet450k.R')
source('~/git/CopyNumber450k/R/CNVObject.R')
source('~/git/CopyNumber450k/R/formatSegments.R')
source('~/git/CopyNumber450k/R/FunNormCN450k.R')


path <- '~/Documents/iChange/data_ETMR'
load('~/git/CopyNumber450k/data/control_RGset.RData')
case_RGset <- read.450k.exp(base = path, targets = read.450k.sheet(path))
RGset <- combine(control_RGset, case_RGset)

CNVobj <- CNVObject(RGset)
rm(RGset, control_RGset, case_RGset)

CNVobj <- filterSNPProbes(CNVobj)

CNVobj <- normalize(CNVobj)

CNVobj <- buildSegments(CNVobj)

CNVobj <- createFilters(CNVobj)

plot(CNVobj)
