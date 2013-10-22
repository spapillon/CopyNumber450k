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


source('~/git/CopyNumber450k/R/extractFromRGSet450k.R')

source('~/git/CopyNumber450k/R/CNVObject.R')
source('~/git/CopyNumber450k/R/CNVObject_initialization.R')
source('~/git/CopyNumber450k/R/CNVObject_accession.R')
source('~/git/CopyNumber450k/R/CNVObject_replacement.R')
source('~/git/CopyNumber450k/R/CNVObject_comparison.R')
source('~/git/CopyNumber450k/R/CNVObject_plotting.R')

source('~/git/CopyNumber450k/R/formatSegments.R')
source('~/git/CopyNumber450k/R/functionalNormalization.R')
source('~/git/CopyNumber450k/R/quantileNormalization.R')
source('~/git/CopyNumber450k/R/subgroupDifferenceCNVByType.R')


path <- '~/Documents/iChange/data_ETMR'
load('~/git/CopyNumber450k/data/control_RGset.RData')

#spapillon-specific
#case_RGset <- read.450k.exp(base = path, targets = read.450k.sheet(path))
#RGset <- combine(control_RGset, case_RGset)

#ndejay-specific
pData(control_RGset)$Sample_Group[1:10] <- "case"
RGset <- control_RGset

CNVobj <- CNVObject(RGset)
sampleSexes(CNVobj) <- predictSex(CNVobj)
sampleGroups(CNVobj) <- pData(RGset)$Sample_Group
sampleNames(CNVobj) <- pData(RGset)$Sample_Name



#rm(RGset, control_RGset, case_RGset)

CNVobj <- filterSNPProbes(CNVobj)

CNVobj_norm <- normalize(CNVobj)
CNVobj_norm2 <- normalize(CNVobj, "quantile")

CNVobj <- buildSegments(CNVobj)
CNVobj_norm <- buildSegments(CNVobj_norm)
CNVobj_norm2 <- buildSegments(CNVobj_norm2)

CNVobj <- createFilters(CNVobj)
CNVobj_norm <- createFilters(CNVobj_norm)
CNVobj_norm2 <- createFilters(CNVobj_norm2)
plot(CNVobj, path='raw/')
plot(CNVobj_norm, path='funnorm/')
plot(CNVobj_norm2, path='quannorm/')

