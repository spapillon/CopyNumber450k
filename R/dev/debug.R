library(minfi)

load("~/git/CopyNumber450k/data/control_RGset.RData")
pData(control_RGset)$Sample_Group[1] <- "case"

source("~/git/CopyNumber450k/R/dev/MethylCNVDataSet.R")

mcds <- new("MethylCNVDataSet", control_RGset) 
