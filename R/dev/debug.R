library(minfi)

source("~/git/CopyNumber450k/R/dev/MethylCNVDataSet.R")

load("~/git/CopyNumber450k/data/control_RGset.RData")
pData(control_RGset)$Sample_Group[1] <- "case"

mcds <- new("MethylCNVDataSet", control_RGset) 
