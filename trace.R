# TODO: Add comment
# 
# Author: spapillo
###############################################################################


# load everything (RGSet)
library(minfi)
source('formatCNVobject.R')
source('extractFromRGSet450k.R')
source('FunNormCN450k.R')

path <- '~/Documents/iChange/data_ETMR'
RGset <- read.450k.exp(base = path, targets = read.450k.sheet(path))
# normalize (JP's new thing)
extracted_data <- extractFromRGSet450k(RGset)

# log2
# segmentation
# formatting
# visualize
