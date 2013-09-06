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
mset <- preprocessRaw(RGset)
cnMatrix <- getMeth(mset) + getUnmeth(mset)
sexes <- predictSex(extractedData = extracted_data, cutoff = 0)
cnMatrix.norm <- normalizeFunNorm450kCN(cnMatrix = cnMatrix, extractedData = extracted_data, predictedSex = sexes)

# filter probes (verify utility)
con.intensity.median <- rowMedians(cnMatrix[, phenotype == "control"])
median.quantiles <- quantile(con.intensity.median, c(0.05, 0.85), na.rm=TRUE)

con.intensity.mad <- apply(cnMatrix[, phenotype == "control"], 1, mad, na.rm=TRUE)
mad.quantiles <- quantile(con.intensity.mad, 0.8, na.rm=TRUE)

usable_probes <- con.intensity.median >= median.quantiles[1] &
	con.intensity.median <= median.quantiles[2] &
	con.intensity.mad <= mad.quantiles

# log2
# segmentation
# formatting
# visualize
