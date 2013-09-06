# TODO: Add comment
# 
# Author: spapillo
###############################################################################


# load everything (RGSet)
library(minfi)
source('~/git/CopyNumber450k/formatCNVobject.R')
source('~/git/CopyNumber450k/extractFromRGSet450k.R')
source('~/git/CopyNumber450k/FunNormCN450k.R')
source('~/git/CopyNumber450k/createCNVObject.R')


path <- '~/Documents/iChange/data_GBM'
RGset <- read.450k.exp(base = path, targets = read.450k.sheet(path))

# normalize (JP's new thing)
extracted_data <- extractFromRGSet450k(RGset)
mset <- preprocessRaw(RGset)
cnMatrix <- getMeth(mset) + getUnmeth(mset)
sexes <- predictSex(extractedData = extracted_data, cutoff = 0)
cnMatrix.norm <- normalizeFunNorm450kCN(cnMatrix = cnMatrix, extractedData = extracted_data, predictedSex = sexes)

rm(sexes)

# define phenotypes

# filter probes (verify utility) based on variance
con.intensity.median <- rowMedians(cnMatrix[, phenotype == "control"])
median.quantiles <- quantile(con.intensity.median, c(0.05, 0.85), na.rm=TRUE)

con.intensity.mad <- apply(cnMatrix[, phenotype == "control"], 1, mad, na.rm=TRUE)
mad.quantiles <- quantile(con.intensity.mad, 0.8, na.rm=TRUE)

usable_probes <- con.intensity.median >= median.quantiles[1] &
	con.intensity.median <= median.quantiles[2] &
	con.intensity.mad <= mad.quantiles

rm(con.intensity.media, median.quantiles, con.intensity.mad, mad.quantiles)

# filter SNP probes
require(GEOquery)
site_annotation = Table(getGEO('GPL13534') )
rownames(site_annotation) <- site_annotation$ID
site_annotation <- site_annotation[rownames(cnMatrix.norm), ]
non_snp_probes <- site_annotation$Probe_SNPs == "" & site_annotation$Probe_SNPs_10  == ""

good_probes <- usable_probes & non_snp_probes
rm(usable_probes, non_snp_probes)

# log2
CNVobject <- createCNVObject(cnMatrix, phenotype, good_probes)


# segmentation
# formatting
# visualize
