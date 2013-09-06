# TODO: Add comment
# 
# Author: spapillo
###############################################################################


createCNVObject <- function(cnMatrix, phenotype, good_probes, site_annotation) {
	library(DNAcopy)
	cnMatrix <- cnMatrix[good_probes, ]
	site_annotation <- site_annotation[good_probes, ]
	control_medians <- rowMeans(cnMatrix[, phenotype == "control"])
	cases_log2 <- log2(cnMatrix[, phenotype != "control"] / control_medians)
	CNA.object <- CNA(cases_log2, ordered(site_annotation$CHR), as.numeric(site_annotation$MAPINFO), data.type="logratio" )
	
	smoothed.CNA.object <- smooth.CNA(CNA.object)
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width=5 ,verbose=1, nperm=10000, alpha=0.01, undo.splits="sdundo", undo.SD=2)
}