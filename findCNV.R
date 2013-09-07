# TODO: Add comment
# 
# Author: spapillo
###############################################################################


findCNV <- function(CNVformat, CNVfilter, CNVs) {
	
	x <- sapply(1:length(CNVformat), function(i) {
						genes <- unique(unlist(strsplit(x=CNVformat[[i]][filter[[i]],'genes'], ";")))
						return(CNVs %in% genes)
					})

	ifelse(is.matrix(x), x <- t(x), x <- as.matrix(x))
	browser()
	rownames(x) <- names(CNVformat)
	colnames(x) <- CNVs
	return(x)
}
