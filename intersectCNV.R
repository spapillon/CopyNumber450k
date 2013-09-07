# TODO: Add comment
# 
# Author: spapillo
###############################################################################



CNVintersect <- function(CNVformat, filter) {
	x <- unlist(sapply(1:length(CNVformat), function(i) {
						unique(unlist(strsplit(x=CNVformat[[i]][filter[[i]],'genes'], ";")))
					}))
	return(sort(table(x), decreasing=T))
}