# TODO: Add comment
# 
# Author: spapillo
###############################################################################



intersectCNV <- function(CNVformat, filter, type="both") {
	x <- unlist(sapply(1:length(CNVformat), function(i) {
						if(type == "both")
							used_CNVs <- filter[[i]]
						else if(type == "gain")
							used_CNVs <- CNVformat[[i]][,'logratio'] > 0 & filter[[i]]
						else if(type == "loss")
							used_CNVs <- CNVformat[[i]][,'logratio'] < 0 & filter[[i]]
						unique(unlist(strsplit(x=CNVformat[[i]][used_CNVs,'genes'], ";")))
					}))
	return(sort(table(x), decreasing=T))
}