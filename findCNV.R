# TODO: Add comment
# 
# Author: spapillo
###############################################################################


findCNV <- function(CNVformat, CNVfilter, CNVs, type="both") {
	
	x <- sapply(1:length(CNVformat), function(i) {
				if(type == "both")
					used_CNVs <- filter[[i]]
				else if(type == "gain")
					used_CNVs <- CNVformat[[i]][,'logratio'] > 0 & filter[[i]]
				else if(type == "loss")
					used_CNVs <- CNVformat[[i]][,'logratio'] < 0 & filter[[i]]
				
				genes <- unique(unlist(strsplit(x=CNVformat[[i]][used_CNVs,'genes'], ";")))
				return(CNVs %in% genes)
			})

	ifelse(is.matrix(x), x <- t(x), x <- as.matrix(x))
	rownames(x) <- names(CNVformat)
	colnames(x) <- CNVs
	return(x)
}
