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

#findCNV <- function(CNVformat, CNVfilter, CNVs, type="both") {
#	if(length(CNVs) != 1)
#		stop("CNVs has to be of length 1")
#	x <- lapply(1:length(CNVformat), function(i) {
#				if(type == "both")
#					used_CNVs <- filter[[i]]
#				else if(type == "gain")
#					used_CNVs <- CNVformat[[i]][,'logratio'] > 0 & filter[[i]]
#				else if(type == "loss")
#					used_CNVs <- CNVformat[[i]][,'logratio'] < 0 & filter[[i]]
#				if(sum(used_CNVs) == 0)
#					return(NULL)
#				cnv_genes <- strsplit(x=CNVformat[[i]][,'genes'], ";")
#				
#				found_CNVs <- sapply(cnv_genes, function(i) sum(i %in% CNVs) > 0)
#				CNVs %in% cnv_genes[found_CNVs]
#				return(CNVformat[[i]][used_CNVs & found_CNVs,c('chrom', 'loc.start', 'loc.end', 'num.mark','logratio', 'adjusted.pvalue')])
#				
#				#genes <- unique(unlist(strsplit(x=CNVformat[[i]][used_CNVs,'genes'], ";")))
#				
#				#return(CNVs %in% genes)
#			})
#	names(x) <- names(CNVformat)#
#	ifelse(is.matrix(x), x <- t(x), x <- as.matrix(x))
#	rownames(x) <- names(CNVformat)
#	colnames(x) <- CNVs
#	return(x)
#}
