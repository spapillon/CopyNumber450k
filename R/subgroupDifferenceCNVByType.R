# TODO: Add comment
# 
# Author: spapillo
###############################################################################

subgroupDifferenceCNVByType <- function(group1_CNVs, group2_CNVs, group1_size, group2_size) 
{
	target_CNVs <- unique(c(names(group1_CNVs), names(group2_CNVs)))
	x <- t(sapply(target_CNVs, function(cnv) {
						if(cnv %in% names(group1_CNVs))
							group1_hit <- group1_CNVs[cnv]
						else
							group1_hit <- 0
						group1_fail <- group1_size - group1_hit
						if(cnv %in% names(group2_CNVs))
							group2_hit <- group2_CNVs[cnv]
						else
							group2_hit <- 0
						group2_fail <- group2_size - group2_hit
						pvalue <-  fisher.test(matrix(c(group1_hit, group2_hit, group1_fail , group2_fail), 2))$p.value
						return(c(group1_hit, group2_hit, pvalue))
					}))
	x <- x[order(x[, 3]), ]
	return(x)
}
