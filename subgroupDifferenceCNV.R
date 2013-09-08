# TODO: Add comment
# 
# Author: spapillo
###############################################################################


subgroupDifferenceCNV <- function(CNVformat, filter, subgroups, pheno)
{
	group1 <- pheno == subgroups[1] & !is.na(pheno)
	group2 <- pheno == subgroups[2] & !is.na(pheno)
	group1_size <- sum(group1)
	group2_size <- sum(group2)
	
	#Gains
	group1_CNVs <- intersectCNV(CNVformat[group1], filter[group1], type="gain")
	group2_CNVs <- intersectCNV(CNVformat[group2], filter[group2], type="gain")
	gains <- subgroupDifferenceCNVByType(group1_CNVs, group2_CNVs, group1_size, group2_size) 
	colnames(gains) <- c(subgroups[1], subgroups[2], "pvalue")
	
	#Losses
	group1_CNVs <- intersectCNV(CNVformat[group1], filter[group1], type="loss")
	group2_CNVs <- intersectCNV(CNVformat[group2], filter[group2], type="loss")
	losses <- subgroupDifferenceCNVByType(group1_CNVs, group2_CNVs, group1_size, group2_size) 
	colnames(losses) <- c(subgroups[1], subgroups[2], "pvalue")
	
	return(list(gains=gains, losses=losses))
}

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
