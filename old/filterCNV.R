# TODO: Add comment
# 
# Author: spapillo
###############################################################################


filterCNV <- function(CNVformat, tick.threshold=50, pvalue.threshold=0.01, breakpoints=c(0,0,0,0)) {
	lapply(CNVformat, function(sample) as.numeric(sample[,'num.mark']) >= tick.threshold &
						as.numeric(sample[,'adjusted.pvalue']) <= pvalue.threshold)
}

