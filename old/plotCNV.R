# TODO: Add comment
# 
# Author: spapillo
###############################################################################


plotCNV <- function(CNVformat, CNVfilter) {
	lapply(1:length(CNVformat), function(i) {
				png(paste(names(CNVformat)[i], ".png", sep=""), height=900, width=1200)
				
				sample <- CNVformat[[i]]
				chromosomes <-  sort(as.numeric(unique(CNVformat[[1]][,'chrom'])[-c(23,24)]))
				site_per_chr <- cumsum(c(0, sapply(chromosomes, function(chr) max(as.numeric(sample[sample[,'chrom']== chr,'loc.end'])))))
				offset <- site_per_chr - min(as.numeric(sample[sample[,'chrom'] == 1, 'loc.start']))
				Xmax <- max(site_per_chr)
				Ymin <- min(as.numeric(sample[,'logratio']))
				Ymax <- max(as.numeric(sample[,'logratio']))
				par(mfrow=c(2,1))
				plotCNVProbes(cases_log2[,i], site_annotation)
				plot(range(0, Xmax), range(Ymin, Ymax), type='n', xaxt='n', xlab="", ylab="", main=names(CNVformat[i]))
				xlabs <- sapply(2:length(site_per_chr), function(j) (site_per_chr[j] - site_per_chr[(j-1)]) / 2 + site_per_chr[(j-1)])
				axis(1, at=xlabs, labels=chromosomes, lty=0)
				abline(v=site_per_chr, lty=3)
				lapply(1:length(chromosomes), function(j) {
							used_segments <- sample[,'chrom'] == chromosomes[j]
							colors <- ifelse(CNVfilter[[i]][used_segments], "red", "black")
							
							starts <- as.numeric(sample[used_segments,'loc.start']) +  offset[j]
							ends <- as.numeric(sample[used_segments,'loc.end']) + offset[j]
							y <- as.numeric(sample[used_segments,'logratio'])
							segments(starts, y, ends, y, col=colors, lwd=2, lty=1)
						})
				dev.off()
			})
}
