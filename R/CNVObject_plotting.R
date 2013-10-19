# TODO: Add comment
# 
# Author: spapillo
###############################################################################



setMethod("plot", signature("CNVObject"), function(x,  y="missing", path=".") {
	segments_list <- segments(x)
	filters_list <- filters(x)
	template_sample <- segments_list[[1]]
	chromosomes <-  unique(template_sample[,'chrom'])
	idx <- which(chromosomes %in% c("X", "Y"))
	chromosomes <- sort(as.numeric(chromosomes[-idx]))
	site_per_chr <- cumsum(c(0, sapply(chromosomes, function(chr) max(as.numeric(template_sample[template_sample[,'chrom']== chr,'loc.end'])))))
	offset <- site_per_chr - min(as.numeric(template_sample[template_sample[,'chrom'] == 1, 'loc.start']))
	Xmax <- max(site_per_chr)
			
	lapply((1:length(segments_list)), function(i) {
				png(paste(path, "/", names(segments_list)[i], ".png", sep=""), height=900, width=1200)				
				sample <- segments_list[[i]]
				Ymin <- min(c(-1, as.numeric(sample[filters_list[[i]],'seg.mean'])))
				Ymax <- max(c(1, as.numeric(sample[filters_list[[i]],'seg.mean'])))
				plot(range(0, Xmax), range(Ymin, Ymax), type='n', xaxt='n', xlab="", ylab="", main=names(segments_list[i]))
				xlabs <- sapply(2:length(site_per_chr), function(j) (site_per_chr[j] - site_per_chr[(j-1)]) / 2 + site_per_chr[(j-1)])
				axis(1, at=xlabs, labels=chromosomes, lty=0)
				abline(v=site_per_chr, lty=3)
				
				lapply(1:length(chromosomes), function(j) {
							used_segments <- sample[,'chrom'] == chromosomes[j]
							colors <- ifelse(filters_list[[i]][used_segments], "red", "black")		
							starts <- as.numeric(sample[used_segments,'loc.start']) +  offset[j]
							ends <- as.numeric(sample[used_segments,'loc.end']) + offset[j]
							y <- as.numeric(sample[used_segments,'seg.mean'])
							graphics::segments(starts, y, ends, y, col=colors, lwd=2, lty=1)
						})
				dev.off()
			})
})

setMethod("plotSex", signature("CNVObject"), function(object) {
	cnQuantiles <- object@RGSetSummary$cnQuantiles
	plot(log2(cnQuantiles$Y[250, ]) - log2(cnQuantiles$X[250, ]), rep(0,length(cnQuantiles$Y[250, ])))
	abline(v=-3, lty=3, col="red")
	text(x=-4, y=0.5, label="Female", col="red")
	text(x=-1, y=0.5, label="Male", col="red")
})
