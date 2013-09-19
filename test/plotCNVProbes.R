# TODO: Add comment
# 
# Author: spapillo
###############################################################################

plotCNVProbes <- function(intensities, site_annotation) {
	chr_colors <- c("black","green")
	chromosomes <- unique(site_annotation$CHR)
	chromosomes <- sort(as.numeric(chromosomes[!chromosomes %in% c("X","Y")]))
	site_per_chr <- cumsum(c(0, sapply(chromosomes, function(chr) max(as.numeric(site_annotation$MAPINFO[site_annotation$CHR == chr])))))
	
	Xmax <- max(site_per_chr)
	Ymin <- min(intensities)
	Ymax <- max(intensities)
	
	plot(range(0, Xmax), range(Ymin, Ymax), type='n', xaxt='n', xlab="", ylab="", main="")
	
	lapply(1:length(chromosomes), function(i) {
				used_probes <- site_annotation$CHR == chromosomes[i]
				coord <- as.numeric(site_annotation$MAPINFO[used_probes]) + site_per_chr[i]
				values <- intensities[used_probes]
				lines(coord, values, type='p', pch=20, col=chr_colors[(i%%2 + 1)], cex=0.5)
			})
	
}
