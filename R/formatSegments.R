formatSegments <- function(CNAobject, sample_intensity, control_intensity, site_annotation, p.adjust.method="bonferroni", verbose=T, plotting=T) {
	all_segments <- CNAobject$output

	segments_per_sample <- lapply(unique(all_segments$ID), function(sample) all_segments[all_segments$ID == sample,-1])
	names(segments_per_sample) <- unique(all_segments$ID)
	
	if(is.vector(sample_intensity))
		sample_intensity <- as.matrix(sample_intensity)
	# Data check
	if(length(segments_per_sample) != ncol(sample_intensity))
		stop("The number of samples in CNAobject must match ncol(sample_intensity)")
	
	if(nrow(sample_intensity) != nrow(control_intensity))
		stop("sample_intensity and control_intensity must have the same number of rows")
	
	if(is.null(rownames(site_annotation)))
		rownames(site_annotation) <- as.character(site_annotation$ID)

	# Ensure row consistent ordering
	control_intensity <- control_intensity[rownames(sample_intensity), ]
	site_annotation <- site_annotation[rownames(sample_intensity), ]

	x <- lapply(1:length(segments_per_sample), function(i) {
		if(plotting)
			pdf(paste(names(segments_per_sample)[i], ".pdf", sep=""))
		result <- t(apply(segments_per_sample[[i]], 1, function(cnv) {
				# Extract probes
				probes <- site_annotation$CHR == cnv['chrom'] &
					as.numeric(site_annotation$MAPINFO) >= as.numeric(cnv['loc.start']) &
					as.numeric(site_annotation$MAPINFO) <= as.numeric(cnv['loc.end'])
									
				# Compute segment values
				control_int_sum <- colSums(control_intensity[probes,])
				control_mean <- mean(control_int_sum)
				control_sd <- sd(control_int_sum)
				sample_sum <- sum(sample_intensity[probes,i])
									
				# Compute p-value (2 sided t-test)
				z_score <- (sample_sum - control_mean) / control_sd
				p_value <-  2 * pnorm(-abs(z_score))
				if(plotting)
				{
					window <- range(c(control_int_sum, sample_sum))
					title <- paste("chr ", cnv['chrom'], " :  ", cnv['loc.start'], "-", cnv['loc.end'], sep="")
					hist(control_int_sum, main=title, xlim=window, breaks=20)
					abline(v=sample_sum, lty=3, lwd=2, col="red")
				}

				# Extract genes in the segment
				genes <- as.character(site_annotation$UCSC_RefGene_Name[probes])
				genes <- unique(unlist(strsplit(x=genes, ";")))
				genes <- paste(genes, collapse=";")
				l <- as.numeric(cnv['loc.end']) - as.numeric(cnv['loc.start'])
				return(c(cnv, seg.length=l, pvalue=p_value, genes=genes, ctrl.mean=control_mean, ctrl.sd=control_sd, sample.value=sample_sum, z=z_score))
			}))

		if(plotting)
			dev.off()
		# pvalue correction for multiple testing
		result <- as.data.frame(result, stringsAsFactors=F)
		result$adjusted.pvalue <- p.adjust(result$pvalue, method=p.adjust.method)

		if(verbose)
			message(paste("Processed", names(segments_per_sample)[i]))
		return(result)
	})

	names(x) <- names(segments_per_sample)
	return(x)
}
