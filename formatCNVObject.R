formatCNVObject <- function(CNAobject, sample_intensity, control_intensity, site_annotation, p.adjust.method="bonferroni", verbose=T) {
	all_segments <- CNAobject$output
	
	segments_per_sample <- lapply(unique(all_segments$ID), function(sample) all_segments[all_segments$ID == sample,-1])
	names(segments_per_sample) <- unique(all_segments$ID)
	
	# Data check
	if(length(segments_per_sample) != ncol(sample_intensity))
		stop("The number of samples in CNAobject must match ncols(sample_intensity)")
	
	if(nrow(sample_intensity) != nrow(control_intensity))
		stop("sample_intensity and control_intensity must have the same number of rows")
	
	if(is.null(rownames(site_annotation)))
		rownames(site_annotation) <- site_annotation$ID
	
	# Ensure row consistent ordering
	control_intensity <- control_intensity[rownames(sample_intensity), ]
	site_annotation <- site_annotation[rownames(sample_intensity), ]
	
	x <- lapply(1:length(segments_per_sample), function(i) {
			result <- t(apply(segments_per_sample[[i]], 1, function(cnv) {
				# Extract probes
				probes <- site_annotation$CHR == cnv['chrom'] &
					as.numeric(site_annotation$MAPINFO) >= as.numeric(cnv['loc.start']) &
					as.numeric(site_annotation$MAPINFO) <= as.numeric(cnv['loc.end'])
									
				# Compute segment values
				segment_control_int_sum <- colSums(control_intensity[probes,])
				segment_sample_sum <- sum(sample_intensity[probes,i])
				segment_log_ratio <- log2(segment_sample_sum / mean(segment_control_int_sum))
									
				# Compute p-value (2 sided t-test)
				segment_z_score <- (segment_sample_sum - mean(segment_control_int_sum)) / sd(segment_control_int_sum)
				segment_p_value <-  2 * pnorm(-abs(segment_z_score))
									
				# Extract genes in the segment
				genes <- as.character(site_annotation$UCSC_RefGene_Name[probes])
				genes <- unique(unlist(strsplit(x=genes, ";")))
				genes <- paste(genes, collapse=";")
				l <- as.numeric(cnv['loc.end']) - as.numeric(cnv['loc.start'])

				return(c(cnv, seg.length=l, logratio=segment_log_ratio, pvalue=segment_p_value, genes=genes))
			}))
				
			# pvalue correction for multiple testing
			result <- as.data.frame(result)
			result$adjusted.pvalue <- p.adjust(result$pvalue, method=p.adjust.method)
				
			if(verbose)
				message(paste("Processed", names(segments_per_sample)[i]))
			return(result)
		})
	names(x) <- names(CNVobject)
	return(x)
}