# TODO: Add comment
# 
# Author: spapillo
###############################################################################


setClass("CNVObject", representation(RGSetSummary = "list", segments = "list", filters = "list", sample_groups = "character", sample_sexes = "character",
				sample_names = "character" ,used_probes = "logical", probe_annotation = "data.frame", intensity_matrix = "matrix", is_normalized = "logical"),
				contains=c("RGChannelSet"))

CNVObject <- function(RGset) {
	return(new("CNVObject", RGset))
}

setMethod("initialize", signature("CNVObject"), function(.Object, RGset) {
	
	message("Extracting RGset data")
	.Object@RGSetSummary <- extractFromRGSet450k(RGset)
	mset <- preprocessRaw(RGset)

	message("Building intensity matrix")
	intensity_matrix <- getMeth(mset) + getUnmeth(mset)
	.Object@intensity_matrix <- intensity_matrix
	
	message("Loading probe annotation")
	# This here, sucks
	require(GEOquery)
	annotation <- Table(getGEO('GPL13534'))
	# until here
	rownames(annotation) <- as.character(annotation$ID)
	annotation <- annotation[as.character(rownames(intensityMatrix(.Object))),]
	.Object@probe_annotation <- annotation
	.Object@used_probes <- rep(TRUE, nrow(annotation))
	sampleSexes(.Object) <- predictSex(.Object)
	.Object@is_normalized <- FALSE

	.Object
})

setMethod("filterSNPProbes", signature("CNVObject"), function(object) {
	annotation <- probesAnnotation(object)
	usedProbes(object) <- usedProbes(object) & 
						annotation$Probe_SNPs == "" & 
						annotation$Probe_SNPs_10 == ""
	object
})

setMethod("filterVariantProbes", signature("CNVObject"), function(object, variance_centile) {
	if(!is.numeric(variance_centile) || length(variance_centile) != 1 || variance_centile <= 0 | variance_centile >= 1)
		stop("The variance_centile argument is the % of the most variant sites you want to remove. It should be ]0,1[")
	controls <- sampleGroups(object) == "control"
	probe_variance <- apply(intensityMatrix(object)[ , controls], 1, var, na.rm=T)
	threshold <- quantile(probe_variance, probs=variance_centile)
	usedProbes(object) <- usedProbes(object) & probe_variance < threshold
	object		
})


setMethod("predictSex", signature("CNVObject"), function(object, threshold) {
	cnQuantiles <- object@RGSetSummary$cnQuantiles
	diffs <- log2(cnQuantiles$Y[250, ]) - log2(cnQuantiles$X[250, ])
	predicted_sexes <- ifelse(diffs <= threshold, "Female", "Male")
	return(predicted_sexes)
})

setMethod("normalize", signature("CNVObject"), function(object, type) {
	if(isNormalized(object))
		stop("This object has already been normalized.")
	
	method <- match.arg(type)
	if(!method %in% c("functional", "quantile"))
		stop("Valid normalization type are {default, quantile}")

	if(method == "functional") {
		object@intensity_matrix  <- functionalNormalization(cnMatrix = intensityMatrix(object)[usedProbes(object), ], 
				extractedData = RGSetSummary(object), predictedSex = sampleSexes(object))
	} else if(method == "quantile") {
		object@intensity_matrix  <- quantileNormalization(cnMatrix = intensityMatrix(object)[usedProbes(object), ], 
				predictedSex = sampleSexes(object))
	}
	
	# The returned matrix has dropped the FALSE usedProbes, remove other objects accordingly
	probesAnnotation(object) <- probesAnnotation(object)[usedProbes(object), ]
	usedProbes(object) <- rep(TRUE, sum(usedProbes(object)))
	object@is_normalized <- TRUE
	object
})

setMethod("buildSegments", signature("CNVObject"), function(object, verbose) {
	if(length(segments(object)) > 0)
		stop("This object has already been segmented")
	used_probes <- usedProbes(object)
	cnMatrix <- intensityMatrix(object)[used_probes, ]
	annotation <- probesAnnotation(object)[used_probes, ]
	sample_groups <- sampleGroups(object)
	sample_names <- sampleNames(object)[sample_groups != "control"]
	
	control_intensity <- cnMatrix[, sample_groups == "control"]
	case_intensity <- as.matrix(cnMatrix[, sample_groups != "control"])
  
  control_sexes <- sampleSexes(object)[sample_groups == "control"]
  case_sexes <- sampleSexes(object)[sample_groups != "control"]
	
	control_medians <- apply(control_intensity, 1, median, na.rm=T)
	cases_log2 <- log2(case_intensity / control_medians)

	if(is.vector(cases_log2))
		cases_log2 <- as.matrix(cases_log2)

	require(DNAcopy)
	CNA.object <- CNA(cases_log2, ordered(annotation$CHR), as.numeric(annotation$MAPINFO),
			data.type="logratio", sampleid=sample_names)
	smoothed.CNA.object <- smooth.CNA(CNA.object)
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width=5 ,verbose=1, nperm=10000, 
			alpha=0.01, undo.splits="sdundo", undo.SD=2)

  ### formatSegment
  all_segments <- segment.smoothed.CNA.object$output
  
  segments_per_sample <- lapply(unique(all_segments$ID), function(sample) all_segments[all_segments$ID == sample,-1])
  names(segments_per_sample) <- unique(all_segments$ID)
  plotting <- TRUE
  p.adjust.method <- "bonferroni"
  x <- lapply(1:length(segments_per_sample), function(i) {
        if(plotting)
          pdf(paste(path, "/", names(segments_per_sample)[i], ".pdf", sep=""))
        result <- t(apply(segments_per_sample[[i]], 1, function(cnv) {
                  # When assessing sexual chromosomes, consider controls
                  # of the same sex as the sample.
                  if(cnv['chrom'] %in% c("X", "Y"))
                    used_controls <- control_sexes == case_sexes[i]
                  else
                    used_controls <- rep(TRUE, ncol(control_intensity))
                  # Extract probes
                  probes <- probesAnnotation(object)$CHR == cnv['chrom'] &
                      as.numeric(probesAnnotation(object)$MAPINFO) >= as.numeric(cnv['loc.start']) &
                      as.numeric(probesAnnotation(object)$MAPINFO) <= as.numeric(cnv['loc.end'])
                  
                  # Compute segment values
                  control_int_sum <- colSums(control_intensity[probes, used_controls])
                  control_mean <- mean(control_int_sum)
                  control_sd <- sd(control_int_sum)
                  sample_sum <- sum(case_intensity[probes,i])
                  
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
                  genes <- as.character(probesAnnotation(object)$UCSC_RefGene_Name[probes])
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
  object@segments <- x
	object
})

setMethod("createFilters", signature("CNVObject"), function(object, tick.threshold, pvalue.threshold, breakpoints) {
	if(length(filters(object)) > 0)
		warning("Previous filters existed, overwriting with new ones")
	filters <- lapply(segments(object), function(sample) as.numeric(sample[,'num.mark']) >= tick.threshold &
														as.numeric(sample[,'adjusted.pvalue']) <= pvalue.threshold)
	object@filters <- filters
	object
})

setMethod("show", signature("CNVObject"), function(object) {
	str(object)
})