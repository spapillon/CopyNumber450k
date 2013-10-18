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
	if(!type %in% c("default", "quantile"))
		stop("Valid normalization type are {default, quantile}")
	
	method <- match.arg(type)
	if(method == "default") {
		object@intensity_matrix  <- normalizeFunNorm450kCN(cnMatrix = intensityMatrix(object)[usedProbes(object), ], 
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
	case_intensity <- cnMatrix[, sample_groups != "control"]
	
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

	object@segments <- formatSegments(segment.smoothed.CNA.object, case_intensity, control_intensity, 
			probesAnnotation(object) , verbose=verbose)
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