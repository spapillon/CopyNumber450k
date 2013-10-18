# TODO: Add comment
# 
# Author: spapillo
###############################################################################


setClass("CNVObject", representation(RGSetSummary = "list", segments = "list", filters = "list", sample_groups = "character", samples_sexes = "numeric",
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
	colnames(intensity_matrix) <- pData(RGset)$Sample_Name
	.Object@intensity_matrix <- intensity_matrix

	sampleGroups(.Object) <- pData(RGset)$Sample_Group
	sampleNames(.Object) <- pData(RGset)$Sample_Name
	
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

setMethod("normalize", signature("CNVObject"), function(object, sex_cutoff) {
	if(isNormalized(object))
		stop("This object has already been normalized.")

	object@intensity_matrix  <- normalizeFunNorm450kCN(cnMatrix = intensityMatrix(object)[usedProbes(object), ], 
			extractedData = RGSetSummary(object), predictedSex = sampleSexes(object))
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


setMethod("findCNV", signature("CNVObject"), function(object, CNVs, type) {
	if(!is.character(CNVs))
		stop("The CNVs arguments needs to be a vector of characters")
	segments_list <- segments(object)
	filters_list <- filters(object)
	x <- sapply(1:length(segments_list), function(i) {
				if(type == "both")
					used_CNVs <- filters_list[[i]]
				else if(type == "gain")
					used_CNVs <- segments_list[[i]][,'logratio'] > 0 & filters_list[[i]]
				else if(type == "loss")
					used_CNVs <- segments_list[[i]][,'logratio'] < 0 & filters_list[[i]]
				
				genes <- unique(unlist(strsplit(x=segments_list[[i]][used_CNVs,'genes'], ";")))
				return(CNVs %in% genes)
			})

	ifelse(is.matrix(x), x <- t(x), x <- as.matrix(x))
	rownames(x) <- names(segments_list)
	colnames(x) <- CNVs
	return(x)
})

setMethod("intersectCNV", signature("CNVObject"), function(object, sample_indices, type) {
	sample_count <- length(segments(object))
	if(missing(sample_indices))
		sample_indices <- 1:sample_count
	if(!is.numeric(sample_indices) || length(setdiff(sample_indices, 1:sample_count)) > 0)
		stop("The sample_indices argument needs to be a integer vector with values corresponding to existing samples indices")

	segments_list <- segments(object)[sample_indices]
	filters_list <- filters(object)[sample_indices]
	x <- unlist(sapply(1:length(segments_list), function(i) {
			if(type == "both")
				used_CNVs <- filters_list[[i]]
			else if(type == "gain")
				used_CNVs <- segments_list[[i]][,'logratio'] > 0 & filters_list[[i]]
			else if(type == "loss")
				used_CNVs <- segments_list[[i]][,'logratio'] < 0 & filters_list[[i]]
			unique(unlist(strsplit(x=segments_list[[i]][used_CNVs,'genes'], ";")))
			}))
	return(sort(table(x), decreasing=T))
})

setMethod("subgroupDifference", signature("CNVObject"), function(object, group1_indices, group2_indices) {
	sample_count <- length(segments(object))
	if(!is.numeric(group1_indices) || length(setdiff(group1_indices, 1:sample_count)) > 0 )
		stop("The group1_indices argument needs to be a integer vector with values corresponding to existing samples indices")		
	if(!is.numeric(group2_indices) || length(setdiff(group1_indices, 1:sample_count)) > 0 )
		stop("The group2_indices argument needs to be a integer vector with values corresponding to existing samples indices")
	if(length(intersect(group1_indices, group2_indices)) > 0)
		warning("Some samples are present in both subgroups")

	group1_size <- length(group1_indices)
	group2_size <- length(group2_indices)

	#Gains
	group1_CNVs <- intersectCNV(object, group1_indices, type="gain")
	group2_CNVs <- intersectCNV(object, group2_indices, type="gain")
	gains <- subgroupDifferenceCNVByType(group1_CNVs, group2_CNVs, group1_size, group2_size) 
	colnames(gains) <- c("Group 1", "Group 2", "pvalue")

	#Losses
	group1_CNVs <- intersectCNV(object, group1_indices, type="loss")
	group2_CNVs <- intersectCNV(object, group2_indices, type="loss")
	losses <- subgroupDifferenceCNVByType(group1_CNVs, group2_CNVs, group1_size, group2_size) 
	colnames(losses) <- c("Group 1", "Group 2", "pvalue")

	return(list(gains=gains, losses=losses))
})
