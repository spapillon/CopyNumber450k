setReplaceMethod("probesAnnotation", signature("CNVObject"), function(object, value) {
	if(!is.data.frame(value))
		stop("Input parameter needs to be of type data.frame")
	if(length(setdiff(c("MAPINFO", "CHR", "Probe_SNPs", "Probe_SNPs_10"), colnames(value))) > 0)
		stop("Colnames of input parameter needs to contain {\"MAPINFO\", \"CHR\", \"Probe_SNPs\", \"Probe_SNPs_10\"}")
	object@probe_annotation <- value
	object
})

setReplaceMethod("sampleGroups", signature("CNVObject"), function(object, value) {
	if(!is.character(value) || length(value) != ncol(intensityMatrix(object)) || length(setdiff(c("control"), value)) > 0)
		stop("Input parameter needs to be a character vector of length ncol(intensity_matrix) with at least one \"control\" entry")
	object@sample_groups <- value
	object
})

setReplaceMethod("sampleSexes", signature("CNVObject"), function(object, value) {
	if(!is.numeric(value) || length(value) != ncol(intensityMatrix(object)) || length(setdiff(c(1,2), value)) > 0)
		stop("Input parameter needs to be a numerci vector of length ncol(intensity_matrix) containing values {1,2} (male, female)")
	predicted_values <- predictSex(object@RGSetSummary, -3)
	if(sum(predicted_values != value))
	{
		warning(paste("According to the X and Y chromosome intensities, it seems you have entered the wrong sex for samples", 
						paste(sampleNames[predicted_values != value], sep=" ,"), sep=""))
	}
	object@sample_sexes <- value
	object
})

setReplaceMethod("sampleNames", signature("CNVObject"), function(object, value) {
	if(!is.character(value) || length(value) != ncol(intensityMatrix(object)))
		stop("Input parameter needs to be a character vector of length ncol(intensity_matrix)")
	object@sample_names <- value
	object
})

setReplaceMethod("usedProbes", signature("CNVObject"), function(object, value) {
	if(!is.logical(value) || length(value) != nrow(intensityMatrix(object)))
		stop("Input parameter needs to be a logical vector of length nrow(intensity_matrix)")
	object@used_probes <- value
	object
})
