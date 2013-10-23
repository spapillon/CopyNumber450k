# Replacement methods

setReplaceMethod("probesAnnotation", signature("CNVObject"), function(object, annotation) {
    if (!is.data.frame(annotation)) {
        stop("Expected argument [probe] `annotation` to be a data.frame.")
    } else if (length(setdiff(c("MAPINFO", "CHR", "Probe_SNPs", "Probe_SNPs_10"), colnames(annotation))) > 
        0) {
        stop("Expected argument [probe] `annotation` colnames to contain {\"MAPINFO\", \"CHR\", \"Probe_SNPs\", \"Probe_SNPs_10\"}.")
    }
    
    object@probe_annotation <- annotation
    object
})

setReplaceMethod("usedProbes", signature("CNVObject"), function(object, probes) {
    if (!is.logical(probes)) {
        stop("Expected argument `probes` to be a logical vector.")
    } else if (length(probes) != nrow(intensityMatrix(object))) {
        stop("Expected argument `probes` to be of same length as number of probes. Each probe must be given a logical status {\"TRUE\", \"FALSE\"}.")
    }
    
    object@used_probes <- probes
    object
})

setReplaceMethod("sampleGroups", signature("CNVObject"), function(object, groups) {
    if (!is.character(groups)) {
        stop("Expected argument `groups` to be a character vector.")
    } else if (length(groups) != ncol(intensityMatrix(object))) {
        stop("Expected argument `groups` to be of same length as number of samples. Each sample must be assigned to a group.")
    } else if (!any(groups == "control")) {
        stop("Expected argument `groups` to contain at least one \"control\" element. Atleast one sample must be of group \"control\".")
    }
    
    object@sample_groups <- groups
    object
})

setReplaceMethod("sampleSexes", signature("CNVObject"), function(object, sexes) {
    if (!is.character(sexes)) {
        stop("Expected argument `sexes` to be a character vector.")
    } else if (length(sexes) != ncol(intensityMatrix(object))) {
        stop("Expected argument `sexes` to be of same length as number of samples. Each sample must be attributed a sex.")
    } else if (length(setdiff(sexes, c("Male", "Female"))) > 0) {
        stop("Expected argument `sexes` elements to be {\"Male\", \"Female\"}. No other sex can be entered.")
    }
    
    predicted_values <- predictSex(object)
    if (sum(predicted_values != sexes) > 0) {
        warning(paste("According to the X and Y chromosome intensities, it seems you have entered the wrong sex for samples", 
            paste(sampleNames(object)[predicted_values != sexes], collapse = ", "), 
            sep = ": "))
    }
    
    object@sample_sexes <- sexes
    object
})

setReplaceMethod("sampleNames", signature("CNVObject"), function(object, names) {
    if (!is.character(names)) {
        stop("Expected argument `names` to be a character vector.")
    } else if (length(names) != ncol(intensityMatrix(object))) {
        stop("Expected argument `names` to be of same length as number of samples. Each sample must be attributed a name.")
    }
    
    # TODO: Unique names?
    
    object@sample_names <- names
    object
})

setReplaceMethod("sampleChipRows", signature("CNVObject"), function(object, chipRows) {
    object@sample_chip_rows <- as.character(chipRows)
    object
})

setReplaceMethod("sampleChipColumns", signature("CNVObject"), function(object, chipColumns) {
    object@sample_chip_columns <- as.character(chipColumns)
    object
})

setReplaceMethod("sampleChipIDs", signature("CNVObject"), function(object, chipIDs) {
    object@sample_chip_ids <- as.character(chipIDs)
    object
})
 
