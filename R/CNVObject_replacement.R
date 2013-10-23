# Replacement methods

setReplaceMethod("probesAnnotation", signature("CNVObject"), function(object, value) {
    if (!is.data.frame(value)) {
        stop("Expected parameter value to be a data.frame.")
    }
    
    if (length(setdiff(c("MAPINFO", "CHR", "Probe_SNPs", "Probe_SNPs_10"), colnames(value))) > 
        0) {
        stop("Expected parameter value's colnames to contain {\"MAPINFO\", \"CHR\", \"Probe_SNPs\", \"Probe_SNPs_10\"}.")
    }
    
    object@probe_annotation <- value
    object
})

setReplaceMethod("sampleGroups", signature("CNVObject"), function(object, value) {
    if (!is.character(value)) {
        stop("Expected parameter value to be a character vector.")
    } else if (length(value) != ncol(intensityMatrix(object))) {
        stop("Expected parameter value to be of length equal to the number of samples.")
    } else if (!any(value == "control")) {
        # TODO: Use more informative desc
        stop("Expected parameter value to contain at least one \"control\" entry.")
        # ---
    }
    
    object@sample_groups <- value
    object
})

setReplaceMethod("sampleSexes", signature("CNVObject"), function(object, value) {
    if (!is.character(value)) {
        stop("Expected parameter value to be a character vector.")
    } else if (length(value) != ncol(intensityMatrix(object))) {
        stop("Expected parameter value to be of length equal to the number of samples.")
    } else if (length(setdiff(value, c("Male", "Female"))) > 0) {
        # TODO: More descriptive
        stop("Expected parameter value to contain values {\"Male\",\"Female\"}.")
        # ---
    }
    
    predicted_values <- predictSex(object)
    if (sum(predicted_values != value) > 0) {
        warning(paste("According to the X and Y chromosome intensities, it seems you have entered the wrong sex for samples", 
            paste(sampleNames(object)[predicted_values != value], collapse = ", "), 
            sep = ": "))
    }
    
    object@sample_sexes <- value
    object
})

setReplaceMethod("sampleNames", signature("CNVObject"), function(object, value) {
    if (!is.character(value)) {
        stop("Expected parameter value to be a character vector.")
    } else if (length(value) != ncol(intensityMatrix(object))) {
        stop("Expected parameter value to be of length equal to the number of samples.")
    }
    
    object@sample_names <- value
    object
})

setReplaceMethod("usedProbes", signature("CNVObject"), function(object, value) {
    if (!is.logical(value)) {
        stop("Expected parameter value to be a logical vector.")
    } else if (length(value) != nrow(intensityMatrix(object))) {
        stop("Expected parameter value to be of length equal to the number of probes.")
    }
    
    object@used_probes <- value
    object
})

setReplaceMethod("sampleChipRows", signature("CNVObject"), function(object, value) {
    object@sample_chip_rows <- as.character(value)
    object
})

setReplaceMethod("sampleChipColumns", signature("CNVObject"), function(object, value) {
    object@sample_chip_columns <- as.character(value)
    object
})

setReplaceMethod("sampleChipIDs", signature("CNVObject"), function(object, value) {
    object@sample_chip_ids <- as.character(value)
    object
}) 
