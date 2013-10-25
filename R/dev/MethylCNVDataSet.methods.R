setGeneric("predictSampleSexes", function(object, threshold = -3) standardGeneric("predictSampleSexes"))

# Returns list containing the sex of each sample via sex chromosome methylation
# intensity-based prediction.
setMethod("predictSampleSexes", signature("MethylCNVDataSet"), function(object, threshold) {
    cnQuantiles <- object@summary$cnQuantiles
    diffs <- log2(cnQuantiles$Y[250, ]) - log2(cnQuantiles$X[250, ])
    predicted_sexes <- ifelse(diffs <= threshold, "Female", "Male")
})



setGeneric("normalize", function(object, type = c("functional", "quantile")) {
    standardGeneric("normalize")
})

# Returns a new MethylCNVDataSet object whose intensities have been normalized to
# the user-specified method.
setMethod("normalize", signature("MethylCNVDataSet"), function(object, type) {
    method <- match.arg(type)
    
    if (!method %in% c("functional", "quantile")) {
        stop("Expected argument [normalization] type to be {default, quantile}.")
    }
    
    sexes <- pData(object)$sex
    featuresUsed <- fData(object)$isUsed
    intensities <- assayData(object)$intensity[featuresUsed, ]
    
    if (method == "functional") {
        fData(object)$intensity <- functionalNormalization(cnMatrix = intensities, 
            extractedData = object@summary, predictedSex = sexes)
    } else if (method == "quantile") {
        fData(object)$intensity <- quantileNormalization(cnMatrix = intensities, 
            predictedSex = sexes)
    }
    
    # The returned matrix has dropped the FALSE usedProbes, remove other objects
    # accordingly.
    
    # TODO: Annotation? What happens to these? lol
    
    # probesAnnotation(object) <- probesAnnotation(object)[usedProbes(object), ]
    
    # usedProbes(object) <- rep(TRUE, sum(usedProbes(object)))
    
    # object@is_normalized <- TRUE
    
    object
})



setGeneric("segmentize", function(object, verbose = TRUE, p.adjust.method = "bonferroni", 
    plotting = FALSE) standardGeneric("segmentize"))

# Returns a new MethylCNVDataSet object. TODO: Fix this lol
setMethod("segmentize", signature("MethylCNVDataSet"), function(object, verbose, 
    p.adjust.method, plotting) {
    if (length(segments(object)) > 0) {
        stop("Object has already been segmented.")
    }
    
    featuresUsed <- fData(object)$isUsed
    intensities <- assayData(object)$intensity[featuresUsed, ]
    
    groups <- pData(object)$groups
    sexes <- pData(object)$sex
    sampleNames <- sampleNames(object)[groups != "control"]
    
    # TODO: Annotation?
    filteredAnnotations <- probesAnnotation(object)[featuresUsed, ]
    allAnnotations <- probesAnnotation(object)
    
    control_intensity <- intensities[, groups == "control"]
    case_intensity <- as.matrix(intensities[, groups != "control"])
    
    control_sexes <- sexes[groups == "control"]
    case_sexes <- sexes[groups != "control"]
    
    control_medians <- apply(control_intensity, 1, median, na.rm = T)
    cases_log2 <- log2(case_intensity/control_medians)
    
    if (is.vector(cases_log2)) {
        cases_log2 <- as.matrix(cases_log2)
    }
    
    # TODO: Send this to the DESCRIPTION file ?
    require(DNAcopy)
    # ---
    
    CNA.object <- CNA(cases_log2, ordered(filteredAnnotations$CHR), as.numeric(filteredAnnotations$MAPINFO), 
        data.type = "logratio", sampleid = sampleNames)
    smoothed.CNA.object <- smooth.CNA(CNA.object)
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width = 5, verbose = 1, 
        nperm = 10000, alpha = 0.01, undo.splits = "sdundo", undo.SD = 2)
    
    # formatSegment
    all_segments <- segment.smoothed.CNA.object$output
    
    segments_per_sample <- lapply(unique(all_segments$ID), function(sample) all_segments[all_segments$ID == 
        sample, -1])
    names(segments_per_sample) <- unique(all_segments$ID)
    
    x <- lapply(1:length(segments_per_sample), function(i) {
        if (plotting) {
            pdf(paste(path, "/", names(segments_per_sample)[i], ".pdf", sep = ""))
        }
        
        result <- t(apply(segments_per_sample[[i]], 1, function(cnv) {
            # When assessing sexual chromosomes, only consider controls of the same sex as
            # the sample.
            if (cnv["chrom"] %in% c("X", "Y")) {
                used_controls <- control_sexes == case_sexes[i]
            } else {
                used_controls <- rep(TRUE, ncol(control_intensity))
            }
            
            # Extract probes
            probes <- allAnnotations$CHR == cnv["chrom"] & as.numeric(allAnnotations$MAPINFO) >= 
                as.numeric(cnv["loc.start"]) & as.numeric(allAnnotations$MAPINFO) <= 
                as.numeric(cnv["loc.end"])
            
            # Compute segment values
            control_int_sum <- colSums(control_intensity[probes, used_controls])
            control_mean <- mean(control_int_sum)
            control_sd <- sd(control_int_sum)
            sample_sum <- sum(case_intensity[probes, i])
            
            # Compute p-value (2 sided t-test)
            z_score <- (sample_sum - control_mean)/control_sd
            p_value <- 2 * pnorm(-abs(z_score))
            if (plotting) {
                window <- range(c(control_int_sum, sample_sum))
                title <- paste("chr ", cnv["chrom"], " :  ", cnv["loc.start"], "-", 
                  cnv["loc.end"], sep = "")
                hist(control_int_sum, main = title, xlim = window, breaks = 20)
                abline(v = sample_sum, lty = 3, lwd = 2, col = "red")
            }
            
            # Extract genes in the segment
            genes <- as.character(allAnnotations$UCSC_RefGene_Name[probes])
            genes <- unique(unlist(strsplit(x = genes, ";")))
            genes <- paste(genes, collapse = ";")
            l <- as.numeric(cnv["loc.end"]) - as.numeric(cnv["loc.start"])
            c(cnv, seg.length = l, pvalue = p_value, genes = genes, ctrl.mean = control_mean, 
                ctrl.sd = control_sd, sample.value = sample_sum, z = z_score)
        }))
        
        if (plotting) {
            dev.off()
        }
        
        # pvalue correction for multiple testing
        result <- as.data.frame(result, stringsAsFactors = F)
        result$adjusted.pvalue <- p.adjust(result$pvalue, method = p.adjust.method)
        
        if (verbose) {
            message(paste("Processed", names(segments_per_sample)[i]))
        }
        
        result
    })
    
    names(x) <- names(segments_per_sample)
    object@segments <- x
    object
}) 
