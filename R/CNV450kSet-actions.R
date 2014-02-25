################################################################################ 

setMethod("dropSNPprobes", signature("CNV450kSet"), function(object, maf_threshold) {
    annotation <- fData(object)
    
    # c('Probe_rs', 'Probe_maf', 'CpG_rs', 'CpG_maf', 'SBE_rs', 'SBE_maf')
    usedProbes <- (is.na(annotation[, "Probe_rs"]) | as.numeric(annotation[, "Probe_maf"]) <= 
        maf_threshold) & (is.na(annotation[, "CpG_rs"]) | as.numeric(annotation[, 
        "CpG_maf"]) <= maf_threshold) & (is.na(annotation[, "SBE_rs"]) | as.numeric(annotation[, 
        "SBE_maf"]) <= maf_threshold)
    
    fData(object) <- fData(object)[usedProbes, ]
    intensities <- assayData(object)$intensity[usedProbes, ]
    assayData(object) <- assayDataNew(storage.mode = "lockedEnvironment", intensity = intensities)
    object
})

################################################################################ 

# Returns a new CNV450kSet object whose intensities have been normalized to the
# user-specified method.
setMethod("normalize", signature("CNV450kSet"), function(object, type = c("functional", 
    "quantile")) {
    method <- match.arg(type)
    
    if (!method %in% c("functional", "quantile")) {
        stop("Argument [normalization] type must be {default, quantile}.")
    }
    
    intensities <- assayData(object)$intensity
    manifest <- getManifest(object)
    
    if (method == "functional") {
        new_intensities <- functionalNormalization(cnMatrix = intensities, extractedData = getSummary(object), 
            manifest = manifest)
    } else if (method == "quantile") {
        new_intensities <- quantileNormalization(cnMatrix = intensities, manifest = manifest)
    }
    
    assayData(object) <- assayDataNew(storage.mode = "lockedEnvironment", intensity = new_intensities)
    object
})

################################################################################ 

# Returns a new CNV450kSet object.
setMethod("segmentize", signature("CNV450kSet"), function(object, verbose, p.adjust.method, 
    min.width, nperm, alpha, undo.splits, undo.SD, trim) {
    if (length(getSegments(object)) != 0 && verbose) {
        warning("Object has already been segmentized.")
    }
    
    intensities <- assayData(object)$intensity
    
    groups <- pData(object)$Sample_Group
    sexes <- pData(object)$Sample_Sex
    sampleNames <- sampleNames(object)[groups != "control"]
    
    annotation <- fData(object)
    
    control_intensity <- intensities[, groups == "control"]
    case_intensity <- as.matrix(intensities[, groups != "control"])
    
    consider_sexes <- !any(is.na(sexes))
    if (consider_sexes) {
        control_sexes <- sexes[groups == "control"]
        case_sexes <- sexes[groups != "control"]
    }
    
    control_medians <- apply(control_intensity, 1, median, na.rm = TRUE)
    cases_log2 <- log2(case_intensity/control_medians)
    
    if (is.vector(cases_log2)) {
        cases_log2 <- as.matrix(cases_log2)
    }
    
    CNA.object <- CNA(cases_log2, ordered(annotation$chr, levels = c(paste("chr", 
        1:22, sep = ""), "chrX", "chrY")), as.numeric(annotation$pos), data.type = "logratio", 
        sampleid = sampleNames)
    smoothed.CNA.object <- smooth.CNA(CNA.object, trim = trim)
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width = min.width, 
        verbose = verbose, nperm = nperm, alpha = alpha, undo.splits = undo.splits, 
        undo.SD = undo.SD, trim = trim)
    
    # formatSegment
    all_segments <- segment.smoothed.CNA.object$output
    
    segments_per_sample <- lapply(unique(all_segments$ID), function(sample) all_segments[all_segments$ID == 
        sample, -1])
    names(segments_per_sample) <- unique(all_segments$ID)
    
    x <- lapply(1:length(segments_per_sample), function(i) {
        
        result <- t(apply(segments_per_sample[[i]], 1, function(cnv) {
            # When assessing sexual chromosomes, only consider controls of the same sex as
            # the sample.
            
            if (cnv["chrom"] %in% c("chrX", "chrY") & consider_sexes) {
                used_controls <- control_sexes == case_sexes[i]
            } else {
                used_controls <- rep(TRUE, ncol(control_intensity))
            }
            
            # Extract probes
            probes <- annotation$chr == cnv["chrom"] & as.numeric(annotation$pos) >= 
                as.numeric(cnv["loc.start"]) & as.numeric(annotation$pos) <= as.numeric(cnv["loc.end"])
            
            # Compute segment values
            if (sum(probes) < 2 || sum(used_controls) < 2) {
                control_int_sum <- sum(control_intensity[probes, used_controls])
            } else {
                control_int_sum <- colSums(control_intensity[probes, used_controls])
            }
            control_mean <- mean(control_int_sum)
            control_sd <- sd(control_int_sum)
            sample_sum <- sum(case_intensity[probes, i])
            
            # Compute p-value (2 sided t-test)
            z_score <- (sample_sum - control_mean)/control_sd
            p_value <- 2 * pnorm(-abs(z_score))
            
            # Extract genes in the segment
            genes <- as.character(annotation$UCSC_RefGene_Name[probes])
            genes <- unique(unlist(strsplit(x = genes, ";")))
            genes <- paste(genes, collapse = ";")
            l <- as.numeric(cnv["loc.end"]) - as.numeric(cnv["loc.start"])
            c(cnv, seg.length = l, pvalue = p_value, genes = genes, ctrl.mean = control_mean, 
                ctrl.sd = control_sd, sample.value = sample_sum, z = z_score)
        }))
        
        # pvalue correction for multiple testing
        result <- as.data.frame(result, stringsAsFactors = FALSE)
        result$adjusted.pvalue <- p.adjust(result$pvalue, method = p.adjust.method)
        if (verbose > 0) {
            message(paste("Processed", names(segments_per_sample)[i]))
        }
        
        result
    })
    
    names(x) <- names(segments_per_sample)
    object@segments <- x
    computeSignificance(object)
})


################################################################################ 

# Returns a new CNV450kSet object.
setMethod("computeSignificance", signature("CNV450kSet"), function(object, p.value.threshold, 
    num.mark.threshold) {
    current_segments <- getSegments(object)
    new_segments <- lapply(current_segments, function(sample) {
        significant <- as.numeric(sample$adjusted.pvalue) <= p.value.threshold & 
            as.numeric(sample$num.mark) >= num.mark.threshold
        sample$isSignificant <- significant
        sample
    })
    object@segments <- new_segments
    object
})

################################################################################  
