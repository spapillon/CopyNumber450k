################################################################################ 

setGeneric("predictSampleSexes", function(object, threshold = -3) standardGeneric("predictSampleSexes"))

# Returns list containing the sex of each sample via sex chromosome methylation
# intensity-based prediction.
setMethod("predictSampleSexes", signature("MethylCNVDataSet"), function(object, threshold) {
    cnQuantiles <- getSummary(object)$cnQuantiles
    diffs <- log2(cnQuantiles$Y[250, ]) - log2(cnQuantiles$X[250, ])
    predicted_sexes <- ifelse(diffs <= threshold, "Female", "Male")
    predicted_sexes
})

################################################################################ 

setGeneric("dropSNPprobes", function(object, maf_threshold = 0) standardGeneric("dropSNPprobes"))

setMethod("dropSNPprobes", signature("MethylCNVDataSet"), function(object, maf_threshold) {
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

setGeneric("normalize", function(object, type = c("functional", "quantile")) {
    standardGeneric("normalize")
})

# Returns a new MethylCNVDataSet object whose intensities have been normalized to
# the user-specified method.
setMethod("normalize", signature("MethylCNVDataSet"), function(object, type) {
    method <- match.arg(type)
    
    if (!method %in% c("functional", "quantile")) {
        stop("Argument [normalization] type must be {default, quantile}.")
    }
    
    sexes <- pData(object)$Sample_Sex
    intensities <- assayData(object)$intensity
    annotation <- annotation(object)
    manifest <- getManifest(object)
    
    if (method == "functional") {
        new_intensities <- functionalNormalization(cnMatrix = intensities, extractedData = getSummary(object), 
            annotation = annotation, manifest = manifest, predictedSex = sexes)
    } else if (method == "quantile") {
        new_intensities <- quantileNormalization(cnMatrix = intensities, annotation = annotation, 
            manifest = manifest, predictedSex = sexes)
    }
    
    assayData(object) <- assayDataNew(storage.mode = "lockedEnvironment", intensity = new_intensities)
    object
})

################################################################################ 

setGeneric("segmentize", function(object, verbose = TRUE, p.adjust.method = "bonferroni", 
    plotting = FALSE) standardGeneric("segmentize"))

# Returns a new MethylCNVDataSet object.
setMethod("segmentize", signature("MethylCNVDataSet"), function(object, verbose, 
    p.adjust.method, plotting) {
    
    intensities <- assayData(object)$intensity
    
    groups <- pData(object)$Sample_Group
    sexes <- pData(object)$Sample_Sex
    sampleNames <- sampleNames(object)[groups != "control"]
    
    annotation <- fData(object)
    
    control_intensity <- intensities[, groups == "control"]
    case_intensity <- as.matrix(intensities[, groups != "control"])
    
    control_sexes <- sexes[groups == "control"]
    case_sexes <- sexes[groups != "control"]
    
    control_medians <- apply(control_intensity, 1, median, na.rm = TRUE)
    cases_log2 <- log2(case_intensity/control_medians)
    
    if (is.vector(cases_log2)) {
        cases_log2 <- as.matrix(cases_log2)
    }
    
    # TODO: Send this to the DESCRIPTION file ?
    require(DNAcopy)
    # ---
    
    CNA.object <- CNA(cases_log2, ordered(annotation$chr, levels = c(paste("chr", 
        1:22, sep = ""), "chrX", "chrY")), as.numeric(annotation$pos), data.type = "logratio", 
        sampleid = sampleNames)
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
            pdf(paste(names(segments_per_sample)[i], ".pdf", sep = ""))
        }
        
        result <- t(apply(segments_per_sample[[i]], 1, function(cnv) {
            # When assessing sexual chromosomes, only consider controls of the same sex as
            # the sample.
            if (cnv["chrom"] %in% c("chrX", "chrY")) {
                used_controls <- control_sexes == case_sexes[i]
            } else {
                used_controls <- rep(TRUE, ncol(control_intensity))
            }
            
            # Extract probes
            probes <- annotation$chr == cnv["chrom"] & as.numeric(annotation$pos) >= 
                as.numeric(cnv["loc.start"]) & as.numeric(annotation$pos) <= as.numeric(cnv["loc.end"])
            
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
            genes <- as.character(annotation$UCSC_RefGene_Name[probes])
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

################################################################################  
