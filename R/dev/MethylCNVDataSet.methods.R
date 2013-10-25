################################################################################ 

setGeneric("predictSampleSexes", function(object, threshold = -3) standardGeneric("predictSampleSexes"))

# Returns list containing the sex of each sample via sex chromosome methylation
# intensity-based prediction.
setMethod("predictSampleSexes", signature("MethylCNVDataSet"), function(object, threshold) {
    cnQuantiles <- object@summary$cnQuantiles
    diffs <- log2(cnQuantiles$Y[250, ]) - log2(cnQuantiles$X[250, ])
    predicted_sexes <- ifelse(diffs <= threshold, "Female", "Male")
    predicted_sexes
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
    
    sexes <- pData(object)$sex
    intensities <- assayData(object)$intensity
    
    if (method == "functional") {
        fData(object)$intensity <- functionalNormalization(cnMatrix = intensities, 
            extractedData = object@summary, predictedSex = sexes)
    } else if (method == "quantile") {
        fData(object)$intensity <- quantileNormalization(cnMatrix = intensities, 
            predictedSex = sexes)
    }
    
    object
})

################################################################################ 

setGeneric("segmentize", function(object, verbose = TRUE, p.adjust.method = "bonferroni", 
    plotting = FALSE) standardGeneric("segmentize"))

# Returns a new MethylCNVDataSet object. TODO: Fix this lol
setMethod("segmentize", signature("MethylCNVDataSet"), function(object, verbose, 
    p.adjust.method, plotting) {
    
    intensities <- assayData(object)$intensity
    
    groups <- pData(object)$groups
    sexes <- pData(object)$sex
    sampleNames <- sampleNames(object)[groups != "control"]
    
    # TODO: Annotation?
    annotation <- fData(object)
    
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
    
    CNA.object <- CNA(cases_log2, ordered(annotation$chr, levels=c(paste("chr", 1:22, sep=""), "chrX", "chrY"))
            , as.numeric(annotation$pos), data.type = "logratio", sampleid = sampleNames)
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
                as.numeric(cnv["loc.start"]) & as.numeric(annotation$pos) <= 
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
            #genes <- as.character(allAnnotations$UCSC_RefGene_Name[probes])
            #genes <- unique(unlist(strsplit(x = genes, ";")))
            #genes <- paste(genes, collapse = ";")
            l <- as.numeric(cnv["loc.end"]) - as.numeric(cnv["loc.start"])
            #c(cnv, seg.length = l, pvalue = p_value, genes = genes, ctrl.mean = control_mean, 
            #    ctrl.sd = control_sd, sample.value = sample_sum, z = z_score)
            c(cnv, seg.length = l, pvalue = p_value, ctrl.mean = control_mean, 
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

setGeneric("plotSample", function(object, index, chr, start, end) standardGeneric("plotSample"))

setMethod("plotSample", signature("MethylCNVDataSet"), function(object, index, chr, start, 
            end) {
          sample_segments <- object@segments[[index]]
          sample_filters <- rep(FALSE, nrow(sample_segments))
          sample_name <- sampleNames(object)[index]
          
          if (missing(chr)) {
        # Plotting the whole genome
        chromosomes <- unique(sample_segments[, "chrom"])
        
        
        site_per_chr <- cumsum(c(0, sapply(chromosomes, function(chr) max(as.numeric(sample_segments[sample_segments[, 
                                                  "chrom"] == chr, "loc.end"])))))
        offset <- site_per_chr - min(as.numeric(sample_segments[sample_segments[, 
                                    "chrom"] == 'chr1', "loc.start"]))
        start <- 0
        end <- max(site_per_chr)
        x_axis_type <- 'n'
          } else {
        # Plotting a region
        if (missing(start)) {
          start <- 0
        }
        
        if (missing(end)) {
          end <- max(sample_segments[sample_segments[, "chrom"] == chr, "loc.end"])
        }
        
        chromosomes <- chr
        offset <- 0
        x_axis_type <- NULL
          }
          
          yMin <- min(c(-1, as.numeric(sample_segments[sample_filters, "seg.mean"])))
          yMax <- max(c(1, as.numeric(sample_segments[sample_filters, "seg.mean"])))
          myPlot <- plot(range(start, end), range(yMin, yMax), type = "n", xaxt = x_axis_type, 
                  xlab = "", ylab = "", main = sample_name)
          
          if (missing(chr)) {
        xlabs <- sapply(2:length(site_per_chr), function(j) {
                          ((site_per_chr[j] - site_per_chr[j - 1])/2) + site_per_chr[j - 1]
                    })
        
        axis(1, at = xlabs, labels = chromosomes, lty = 0)
        abline(v = site_per_chr, lty = 3)
          }
          
          lapply(1:length(chromosomes), function(i) {
                    used_segments <- sample_segments[, "chrom"] == chromosomes[i]
                    colors <- ifelse(sample_filters[used_segments], "red", "black")
                    starts <- as.numeric(sample_segments[used_segments, "loc.start"]) + offset[i]
                    ends <- as.numeric(sample_segments[used_segments, "loc.end"]) + offset[i]
                    y <- as.numeric(sample_segments[used_segments, "seg.mean"])
                    graphics::segments(starts, y, ends, y, col = colors, lwd = 2, lty = 1)
              })
          
          myPlot
    })