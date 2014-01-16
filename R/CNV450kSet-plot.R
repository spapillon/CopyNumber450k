################################################################################ 

setMethod("plotSample", signature("CNV450kSet"), function(object, index, chr, start, 
    end, ...) {
    if (length(getSegments(object)) == 0) {
        stop("Object has not been segmentized yet.")
    }
    
    sample_segments <- getSegments(object)[[index]]
    significant_segments <- sample_segments$isSignificant
    sample_name <- names(getSegments(object))[index]
    
    if (missing(chr)) {
        # Plotting the whole genome
        chromosomes <- unique(sample_segments[, "chrom"])
        site_per_chr <- cumsum(c(0, sapply(chromosomes, function(chr) max(as.numeric(sample_segments[sample_segments[, 
            "chrom"] == chr, "loc.end"])))))
        offset <- site_per_chr - min(as.numeric(sample_segments[sample_segments[, 
            "chrom"] == "chr1", "loc.start"]))
        start <- 0
        end <- max(site_per_chr)
        x_axis_type <- "n"
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
    
    yMin <- min(c(-1, as.numeric(sample_segments[significant_segments, "seg.mean"])))
    yMax <- max(c(1, as.numeric(sample_segments[significant_segments, "seg.mean"])))
    myPlot <- plot(range(start, end), range(yMin, yMax), type = "n", xaxt = x_axis_type, 
            ylab="L-value", xlab="", ...)
    
    if (missing(chr)) {
        xlabs <- sapply(2:length(site_per_chr), function(j) {
            ((site_per_chr[j] - site_per_chr[j - 1])/2) + site_per_chr[j - 1]
        })
        axis(1, at = xlabs, labels = chromosomes, lty = 0, las=2, ...)
        abline(v = site_per_chr, lty = 3)
    }
    
    lapply(1:length(chromosomes), function(i) {
        used_segments <- sample_segments[, "chrom"] == chromosomes[i]
        colors <- ifelse(significant_segments[used_segments], "red", "black")
        starts <- as.numeric(sample_segments[used_segments, "loc.start"]) + offset[i]
        ends <- as.numeric(sample_segments[used_segments, "loc.end"]) + offset[i]
        y <- as.numeric(sample_segments[used_segments, "seg.mean"])
        graphics::segments(starts, y, ends, y, col = colors, lwd = 2, lty = 1)
    })
    
    myPlot
})


################################################################################ 

setMethod("getColoring", signature("CNV450kSet"), function(object, color.by, color.function) {
    coloring <- match.arg(color.by)
    
    if (coloring == "array.row") {
        if (!"Array" %in% colnames(pData(object))) {
            stop("Column \"Array\" must be present in pData(object) in order for coloring = \"array.row\" to be used.")
        }
        
        samples <- substr(pData(object)$Array, 1, 3)
    } else if (coloring == "array.col") {
        if (!"Array" %in% colnames(pData(object))) {
            stop("Column \"Array\" must be present in pData(object) in order for coloring = \"array.col\" to be used.")
        }
        
        samples <- substr(pData(object)$Array, 4, 6)
    } else if (coloring == "sample.group") {
        if (!"Sample_Group" %in% colnames(pData(object))) {
            stop("Column \"Sample_Group\" must be present in pData(object) in order for coloring = \"sample.group\" to be used.")
        }
        
        samples <- pData(object)$Sample_Group
    } else if (coloring == "slide") {
        if (!"Slide" %in% colnames(pData(object))) {
            stop("Column \"Slide\" must be present in pData(object) in order for coloring = \"slide\" to be used.")
        }
        
        samples <- pData(object)$Slide
    } else if (coloring == "origin") {
        if (!"Origin" %in% colnames(pData(object))) {
            stop("Column \"Origin\" must be present in pData(object) in order for coloring = \"origin\" to be used.")
        }
        
        samples <- pData(object)$Origin
    } else {
        stop("Argument color.by must be {\"array.row\", \"array.col\", \"sample.group\", \"slide\", \"origin\"}.")
    }
    
    cols <- color.function(length(unique(samples)))
    col_vec <- cols[match(samples, unique(samples))]
    
    list(sample.colors = col_vec, groups = unique(samples), group.colors = cols)
})

################################################################################ 

setMethod("plotDensity", signature("CNV450kSet"), function(object, color.by, color.function, 
    legend.position, ...) {
    
    intensities <- assayData(object)$intensity
    
    coloring <- getColoring(object, color.by, color.function)
    myPlot <- plot(density(intensities), col = coloring$sample.colors[1], ylim = c(0, 
        9e-05), ...)
    sapply(2:ncol(intensities), function(i) lines(density(intensities[, i]), col = coloring$sample.colors[i]))
    
    if (!is.null(legend.position)) {
        legend(legend.position, legend = coloring$groups, fill = coloring$group.colors)
    }
    
    myPlot
})

################################################################################ 

setMethod("plotPCA", signature("CNV450kSet"), function(object, color.by, color.function, 
    legend.position, ...) {
    
    intensities <- assayData(object)$intensity
    
    coloring <- getColoring(object, color.by, color.function)
    pca <- prcomp(t(intensities))
    myPlot <- plot(pca$x, col = coloring$sample.colors, ...)
    
    if (!is.null(legend.position)) {
        legend(legend.position, legend = coloring$groups, fill = coloring$group.colors)
    }
    
    myPlot
})

################################################################################  

setMethod("write.csv", signature("CNV450kSet"), function(object, ...) {
    segment_list <-getSegments(object)
    sample_names <- names(segment_list)
    output <- Reduce(rbind, lapply(1:length(segment_list), function(i) {
                            name <- sample_names[i]
                            datum <- segment_list[[i]]
                            x <- cbind(Sample=rep(name, nrow(datum)), datum)
                            return(x)
                        }))
    write.csv(output, ...)
})
