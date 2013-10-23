#
# Plotting Methods
#

setMethod("plotSample",
          signature("CNVObject"),
          function(object,
                   index,
                   chr,
                   start,
                   end) {
                 
	sample_segments <- segments(object)[[index]]
	sample_filters  <- filters(object)[[index]]
	sample_name     <- sampleNames(object)[index]
	
	if (missing(chr)) { # Plotting the whole genome
    chromosomes <- unique(sample_segments[, 'chrom'])
		idx <- which(chromosomes %in% c("X", "Y"))
		chromosomes <- sort(as.numeric(chromosomes[-idx]))
    # TODO: Refactor this ...
		site_per_chr <- cumsum(c(0, sapply(chromosomes,
                function(chr) max(as.numeric(sample_segments[sample_segments[, 'chrom']== chr,'loc.end'])))))
		offset <- site_per_chr - min(as.numeric(sample_segments[sample_segments[,'chrom'] == 1, 'loc.start']))
		start <- 0
		end <- max(site_per_chr)
	}
  else { # Plotting a region
		if (missing(start)) {
			start <- 0
    }
		
    if (missing(end)) {
			end <- max(sample_segments[sample_segments[,'chrom'] == chr, 'loc.end'])
    }
    
		chromosomes <- chr
		offset <- 0
	}
	
	yMin <- min(c(-1, as.numeric(sample_segments[sample_filters, 'seg.mean'])))
	yMax <- max(c(1, as.numeric(sample_segments[sample_filters, 'seg.mean'])))
	myPlot <- plot(range(start, end),
                 range(yMin, yMax),
                 type = 'n',
                 xaxt = 'n',
                 xlab = "",
                 ylab = "",
                 main = sample_name)
	
	if (missing(chr)) {
    # TODO: Review this formula
		xlabs <- sapply(2:length(site_per_chr), function(j) {
          ((site_per_chr[j] - site_per_chr[j - 1]) / 2) + site_per_chr[j - 1]
    })
		
    axis(1, at = xlabs, labels = chromosomes, lty = 0)
		abline(v = site_per_chr, lty = 3)
	}
  
	lapply(1:length(chromosomes), function(i) { 
		used_segments <- sample_segments[, 'chrom'] == chromosomes[i]
		colors <- ifelse(sample_filters[used_segments], "red", "black")		
		starts <- as.numeric(sample_segments[used_segments, 'loc.start']) + offset[i]
		ends <- as.numeric(sample_segments[used_segments, 'loc.end']) + offset[i]
		y <- as.numeric(sample_segments[used_segments, 'seg.mean'])
		graphics::segments(starts, y, ends, y, col=colors, lwd = 2, lty = 1)
	})

	myPlot
})
			


setMethod("plotSex",
          signature("CNVObject"),
          function(object) {
            
	cnQuantiles <- object@RGSetSummary$cnQuantiles
	plot(log2(cnQuantiles$Y[250, ]) - log2(cnQuantiles$X[250, ]),
       rep(0, length(cnQuantiles$Y[250, ])))
	abline(v = -3, lty = 3, col = "red")
	text(x = -4, y = 0.5, label = "Female", col = "red")
	text(x = -1, y = 0.5, label = "Male", col = "red")
})

setMethod("plotDensity",
    signature("CNVObject"),
    function(object, color.by, color.function) {
  
	col_vec <- getColors(object, color.by, color.function)
  plot(density(int_matrix), col=col_vec[1], ylim=c(0,9e-5))
	sapply(2:ncol(int_matrix), function(i) lines(density(int_matrix[,i]), col=col_vec[i]))	
})

setMethod("getColors",
    signature("CNVObject"),
    function(object, color.by, color.function) {

  coloring <- match.arg(color.by)
  if(coloring == "sentrix.row") {
    samples <- sampleChipRows(object)
	} else if(coloring == "sentrix.col") {
    samples <- sampleChipColumns(object)
  } else if(coloring == "sample.group") {
    samples <- sampleGroups(object)
  } else if(coloring == "chip.id") {
    samples <- sampleChipIDs(object)
  }
  
  cols <- color.function(length(unique(samples)))
	col_vec <- cols[match(samples, unique(samples))]
})


