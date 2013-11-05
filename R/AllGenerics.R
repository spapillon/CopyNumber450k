################################################################################  

setGeneric("getSummary", function(object) standardGeneric("getSummary"))

################################################################################  

setGeneric("getSegments", function(object) standardGeneric("getSegments"))

################################################################################  

setGeneric("findCNV", function(object, gene_names, type = "both") standardGeneric("findCNV"))

################################################################################  

setGeneric("intersectCNV", function(object, sample_indices, type = "both") standardGeneric("intersectCNV"))

################################################################################  

setGeneric("subgroupDifference", function(object, group1_indices, group2_indices) standardGeneric("subgroupDifference"))

################################################################################  

setGeneric("dropSNPprobes", function(object, maf_threshold = 0) standardGeneric("dropSNPprobes"))

################################################################################  

setGeneric("normalize", function(object, ...) standardGeneric("normalize", ...))

################################################################################  

setGeneric("segmentize", function(object, verbose = TRUE, p.adjust.method = "bonferroni") standardGeneric("segmentize"))

################################################################################  

setGeneric("computeSignificance", function(object, p.value.threshold = 0.01, num.mark.threshold = 10)
            standardGeneric("computeSignificance"))

################################################################################  

setGeneric("getColoring", function(object, color.by = c("array.row", "array.col", 
                        "sample.group", "slide", "origin"), color.function = rainbow) standardGeneric("getColoring"))

################################################################################  

setGeneric("plotSample", function(object, index, chr, start, end) standardGeneric("plotSample"))

################################################################################  

setGeneric("plotDensity", function(object, color.by = c("array.row", "array.col", 
                        "sample.group", "slide", "origin"), color.function = rainbow, legend.position = "topright") standardGeneric("plotDensity"))

################################################################################  

setGeneric("plotPCA", function(object, color.by = c("array.row", "array.col", "sample.group", 
                        "slide", "origin"), color.function = rainbow, legend.position = "topright") standardGeneric("plotPCA"))

################################################################################  