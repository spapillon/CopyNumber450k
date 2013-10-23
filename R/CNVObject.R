# Class definition

setClass("CNVObject", representation(RGSetSummary = "list", segments = "list", filters = "list", 
    sample_groups = "character", sample_sexes = "character", sample_names = "character", 
    used_probes = "logical", probe_annotation = "data.frame", intensity_matrix = "matrix", 
    is_normalized = "logical", sample_chip_ids = "character", sample_chip_rows = "character", 
    sample_chip_columns = "character"), contains = c("RGChannelSet"))

CNVObject <- function(RGset) new("CNVObject", RGset)



# Initialization methods

setGeneric("filterSNPProbes", function(object) standardGeneric("filterSNPProbes"))

setGeneric("filterVariantProbes", function(object, variance_centile = 0.95) {
    standardGeneric("filterVariantProbes")
})

setGeneric("predictSex", function(object, threshold = -3) standardGeneric("predictSex"))

setGeneric("normalize", function(object, type = c("functional", "quantile")) {
    standardGeneric("normalize")
})

setGeneric("segmentize", function(object, verbose = T, p.adjust.method = "bonferroni", 
    plotting = F) standardGeneric("segmentize"))

setGeneric("createFilters", function(object, tick.threshold = 10, pvalue.threshold = 0.01, 
    breakpoints = c(0, 0, 0, 0)) {
    standardGeneric("createFilters")
})



# Accession methods

setGeneric("intensityMatrix", function(object) standardGeneric("intensityMatrix"))

setGeneric("RGSetSummary", function(object) standardGeneric("RGSetSummary"))

setGeneric("probesAnnotation", function(object) standardGeneric("probesAnnotation"))

setGeneric("segments", function(object) standardGeneric("segments"))

setGeneric("filters", function(object) standardGeneric("filters"))

setGeneric("sampleGroups", function(object) standardGeneric("sampleGroups"))

setGeneric("sampleSexes", function(object) standardGeneric("sampleSexes"))

setGeneric("usedProbes", function(object) standardGeneric("usedProbes"))

setGeneric("isNormalized", function(object) standardGeneric("isNormalized"))

setGeneric("sampleChipRows", function(object) standardGeneric("sampleChipRows"))

setGeneric("sampleChipColumns", function(object) standardGeneric("sampleChipColumns"))

setGeneric("sampleChipIDs", function(object) standardGeneric("sampleChipIDs"))



# Comparison Methods

setGeneric("findCNV", function(object, CNVs, type = "both") standardGeneric("findCNV"))

setGeneric("intersectCNV", function(object, sample_indices, type = "both") standardGeneric("intersectCNV"))

setGeneric("subgroupDifference", function(object, group1_indices, group2_indices) standardGeneric("subgroupDifference"))



# Plotting Methods

setGeneric("plotSample", function(object, index, chr, start, end) standardGeneric("plotSample"))

setGeneric("plotSex", function(object) standardGeneric("plotSex"))

setGeneric("plotDensity", function(object, color.by = c("sentrix.row", "sentrix.col", 
    "sample.group", "chip.id"), color.function = rainbow, legend.position = "topright") standardGeneric("plotDensity"))

setGeneric("plotPCA", function(object, color.by = c("sentrix.row", "sentrix.col", 
    "sample.group", "chip.id"), color.function = rainbow, legend.position = "topright") standardGeneric("plotPCA"))

setGeneric("getColors", function(object, color.by = c("sentrix.row", "sentrix.col", 
    "sample.group", "chip.id"), color.function = rainbow) standardGeneric("getColors"))



# Replacement methods

setGeneric("probesAnnotation<-", function(object, annotation) standardGeneric("probesAnnotation<-"))

setGeneric("usedProbes<-", function(object, probes) standardGeneric("usedProbes<-"))

setGeneric("sampleGroups<-", function(object, groups) standardGeneric("sampleGroups<-"))

setGeneric("sampleSexes<-", function(object, sexes) standardGeneric("sampleSexes<-"))

setGeneric("sampleNames<-", function(object, names) standardGeneric("standardNames<-"))

setGeneric("sampleChipRows<-", function(object, chipRows) standardGeneric("sampleChipRows<-"))

setGeneric("sampleChipColumns<-", function(object, chipColumns) standardGeneric("sampleChipColumns<-"))

setGeneric("sampleChipIDs<-", function(object, chipIDs) standardGeneric("sampleChipIDs<-")) 
