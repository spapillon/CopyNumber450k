#
# Class definition
#

setClass("CNVObject",
         representation(RGSetSummary     = "list",
                        segments         = "list",
                        filters          = "list",
                        sample_groups    = "character",
                        sample_sexes     = "character",
                        sample_names     = "character",
                        used_probes      = "logical",
                        probe_annotation = "data.frame",
                        intensity_matrix = "matrix",
                        is_normalized    = "logical",
                        sample_chip_ids = "character",
                        sample_chip_rows = "character",
                        sample_chip_columns = "character"),
         contains = c("RGChannelSet"))

CNVObject <- function(RGset) new("CNVObject", RGset)



#
# Initialization methods
#

setGeneric("filterSNPProbes",
           function(object) standardGeneric("filterSNPProbes"))

setGeneric("filterVariantProbes",
           function(object,
                    variance_centile = 0.95) {
             standardGeneric("filterVariantProbes")
           })

setGeneric("normalize",
           function(object,
                    type = c("functional", "quantile")) {
             standardGeneric("normalize")
           })

setGeneric("segmentize",
           function(object,
                    verbose = T,
                    p.adjust.method = "bonferroni",
                    plotting = F)
            standardGeneric("segmentize"))

setGeneric("createFilters",
           function(object,
                    tick.threshold   = 10,
                    pvalue.threshold = 0.01,
                    breakpoints      = c(0, 0, 0, 0)) {
             standardGeneric("createFilters")
           })



#
# Accession methods
#
       
setGeneric("intensityMatrix",
           function(object) standardGeneric("intensityMatrix"))
 
setGeneric("RGSetSummary",
           function(object) standardGeneric("RGSetSummary"))
 
setGeneric("probesAnnotation",
           function(object) standardGeneric("probesAnnotation"))
 
setGeneric("segments",
           function(object) standardGeneric("segments"))
 
setGeneric("filters",
           function(object) standardGeneric("filters"))
 
setGeneric("sampleGroups",
           function(object) standardGeneric("sampleGroups"))
 
setGeneric("sampleSexes",
           function(object) standardGeneric("sampleSexes"))
 
setGeneric("usedProbes",
           function(object) standardGeneric("usedProbes"))
 
setGeneric("isNormalized",
           function(object) standardGeneric("isNormalized"))
       
setGeneric("sampleChipRows",
           function(object) standardGeneric("sampleChipRows"))
       
setGeneric("sampleChipColumns",
           function(object) standardGeneric("sampleChipColumns"))
       
setGeneric("sampleChipIDs",
           function(object) standardGeneric("sampleChipIDs"))


#
# Comparison Methods
#

setGeneric("findCNV",
           function(object,
                    CNVs,
                    type = "both") standardGeneric("findCNV"))

setGeneric("intersectCNV",
           function(object,
                    sample_indices,
                    type = "both") standardGeneric("intersectCNV"))

setGeneric("subgroupDifference",
           function(object,
                    group1_indices,
                    group2_indices) standardGeneric("subgroupDifference"))



#
# Plotting Methods
#

setGeneric("plotSex",
           function(object) standardGeneric("plotSex"))

setGeneric("plotSample",
           function(object,
                    index,
                    chr,
                    start,
                    end) standardGeneric("plotSample"))
            
setGeneric("predictSex",
           function(object,
                    threshold = -3) standardGeneric("predictSex"))
            
setGeneric("plotDensity",
           function(object,
                    color.by = c("sentrix.row", "sentrix.col", "sample.group", "chip.id"),
                    color.function = rainbow) standardGeneric("plotDensity"))
            
setGeneric("getColors",
           function(object,
                    color.by = c("sentrix.row", "sentrix.col", "sample.group", "chip.id"),
                    color.function = rainbow) standardGeneric("getColors"))


#
# Replacement methods
#

setGeneric("probesAnnotation<-",
           function(object, value) standardGeneric("probesAnnotation<-"))

setGeneric("sampleGroups<-",
           function(object, value) standardGeneric("sampleGroups<-"))

setGeneric("sampleSexes<-",
           function(object, value) standardGeneric("sampleSexes<-"))

setGeneric("usedProbes<-",
           function(object, value) standardGeneric("usedProbes<-"))
       
setGeneric("sampleChipRows<-",
           function(object, value) standardGeneric("sampleChipRows"))
       
setGeneric("sampleChipColumns<-",
           function(object, value) standardGeneric("sampleChipColumns"))
       
setGeneric("sampleChipIDs<-",
           function(object, value) standardGeneric("sampleChipIDs"))
