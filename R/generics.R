# TODO: Add comment
# 
# Author: spapillo
###############################################################################

# Object building
setGeneric("filterSNPProbes", function(object) standardGeneric("filterSNPProbes"))
setGeneric("filterVariantProbes", function(object, variance_centile=0.95) standardGeneric("filterVariantProbes"))
setGeneric("normalize", function(object, type=c("functional", "quantile")) standardGeneric("normalize"))
setGeneric("buildSegments", function(object, verbose=T) standardGeneric("buildSegments"))
setGeneric("createFilters", function(object,  tick.threshold=10, pvalue.threshold=0.01, breakpoints=c(0,0,0,0)) standardGeneric("createFilters"))

# Workable methods
setGeneric("findCNV", function(object, CNVs, type="both") standardGeneric("findCNV"))
setGeneric("intersectCNV", function(object, sample_indices, type="both") standardGeneric("intersectCNV"))
setGeneric("subgroupDifference", function(object, group1_indices, group2_indices) standardGeneric("subgroupDifference"))
setGeneric("plotSex", function(object) standardGeneric("plotSex"))
setGeneric("plotSample", function(object, index, chr, start, end) standardGeneric("plotSample"))

setGeneric("predictSex", function(object, threshold=-3) standardGeneric("predictSex"))


# Accession methods
setGeneric("intensityMatrix", function(object) standardGeneric("intensityMatrix"))
setGeneric("RGSetSummary",function(object) standardGeneric("RGSetSummary"))
setGeneric("probesAnnotation", function(object) standardGeneric("probesAnnotation"))
setGeneric("segments", function(object) standardGeneric("segments"))
setGeneric("filters", function(object) standardGeneric("filters"))
setGeneric("sampleGroups",  function(object) standardGeneric("sampleGroups"))
setGeneric("sampleSexes", function(object) standardGeneric("sampleSexes"))

setGeneric("usedProbes", function(object) standardGeneric("usedProbes"))
setGeneric("isNormalized", function(object) standardGeneric("isNormalized"))

# Replacement methods
setGeneric("probesAnnotation<-", function(object, value) standardGeneric("probesAnnotation<-"))
setGeneric("sampleGroups<-", function(object, value) standardGeneric("sampleGroups<-"))
setGeneric("sampleSexes<-", function(object, value) standardGeneric("sampleSexes<-"))

setGeneric("usedProbes<-",  function(object, value) standardGeneric("usedProbes<-"))