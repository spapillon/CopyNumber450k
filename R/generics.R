# TODO: Add comment
# 
# Author: spapillo
###############################################################################

# Object building
setGeneric("filterSNPProbes",function(object) standardGeneric("filterSNPProbes"))
setGeneric("normalize", function(object, sex_cutoff=-3) standardGeneric("normalize"))
setGeneric("buildSegments", function(object, verbose=T) standardGeneric("buildSegments"))
setGeneric("createFilters", function(object,  tick.threshold=50, pvalue.threshold=0.01, breakpoints=c(0,0,0,0)) standardGeneric("createFilters"))

# Workable methods
setGeneric("findCNV", function(object, CNVs, type="both") standardGeneric("findCNV"))
setGeneric("intersectCNV", function(object, sample_indices, type="both") standardGeneric("intersectCNV"))
setGeneric("subgroupDifference", function(object, group1_indices, group2_indices) standardGeneric("subgroupDifference"))

# Accession methods
setGeneric("intensityMatrix", function(object) standardGeneric("intensityMatrix"))
setGeneric("RGSetSummary",function(object) standardGeneric("RGSetSummary"))
setGeneric("probesAnnotation", function(object) standardGeneric("probesAnnotation"))
setGeneric("segments", function(object) standardGeneric("segments"))
setGeneric("filters", function(object) standardGeneric("filters"))
setGeneric("sampleGroups",  function(object) standardGeneric("sampleGroups"))
setGeneric("usedProbes", function(object) standardGeneric("usedProbes"))
setGeneric("isNormalized", function(object) standardGeneric("isNormalized"))

# Replacement methods
setGeneric("probesAnnotation<-", function(object, value) standardGeneric("probesAnnotation<-"))
setGeneric("sampleGroups<-", function(object, value) standardGeneric("sampleGroups<-"))
setGeneric("usedProbes<-",  function(object, value) standardGeneric("usedProbes<-"))