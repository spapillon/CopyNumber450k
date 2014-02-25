################################################################################ 

setMethod("getManifest", signature("CNV450kSet"), function(object) {
    object@manifest
})

################################################################################ 

setMethod("getSummary", signature("CNV450kSet"), function(object) {
    object@summary
})

################################################################################ 

setMethod("getSegments", signature("CNV450kSet"), function(object) {
    object@segments
})

################################################################################  
