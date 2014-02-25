################################################################################ 

# Main function call for quantile normalization
quantileNormalization <- function(cnMatrix, manifest) {
    
    probesI <- getProbeInfo(manifest, type = "I")
    probesII <- getProbeInfo(manifest, type = "II")
    probesIGrn <- probesI$Name[probesI$Color == "Grn"]
    probesIRed <- probesI$Name[probesI$Color == "Red"]
    
    uProbeNames <- rownames(cnMatrix)
    uProbesIGrn <- intersect(uProbeNames, probesIGrn)
    uProbesIRed <- intersect(uProbeNames, probesIRed)
    uProbesII <- intersect(uProbeNames, probesII$Name)
    
    II <- match(uProbesII, uProbeNames)
    IRed <- match(uProbesIRed, uProbeNames)
    IGreen <- match(uProbesIGrn, uProbeNames)
    
    cnMatrix[II, ] <- preprocessCore::normalize.quantiles(cnMatrix[II, ])
    cnMatrix[IRed, ] <- preprocessCore::normalize.quantiles(cnMatrix[IRed, ])
    cnMatrix[IGreen, ] <- preprocessCore::normalize.quantiles(cnMatrix[IGreen, ])
    
    message("Quantile Normalization done.")
    
    cnMatrix
}

################################################################################  
