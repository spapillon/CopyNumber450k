################################################################################ 

# Main function call for quantile normalization
quantileNormalization <- function(cnMatrix, annotation, manifest, predictedSex = NULL) {
    
    probesI <- getProbeInfo(manifest, type = "I")
    probesII <- getProbeInfo(manifest, type = "II")
    
    # Chr probes:
    locations <- getLocations(minfi:::.getAnnotationString(annotation))
    autosomal <- names(locations[seqnames(locations) %in% paste0("chr", 1:22)])
    chrY <- names(locations[seqnames(locations) == "chrY"])
    chrX <- names(locations[seqnames(locations) == "chrX"])
    
    probesIGrn <- intersect(probesI$Name[probesI$Color == "Grn"], autosomal)
    probesIRed <- intersect(probesI$Name[probesI$Color == "Red"], autosomal)
    probesII <- intersect(probesII$Name, autosomal)
    
    uProbeNames <- rownames(cnMatrix)
    uProbesIGrn <- intersect(uProbeNames, probesIGrn)
    uProbesIRed <- intersect(uProbeNames, probesIRed)
    uProbesII <- intersect(uProbeNames, probesII)
    uProbesX <- intersect(uProbeNames, chrX)
    uProbesY <- intersect(uProbeNames, chrY)
    
    II <- match(uProbesII, uProbeNames)
    IRed <- match(uProbesIRed, uProbeNames)
    IGreen <- match(uProbesIGrn, uProbeNames)
    X <- match(uProbesX, uProbeNames)
    Y <- match(uProbesY, uProbeNames)
    
    cnMatrix[II, ] <- preprocessCore::normalize.quantiles(cnMatrix[II, ])
    cnMatrix[IRed, ] <- preprocessCore::normalize.quantiles(cnMatrix[IRed, ])
    cnMatrix[IGreen, ] <- preprocessCore::normalize.quantiles(cnMatrix[IGreen, ])
    
    if (!is.null(predictedSex)) {
        cnMatrix[Y, predictedSex == "Male"] <- preprocessCore::normalize.quantiles(cnMatrix[Y, 
            predictedSex == "Male"])
        cnMatrix[Y, predictedSex == "Female"] <- preprocessCore::normalize.quantiles(cnMatrix[Y, 
            predictedSex == "Female"])
        cnMatrix[X, predictedSex == "Male"] <- preprocessCore::normalize.quantiles(cnMatrix[X, 
            predictedSex == "Male"])
        cnMatrix[X, predictedSex == "Female"] <- preprocessCore::normalize.quantiles(cnMatrix[X, 
            predictedSex == "Female"])
    } else {
        cnMatrix[Y, ] <- preprocessCore::normalize.quantiles(cnMatrix[Y, ])
        cnMatrix[X, ] <- preprocessCore::normalize.quantiles(cnMatrix[X, ])
    }
    
    message("Quantile Normalization done.")
    
    cnMatrix
}

################################################################################  
