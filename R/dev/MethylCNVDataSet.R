setClass("MethylCNVDataSet", representation(summary = "list", segments = "list"), 
    contains = "eSet")



setMethod("initialize", signature("MethylCNVDataSet"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
})



MethylCNVDataSetFromRGChannelSet <- function(RGChannelSet) {
    if (!is(RGChannelSet, "RGChannelSet")) {
        stop("Expected argument RGChannelSet to be of type minfi::RGChannelSet-class.")
    }
    
    # Summary of methylation data
    summary <- extractFromRGChannelSet450k(RGChannelSet)
    
    # High-throughput data
    MSet <- preprocessRaw(RGSet)
    intensities <- getMeth(MSet) + getUnmeth(MSet)
    assayData <- assayDataNew(storage.mode = "lockedEnvironment", intensity = intensities)
    sampleNames <- sampleNames(assayData)
    featureNames <- featureNames(assayData)
    
    # Sample covariates
    sexes <- rep(NA, ncol(assayData$intensity))
    names(sexes) <- sampleNames
    
    if (!"Sample_Group" %in% colnames(pData(RGSet))) {
        stop("Expected phenoData(RGChannelSet) to contain the Sample_Group column.")
    } else {
        groups <- pData(RGSet)$Sample_Group
    }
    
    phenoData <- AnnotatedDataFrame(data.frame(sex = sexes, group = groups[sampleNames], 
        row.names = sampleNames))
    
    # Feature covariates
    featuresUsed <- rep(TRUE, length(featureNames))
    
    featureData <- AnnotatedDataFrame(data.frame(isUsed = featuresUsed, row.names = featureNames))
    
    ### Experimental description (MIAxE): experimentData
    
    # Assay description
    annotation <- annotation(RGSet)
    
    ### Equipment-generated variables describing sample phenotypes
    ### (AnnotatedDataFrame-class): protocolData
    
    new("MethylCNVDataSet", summary = summary, assayData = assayData, phenoData = phenoData, 
        featureData = featureData, annotation = annotation)
} 
