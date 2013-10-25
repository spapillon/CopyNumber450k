################################################################################ 

setClass("MethylCNVDataSet", representation(summary = "list", segments = "list"), 
    contains = "eSet")

################################################################################ 

setMethod("initialize", signature("MethylCNVDataSet"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
})

################################################################################ 

MethylCNVDataSetFromRGChannelSet <- function(RGChannelSet) {
    if (!is(RGChannelSet, "RGChannelSet")) {
        stop("Argument RGChannelSet must be of type minfi::RGChannelSet-class.")
    }
    
    # Summary of methylation data
    summary <- extractFromRGChannelSet450k(RGChannelSet)
    
    # High-throughput data
    MSet <- preprocessRaw(RGChannelSet)
    intensities <- getMeth(MSet) + getUnmeth(MSet)
    assayData <- assayDataNew(storage.mode = "lockedEnvironment", intensity = intensities)
    sampleNames <- sampleNames(assayData)
    featureNames <- featureNames(assayData)
    
    # Sample covariates
    sexes <- rep(NA, ncol(assayData$intensity))
    names(sexes) <- sampleNames
    
    if (!"Sample_Group" %in% colnames(pData(RGChannelSet))) {
        stop("Argument RGChannelSet must be presented such that phenoData(RGChannelSet)$Sample_Group exists.")
    } else {
        groups <- pData(RGChannelSet)$Sample_Group
    }
    
    phenoData <- AnnotatedDataFrame(data.frame(sex = sexes[sampleNames], group = groups[sampleNames], 
        row.names = sampleNames))
    
    # Feature covariates
    featureData <- AnnotatedDataFrame(data.frame(getAnnotation(RGChannelSet, what=c("Locations", "SNPs.137CommonSingle"))))
    # featureData <- AnnotatedDataFrame(data.frame(isUsed = featuresUsed, row.names =
    # featureNames))
    
    ### Experimental description (MIAxE): experimentData
    
    # Assay description
    annotation <- annotation(RGChannelSet)
    
    ### Equipment-generated variables describing sample phenotypes
    ### (AnnotatedDataFrame-class): protocolData
    new("MethylCNVDataSet", summary = summary, assayData = assayData, phenoData = phenoData, 
        featureData = featureData, annotation = annotation)
}

################################################################################  
