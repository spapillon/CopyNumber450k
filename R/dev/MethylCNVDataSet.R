setClass("MethylCNVDataSet", contains = "eSet")

setMethod("initialize", signature("MethylCNVDataSet"), function(.Object, RGSet) {
    if (!is(RGSet, "RGChannelSet")) {
        stop("Expected argument RGSet to be of type minfi::RGChannelSet-class.")
    }
    
    # High-throughput data (AssayData-class/list of matrices)
    mSet <- preprocessRaw(RGSet)
    assayData <- assayDataNew(storage.mode = "lockedEnvironment", intensity = getMeth(mSet) + 
        getUnmeth(mSet))
    
    ### Sample covariates (AnnotatedDataFrame-class): phenoData
    
    # Feature covariates (AnnotatedDataFrame-class)
    featureData <- AnnotatedDataFrame(data.frame(isUsed = rep(TRUE, nrow(assayData$intensity))))
    ### fvarMetadata / fvarLabels
    
    ### Experimental description (MIAxE): experimentData
    
    # Assay description (character)
    annotation <- annotation(RGSet)
    
    ### Equipment-generated variables describing sample phenotypes
    ### (AnnotatedDataFrame-class): protocolData
    
    .Object <- callNextMethod(.Object, assayData = assayData, featureData = featureData, 
        annotation)
}) 
