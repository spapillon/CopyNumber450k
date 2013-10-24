setClass("MethylCNVDataSet", contains = "eSet")

setMethod("initialize", signature("MethylCNVDataSet"), function(.Object, RGSet) {
    if (!is(RGSet, "RGChannelSet")) {
        stop("Expected argument RGSet to be of type minfi::RGChannelSet-class.")
    }
    
    # High-throughput data (AssayData-class/list of matrices)
    mSet <- preprocessRaw(RGSet)
    assayData <- assayDataNew(storage.mode = "lockedEnvironment", intensity = getMeth(mSet) + 
        getUnmeth(mSet))
    featureNames <- featureNames(assayData)
    
    ### Sample covariates (AnnotatedDataFrame-class): phenoData
    phenoData <- AnnotatedDataFrame(data.frame(sex = rep(NA, ncol(assayData$intensity))))
    
    # Feature covariates (AnnotatedDataFrame-class)
    featureData <- AnnotatedDataFrame(data.frame(isUsed = rep(TRUE, length(featureNames)), 
        row.names = featureNames))
    
    ### Experimental description (MIAxE): experimentData
    
    # Assay description (character)
    annotation <- annotation(RGSet)
    
    ### Equipment-generated variables describing sample phenotypes
    ### (AnnotatedDataFrame-class): protocolData
    
    .Object <- callNextMethod(.Object, assayData = assayData, featureData = featureData, 
        annotation = annotation)
}) 
