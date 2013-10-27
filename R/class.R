################################################################################ 

setClass("MethylCNVDataSet", representation(summary = "list", segments = "list", 
    manifest = "IlluminaMethylationManifest"), contains = "eSet")

################################################################################ 

setMethod("initialize", signature("MethylCNVDataSet"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
})

################################################################################ 

setMethod("getManifest", signature("MethylCNVDataSet"), function(object) {
    object@manifest
})

################################################################################ 

setGeneric("getSummary", function(object) standardGeneric("getSummary"))

setMethod("getSummary", signature("MethylCNVDataSet"), function(object) {
    object@summary
})

################################################################################ 

setGeneric("getSegments", function(object) standardGeneric("getSegments"))

setMethod("getSegments", signature("MethylCNVDataSet"), function(object) {
    object@segments
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
    
    # Sample covariates - daisy chain the columns
    if ("Sample_Sex" %in% colnames(pData(RGChannelSet))) {
        sexes <- pData(RGChannelSet)$Sample_Sex
    } else {
        sexes <- rep(NA, ncol(assayData$intensity))
    }
    
    phenoData <- data.frame(Sample_Sex = sexes, row.names = sampleNames)
    
    if (!"Sample_Group" %in% colnames(pData(RGChannelSet))) {
        stop("Argument RGChannelSet must be presented such that phenoData(RGChannelSet)$Sample_Group exists.")
    } else {
        phenoData <- data.frame(phenoData, Sample_Group = pData(RGChannelSet)$Sample_Group)
    }
    
    if ("Array" %in% colnames(pData(RGChannelSet))) {
        phenoData <- data.frame(phenoData, Array = pData(RGChannelSet)$Array)
    }
    
    if ("Slide" %in% colnames(pData(RGChannelSet))) {
        phenoData <- data.frame(phenoData, Slide = pData(RGChannelSet)$Slide)
    }
    
    if ("Origin" %in% colnames(pData(RGChannelSet))) {
        phenoData <- data.frame(phenoData, Origin = pData(RGChannelSet)$Origin)
    }
    
    phenoData <- AnnotatedDataFrame(phenoData)
    
    # Feature covariates
    featureData <- as.data.frame(getAnnotation(RGChannelSet, what = c("Locations", 
        "SNPs.137CommonSingle", "Other")))
    rownames(featureData) <- featureNames
    featureData <- AnnotatedDataFrame(featureData)
    
    # TODO: Experimental description (MIAxE): experimentData
    
    # Assay description
    annotation <- annotation(RGChannelSet)
    manifest <- getManifest(RGChannelSet)
    
    # TODO: Equipment-generated variables describing sample phenotypes
    # (AnnotatedDataFrame-class): protocolData
    
    object <- new("MethylCNVDataSet", summary = summary, assayData = assayData, phenoData = phenoData, 
        featureData = featureData, annotation = annotation, manifest = manifest, 
        segments = list())

    # TODO: Weird to put it here... should predictSampleSex just be a general function and not a method?
    pData(object)$Sample_Sex <- predictSampleSexes(object)
    object
}

################################################################################  
