################################################################################ 

setMethod("initialize", signature("CNV450kSet"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
})

################################################################################ 

CNV450kSet <- function(RGChannelSet) {
    if (!is(RGChannelSet, "RGChannelSet")) {
        stop("Argument RGChannelSet must be of type minfi::RGChannelSet-class.")
    }
    
    if ("Sample_Name" %in% colnames(pData(RGChannelSet))) {
        sampleNames(RGChannelSet) <- pData(RGChannelSet)$Sample_Name
    }
    # High-throughput data
    MSet <- preprocessRaw(RGChannelSet)
    intensities <- getMeth(MSet) + getUnmeth(MSet)
    assayData <- assayDataNew(storage.mode = "lockedEnvironment", intensity = intensities)
    sampleNames <- sampleNames(assayData)
    featureNames <- featureNames(assayData)
        
    # Summary of methylation data
    summary <- extractFromRGChannelSet450k(RGChannelSet, intensities)
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
    
    # Experimental description (MIAxE)
    experimentData <- experimentData(RGChannelSet)
    
    # Assay description
    annotation <- annotation(RGChannelSet)
    manifest <- getManifest(RGChannelSet)
    
    # Equipment-generated variables describing sample phenotypes
    protocolData <- protocolData(RGChannelSet)
    
    new("CNV450kSet", summary = summary, assayData = assayData, phenoData = phenoData, 
        featureData = featureData, annotation = annotation, manifest = manifest, 
        experimentData = experimentData, protocolData = protocolData, segments = list())
    
}

################################################################################  
