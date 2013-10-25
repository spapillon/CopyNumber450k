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

    pheno <- data.frame(row.names = sampleNames)
    if (!"Sample_Group" %in% colnames(pData(RGChannelSet))) {
        stop("Argument RGChannelSet must be presented such that phenoData(RGChannelSet)$Sample_Group exists.")
    } else {
        pheno <- data.frame(pheno, Sample_Group = pData(RGChannelSet)$Sample_Group)
    }
    
    if ("Array" %in% colnames(pData(RGChannelSet))) {
        pheno <- data.frame(pheno, Array= pData(RGChannelSet)$Array)
    }
    
    if ("Slide" %in% colnames(pData(RGChannelSet))) {
        pheno <- data.frame(pheno, Slide = pData(RGChannelSet)$Slide)
    }
    
    if ("Origin" %in% colnames(pData(RGChannelSet))) {
        pheno <- data.frame(pheno, Origin = pData(RGChannelSet)$Origin)
    }
    
    if ("Sample_Sex" %in% colnames(pData(RGChannelSet))) {
        sexes <- pData(RGChannelSet)$Sample_Sex
    } else {
        #TODO: use predict sex
        sexes <- rep(NA, ncol(assayData$intensity))
        # ---
    }
    pheno <- data.frame(pheno, Sample_Sex = sexes)
    phenoData <- AnnotatedDataFrame(pheno)
    
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
