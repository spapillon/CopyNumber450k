MethylCNVDataSetFromRGChannelSet <- function(RGChannelSet) {
    if (!is(RGChannelSet, "RGChannelSet")) {
        stop("Expected argument RGChannelSet to be of type minfi::RGChannelSet-class.")
    }
    
    # Summary of methylation data
    summary <- summarizeRGChannelSet(RGChannelSet)
    
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

summarizeRGChannelSet <- function(RGChannelSet) {
    if (!is(RGChannelSet, "RGChannelSet")) {
        stop("Expected argument RGChannelSet to be of type minfi::RGChannelSet-class.")
    }
    
    controlType <- c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "EXTENSION", 
        "HYBRIDIZATION", "NEGATIVE", "NON-POLYMORPHIC", "NORM_A", "NORM_C", "NORM_G", 
        "NORM_T", "SPECIFICITY I", "SPECIFICITY II", "TARGET REMOVAL", "STAINING")
    
    MSet <- preprocessRaw(RGChannelSet)
    intensityRed <- getRed(RGChannelSet)
    intensityGreen <- getGreen(RGChannelSet)
    methylated <- getMeth(MSet)
    unmethylated <- getUnmeth(MSet)
    betaValues <- getBeta(MSet)
    mValues <- getM(MSet)
    intensities <- methylated + unmethylated
    
    # Extraction of the controls
    greenControls = vector("list", length(controlType))
    redControls = vector("list", length(controlType))
    names(greenControls) <- controlType
    names(redControls) <- controlType
    
    for (i in 1:length(controlType)) {
        if (controlType[i] != "STAINING") {
            ctrlAddress <- getControlAddress(RGChannelSet, controlType = controlType[i])
        } else {
            ctrlAddress <- getControlAddress(RGChannelSet, controlType = controlType[i])[c(2, 
                3, 4, 6)]
        }
        redControls[[i]] <- intensityRed[ctrlAddress, ]
        greenControls[[i]] <- intensityGreen[ctrlAddress, ]
    }
    
    # Extraction of undefined negative control probes
    locusNames <- getManifestInfo(RGChannelSet, "locusNames")
    TypeI.Red <- getProbeInfo(RGChannelSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(RGChannelSet, type = "I-Green")
    
    numberQuantiles <- 100
    probs <- 1:numberQuantiles/100
    
    greenOOB <- rbind(getGreen(RGChannelSet)[TypeI.Red$AddressA, ], getGreen(RGChannelSet)[TypeI.Red$AddressB, 
        ])
    redOOB <- rbind(getRed(RGChannelSet)[TypeI.Green$AddressA, ], getRed(RGChannelSet)[TypeI.Green$AddressB, 
        ])
    
    greenOOB <- apply(greenOOB, 2, function(x) quantile(x, probs = probs, na.rm = T))
    redOOB <- apply(redOOB, 2, function(x) quantile(x, probs = probs, na.rm = T))
    oob <- list(greenOOB = greenOOB, redOOB = redOOB)
    
    # Defining the Type I, II Green and II Red probes;
    probesI <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I")
    probesII <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "II")
    
    # Chr probes.
    locations <- getLocations(IlluminaHumanMethylation450kannotation.ilmn.v1.2) # TODO: getLocation(RGSet) does not work
    autosomal <- names(locations[seqnames(locations) %in% paste0("chr", 1:22)])
    chrY <- names(locations[seqnames(locations) == "chrY"])
    chrX <- names(locations[seqnames(locations) == "chrX"])
    
    probesIGrn <- intersect(probesI$Name[probesI$Color == "Grn"], autosomal)
    probesIRed <- intersect(probesI$Name[probesI$Color == "Red"], autosomal)
    probesII <- intersect(probesII$Name, autosomal)
    
    uProbeNames <- rownames(betaValues)
    uProbesIGrn <- intersect(uProbeNames, probesIGrn)
    uProbesIRed <- intersect(uProbeNames, probesIRed)
    uProbesII <- intersect(uProbeNames, probesII)
    uProbesX <- intersect(uProbeNames, chrX)
    uProbesY <- intersect(uProbeNames, chrY)
    indicesIGrn <- match(uProbesIGrn, uProbeNames)
    indicesIRed <- match(uProbesIRed, uProbeNames)
    indicesII <- match(uProbesII, uProbeNames)
    indicesX <- match(uProbesX, uProbeNames)
    indicesY <- match(uProbesY, uProbeNames)
    
    indList <- list(indicesIGrn, indicesIRed, indicesII, indicesX, indicesY)
    names(indList) <- c("IGrn", "IRed", "II", "X", "Y")
    
    # Quantile extraction
    mQuantiles <- vector("list", 5)
    betaQuantiles <- vector("list", 5)
    methQuantiles <- vector("list", 5)
    unmethQuantiles <- vector("list", 5)
    cnQuantiles <- vector("list", 5)
    names(mQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    names(betaQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    names(methQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    names(unmethQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    names(cnQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    
    nq <- 500
    probs <- seq(0, 1, 1/(nq - 1))
    
    for (i in 1:5) {
        mQuantiles[[i]] <- apply(mValues[indList[[i]], ], 2, function(x) quantile(x, 
            probs = probs, na.rm = T))
        betaQuantiles[[i]] <- apply(betaValues[indList[[i]], ], 2, function(x) quantile(x, 
            probs = probs, na.rm = T))
        methQuantiles[[i]] <- apply(methylated[indList[[i]], ], 2, function(x) quantile(x, 
            probs = probs, na.rm = T))
        unmethQuantiles[[i]] <- apply(unmethylated[indList[[i]], ], 2, function(x) quantile(x, 
            probs = probs, na.rm = T))
        cnQuantiles[[i]] <- apply(intensities[indList[[i]], ], 2, function(x) quantile(x, 
            probs = probs, na.rm = T))
    }
    
    medianXU <- unmethQuantiles$X[250, ]
    medianXM <- methQuantiles$X[250, ]
    medianYU <- unmethQuantiles$Y[250, ]
    medianYM <- methQuantiles$Y[250, ]
    
    XYMedians <- list(medianXU = medianXU, medianXM = medianXM, medianYU = medianYU, 
        medianYM = medianYM)
    
    list(mQuantiles = mQuantiles, betaQuantiles = betaQuantiles, methQuantiles = methQuantiles, 
        unmethQuantiles = unmethQuantiles, cnQuantiles = cnQuantiles, greenControls = greenControls, 
        redControls = redControls, XYMedians = XYMedians, oob = oob)
} 
