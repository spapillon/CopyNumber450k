################################################################################ 

# Extract various statistics from channel and methylation values contained in
# RGChannelSet
extractFromRGChannelSet450k <- function(RGset, intensities) {
    controlType <- c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "EXTENSION", 
        "HYBRIDIZATION", "NEGATIVE", "NON-POLYMORPHIC", "NORM_A", "NORM_C", "NORM_G", 
        "NORM_T", "SPECIFICITY I", "SPECIFICITY II", "TARGET REMOVAL", "STAINING")
    
    #MSet <- preprocessRaw(RGset)
    intensityRed <- getRed(RGset)
    intensityGreen <- getGreen(RGset)
    #methylated <- getMeth(MSet)
    #unmethylated <- getUnmeth(MSet)
    #betaValues <- getBeta(MSet)
    #mValues <- getM(MSet)
    #intensities <- methylated + unmethylated
    
    # Extraction of the controls
    greenControls = vector("list", length(controlType))
    redControls = vector("list", length(controlType))
    names(greenControls) <- controlType
    names(redControls) <- controlType
    
    for (i in 1:length(controlType)) {
        if (controlType[i] != "STAINING") {
            ctrlAddress <- getControlAddress(RGset, controlType = controlType[i])
        } else {
            ctrlAddress <- getControlAddress(RGset, controlType = controlType[i])[c(2, 
                3, 4, 6)]
        }
        redControls[[i]] <- intensityRed[ctrlAddress, ]
        greenControls[[i]] <- intensityGreen[ctrlAddress, ]
    }
    
    # Extraction of undefined negative control probes
    locusNames <- getManifestInfo(RGset, "locusNames")
    TypeI.Red <- getProbeInfo(RGset, type = "I-Red")
    TypeI.Green <- getProbeInfo(RGset, type = "I-Green")
    
    numberQuantiles <- 100
    probs <- 1:numberQuantiles/100
    
    greenOOB <- rbind(intensityGreen[TypeI.Red$AddressA, ], intensityGreen[TypeI.Red$AddressB, 
        ])
    redOOB <- rbind(intensityRed[TypeI.Green$AddressA, ], intensityRed[TypeI.Green$AddressB, 
        ])
    
    greenOOB <- apply(greenOOB, 2, function(x) quantile(x, probs = probs, na.rm = T))
    redOOB <- apply(redOOB, 2, function(x) quantile(x, probs = probs, na.rm = T))
    oob <- list(greenOOB = greenOOB, redOOB = redOOB)
    
    # Defining the Type I, II Green and II Red probes;
    probesI <- getProbeInfo(RGset, type = "I")
    probesII <- getProbeInfo(RGset, type = "II")
    
    # Chr probes.
    locations <- getLocations(RGset)
    
    # This is a hack, I coerce the Rle object (seqnames(locations)) into a vector
    # since the base::%in% method has issues dispatching to the correct match()
    # method in the package NAMESPACE autosomal <-
    # names(locations[seqnames(locations) %in% paste0('chr', 1:22)])
    probe_names <- names(locations)
    locations <- as.vector(data.frame(seqnames(locations))[, 1])
    names(locations) <- probe_names
    autosomal <- names(locations)[locations %in% paste0("chr", 1:22)]
    chrY <- names(locations)[locations == "chrY"]
    chrX <- names(locations)[locations == "chrX"]
    # End hack
    
    probesIGrn <- intersect(probesI$Name[probesI$Color == "Grn"], autosomal)
    probesIRed <- intersect(probesI$Name[probesI$Color == "Red"], autosomal)
    probesII <- intersect(probesII$Name, autosomal)
    
    uProbeNames <- rownames(intensities)
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
    #mQuantiles <- vector("list", 5)
    #betaQuantiles <- vector("list", 5)
    #methQuantiles <- vector("list", 5)
    #unmethQuantiles <- vector("list", 5)
    cnQuantiles <- vector("list", 5)
    #names(mQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    #names(betaQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    #names(methQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    #names(unmethQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    names(cnQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    
    nq <- 500
    probs <- seq(0, 1, 1/(nq - 1))
    
    for (i in 1:5) {
        #mQuantiles[[i]] <- apply(mValues[indList[[i]], ], 2, function(x) quantile(x, 
        #    probs = probs, na.rm = T))
        #betaQuantiles[[i]] <- apply(betaValues[indList[[i]], ], 2, function(x) quantile(x, 
        #    probs = probs, na.rm = T))
        #methQuantiles[[i]] <- apply(methylated[indList[[i]], ], 2, function(x) quantile(x, 
        #    probs = probs, na.rm = T))
        #unmethQuantiles[[i]] <- apply(unmethylated[indList[[i]], ], 2, function(x) quantile(x, 
        #    probs = probs, na.rm = T))
        cnQuantiles[[i]] <- apply(intensities[indList[[i]], ], 2, function(x) quantile(x, 
            probs = probs, na.rm = T))
    }
    
    #medianXU <- unmethQuantiles$X[250, ]
    #medianXM <- methQuantiles$X[250, ]
    #medianYU <- unmethQuantiles$Y[250, ]
    #medianYM <- methQuantiles$Y[250, ]
    
    #XYMedians <- list(medianXU = medianXU, medianXM = medianXM, medianYU = medianYU, 
    #    medianYM = medianYM)
    
    #list(mQuantiles = mQuantiles, betaQuantiles = betaQuantiles, methQuantiles = methQuantiles, 
    #    unmethQuantiles = unmethQuantiles, cnQuantiles = cnQuantiles, greenControls = greenControls, 
    #    redControls = redControls, XYMedians = XYMedians, oob = oob)
    list(cnQuantiles = cnQuantiles, greenControls = greenControls, redControls = redControls,  oob = oob)
}

################################################################################  
