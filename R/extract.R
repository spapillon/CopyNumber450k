################################################################################ 

extractControls <- function(RGset, fun) {
    controlType <- c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "EXTENSION", 
        "HYBRIDIZATION", "NEGATIVE", "NON-POLYMORPHIC", "NORM_A", "NORM_C", "NORM_G", 
        "NORM_T", "SPECIFICITY I", "SPECIFICITY II", "TARGET REMOVAL", "STAINING")
    x <- lapply(controlType, function(i) {
      ctrlAddress <- getControlAddress(RGset, controlType = i)
      fun(RGset)[ctrlAddress, ]
    })
    names(x) <- controlType
    x
}

################################################################################ 

# Extract various statistics from channel and methylation values contained in
# RGChannelSet
extractFromRGChannelSet450k <- function(RGset, intensities) {
    # Extraction of the controls
    greenControls <- extractControls(RGset, getGreen)
    redControls <- extractControls(RGset, getRed)
    
    # Extraction of undefined negative control probes
    #locusNames <- getManifestInfo(RGset, "locusNames")
    TypeI.Red <- getProbeInfo(RGset, type = "I-Red")
    TypeI.Green <- getProbeInfo(RGset, type = "I-Green")
    
    numberQuantiles <- 100
    probs <- 1:numberQuantiles/100

    greenOOB <- getGreen(RGset)[c(TypeI.Red$AddressA, TypeI.Red$AddressB), ]
    greenOOB <- apply(greenOOB, 2, function(x) quantile(x, probs = probs, na.rm = TRUE))

    redOOB <- getRed(RGset)[c(TypeI.Green$AddressA, TypeI.Green$AddressB), ]
    redOOB <- apply(redOOB, 2, function(x) quantile(x, probs = probs, na.rm = TRUE))

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
    
    uProbeNames <- rownames(intensities)
    indList <- list("IGrn" =  match(intersect(probesI$Name[probesI$Color == "Grn"], autosomal), uProbeNames), 
            "IRed" = match(intersect(probesI$Name[probesI$Color == "Red"], autosomal), uProbeNames), 
            "II" = match(intersect(probesII$Name, autosomal), uProbeNames), 
            "X" = match(chrX, uProbeNames), 
            "Y" = match(chrY, uProbeNames))
    
    # Quantile extraction
    cnQuantiles <- vector("list", 5)
    names(cnQuantiles) <- c("IGrn", "IRed", "II", "X", "Y")
    
    nq <- 500
    probs <- seq(0, 1, 1/(nq - 1))
    
    for (i in 1:5) {
        cnQuantiles[[i]] <- apply(intensities[indList[[i]], ], 2, 
                function(x) quantile(x, probs = probs, na.rm = TRUE))
    }

    list(cnQuantiles = cnQuantiles, greenControls = greenControls, 
         redControls = redControls,  oob = oob)
}

################################################################################  
