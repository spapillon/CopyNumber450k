################################################################################ 

# Extract various statistics from channel and methylation values contained in
# RGChannelSet
extractFromRGChannelSet450k <- function(RGset, intensities) {
    controlType <- c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II", "EXTENSION", 
        "HYBRIDIZATION", "NEGATIVE", "NON-POLYMORPHIC", "NORM_A", "NORM_C", "NORM_G", 
        "NORM_T", "SPECIFICITY I", "SPECIFICITY II", "TARGET REMOVAL", "STAINING")
    
    intensityRed <- getRed(RGset)
    intensityGreen <- getGreen(RGset)
    
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
    #locusNames <- getManifestInfo(RGset, "locusNames")
    TypeI.Red <- getProbeInfo(RGset, type = "I-Red")
    TypeI.Green <- getProbeInfo(RGset, type = "I-Green")
    
    numberQuantiles <- 100
    probs <- 1:numberQuantiles/100
      
    colQuantiles <- function(x, probs) apply(x, 2, function(col) quantile(col, probs = probs, na.rm = T))
    
    oob <- list(greenOOB = colQuantiles(intensityGreen[c(TypeI.Red$AddressA, TypeI.Red$AddressB), ], probs), 
        redOOB = colQuantiles(intensityRed[c(TypeI.Green$AddressA, TypeI.Green$AddressB), ], probs))
    
    
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
    cnQuantiles <- lapply(1:5, function(i) colQuantiles(intensities[indList[[i]], ], probs))

    
    print(sort(sapply(ls(),function(x){object.size(get(x))})))		
    list(cnQuantiles = cnQuantiles, greenControls = greenControls, redControls = redControls,  oob = oob)
}

################################################################################  
