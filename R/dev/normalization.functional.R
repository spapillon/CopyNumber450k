################################################################################ 

# Build the matrix of control probes intensities for the correction of the
# quantile distributions Extraction of the Control matrix
buildControlMatrix450k <- function(extractedData) {
    greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
    controlNames <- names(greenControls)
    
    # Dye bias: (only for Type II probes)
    index <- match("NEGATIVE", controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    dyebiasMatrix <- log2(greenControls.current/redControls.current)
    dyebias <- apply(dyebiasMatrix, 2, median)
    
    # Bisulfite conversion extraction for probe type II:
    index <- match("BISULFITE CONVERSION II", controlNames)
    redControls.current <- redControls[[index]]
    bisulfite2 <- colMeans(redControls.current)
    
    # Bisulfite conversion extraction for probe type I:
    index <- match("BISULFITE CONVERSION I", controlNames)
    redControls.current <- redControls[[index]][7:9, ]
    greenControls.current <- redControls[[index]][1:3, ]
    bisulfite1 <- colMeans(redControls.current + greenControls.current)
    
    # Staining
    index <- match("STAINING", controlNames)
    sg <- greenControls[[index]][3, ]
    sr <- redControls[[index]][1, ]
    
    # Extension
    index <- match("EXTENSION", controlNames)
    redControls.current <- redControls[[index]]
    greenControls.current <- greenControls[[index]]
    extr <- redControls.current[1:2, ]
    extg <- greenControls.current[3:4, ]
    
    # Hybridization should be monitored only in the green channel
    index <- match("HYBRIDIZATION", controlNames)
    h1 <- greenControls[[index]][1, ]
    h2 <- greenControls[[index]][2, ]
    h3 <- greenControls[[index]][3, ]
    
    # Target removal should be low compared to hybridization probes
    index <- match("TARGET REMOVAL", controlNames)
    tar <- greenControls[[index]]
    
    # Non-polymorphic probes
    index <- match("NON-POLYMORPHIC", controlNames)
    npr <- redControls[[index]][1:2, ]
    npg <- greenControls[[index]][3:4, ]
    
    # Specificity II
    index <- match("SPECIFICITY II", controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    spec2g <- colMeans(greenControls.current)
    spec2r <- colMeans(redControls.current)
    spec2ratio <- spec2g/spec2r
    spec2g <- greenControls.current
    spec2r <- redControls.current
    
    # Specificity I
    index <- match("SPECIFICITY I", controlNames)
    greenControls.current <- greenControls[[index]][1:3, ]
    redControls.current <- redControls[[index]][7:9, ]
    spec1g <- greenControls.current
    spec1r <- redControls.current
    greenControls.current <- greenControls[[index]][1:3, ]
    redControls.current <- redControls[[index]][1:3, ]
    ratio1 <- colMeans(redControls.current)/colMeans(greenControls.current)
    greenControls.current <- greenControls[[index]][7:9, ]
    redControls.current <- redControls[[index]][7:9, ]
    ratio2 <- colMeans(greenControls.current)/colMeans(redControls.current)
    spec1ratio <- (ratio1 + ratio2)/2
    
    # Normalization probes:
    index <- match(c("NORM_A"), controlNames)
    normA <- colMeans(redControls[[index]])
    index <- match(c("NORM_T"), controlNames)
    normT <- colMeans(redControls[[index]])
    index <- match(c("NORM_C"), controlNames)
    normC <- colMeans(greenControls[[index]])
    index <- match(c("NORM_G"), controlNames)
    normG <- colMeans(greenControls[[index]])
    
    dyebias2 <- (normC + normG)/(normA + normT)
    
    model.matrix <- cbind(bisulfite1, bisulfite2, t(extg), t(extr), h1, h2, h3, sg, 
        sr, t(npg), t(npr), t(tar), t(spec1g), t(spec1r), t(spec2g), t(spec2r), ratio1, 
        spec1ratio, spec2ratio, ratio2, normA, normC, normT, normG, dyebias2)
    
    oobG <- extractedData$oob$greenOOB
    oobR <- extractedData$oob$redOOB
    model.matrix <- cbind(model.matrix, t(oobG[c(1, 50, 99), ]), oobG[50, ]/oobR[50, 
        ])
    
    # Imputation
    for (colindex in 1:ncol(model.matrix)) {
        column <- model.matrix[, colindex]
        column[is.na(column)] <- mean(column, na.rm = T)
        model.matrix[, colindex] <- column
    }
    # Scaling
    model.matrix <- scale(model.matrix)
    
    # Fixing outliers
    model.matrix[model.matrix > 3] <- 3
    model.matrix[model.matrix < (-3)] <- -3
    model.matrix
}

################################################################################ 

# Return the corrected quantile distributions
returnFit <- function(model.matrix, quantiles, nPCs) {
    quantiles[1, ] <- 0
    quantiles[500, ] <- quantiles[499, ] + 1000
    newQuantiles <- quantiles
    model.matrix <- prcomp(model.matrix)$x[, 1:nPCs]
    
    meanFunction <- apply(quantiles, 1, mean)
    res <- quantiles - meanFunction
    
    
    ### Estimation of UD(s)
    for (i in 1:nrow(quantiles)) {
        model <- lm(res[i, ] ~ model.matrix - 1)
        res[i, ] <- res[i, ] - model$fitted.values
    }
    
    newQuantiles <- meanFunction + res
    newQuantiles
}

################################################################################ 

# Assumes that the predicted sex is 1 and 2. The sex prediction function used by
# default respects this Return the corrected quantile distributions for the
# X-chromosome intensities
returnFitX <- function(model.matrix, quantiles, nPCs, sex) {
    
    quantiles1 <- quantiles[, sex == "Male"]
    model.matrix1 <- model.matrix[sex == "Male", ]
    
    newQuantiles1 <- returnFit(model.matrix = model.matrix1, quantiles = quantiles1, 
        nPCs = nPCs)
    
    quantiles2 <- quantiles[, sex == "Female"]
    model.matrix2 <- model.matrix[sex == "Female", ]
    
    newQuantiles2 <- returnFit(model.matrix = model.matrix2, quantiles = quantiles2, 
        nPCs = nPCs)
    
    newQuantiles <- quantiles
    newQuantiles[, sex == "Male"] <- newQuantiles1
    newQuantiles[, sex == "Female"] <- newQuantiles2
    newQuantiles
}

################################################################################ 

# Normalize a matrix of intensities with the corrected quantile distributions
normalizeByType <- function(intMatrix, newQuantiles) {
    normMatrix <- matrix(NA, nrow(intMatrix), ncol(intMatrix))
    library(preprocessCore)
    n <- 500
    
    normMatrix <- sapply(1:ncol(intMatrix), function(i) {
        crtColumn <- intMatrix[, i]
        crtColumn.reduced <- crtColumn[!is.na(crtColumn)]
        
        # Generation of the corrected intensities:
        target <- sapply(1:499, function(j) {
            start <- newQuantiles[j, i]
            end <- newQuantiles[j + 1, i]
            sequence <- seq(start, end, (end - start)/n)[-n]
        })
        
        target <- as.vector(target)
        result <- normalize.quantiles.use.target(matrix(crtColumn.reduced), target)
        message(paste("Sample ", i, " is normalized"))
        result
    })
}

################################################################################ 

# Main function call for normalization
functionalNormalization <- function(cnMatrix, extractedData, nPCs = 4, predictedSex) {
    probesI <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I")
    probesII <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "II")
    
    # Chr probes:
    locations <- getLocations(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
    autosomal <- names(locations[seqnames(locations) %in% paste0("chr", 1:22)])
    chrY <- names(locations[seqnames(locations) == "chrY"])
    chrX <- names(locations[seqnames(locations) == "chrX"])
    
    probesIGrn <- intersect(probesI$Name[probesI$Color == "Grn"], autosomal)
    probesIRed <- intersect(probesI$Name[probesI$Color == "Red"], autosomal)
    probesII <- intersect(probesII$Name, autosomal)
    
    uProbeNames <- rownames(cnMatrix)
    uProbesIGrn <- intersect(uProbeNames, probesIGrn)
    uProbesIRed <- intersect(uProbeNames, probesIRed)
    uProbesII <- intersect(uProbeNames, probesII)
    uProbesX <- intersect(uProbeNames, chrX)
    uProbesY <- intersect(uProbeNames, chrY)
    
    II <- match(uProbesII, uProbeNames)
    IRed <- match(uProbesIRed, uProbeNames)
    IGreen <- match(uProbesIGrn, uProbeNames)
    X <- match(uProbesX, uProbeNames)
    Y <- match(uProbesY, uProbeNames)
    types <- c("II", "IRed", "IGreen", "X", "Y")
    indList <- list(II, IRed, IGreen, X, Y)
    
    cnList <- list(extractedData$cnQuantiles$II, extractedData$cnQuantiles$IRed, 
        extractedData$cnQuantiles$IGrn, extractedData$cnQuantiles$X, extractedData$cnQuantiles$Y)
    
    model.matrix <- buildControlMatrix450k(extractedData)
    
    messages <- c("type II", "type I Green", "type I Red")
    message("Normalization of the autosomal probes...")
    
    for (i in 1:3) {
        message(paste0("Normalization of the ", messages[i], " probes..."))
        message("Generating the adjusted quantile distributions...")
        
        newQuantiles <- returnFit(model.matrix = model.matrix, quantiles = cnList[[i]], 
            nPCs = nPCs)
        
        message("Normalizing subset of probes...")
        cnMatrix[indList[[i]], ] <- normalizeByType(cnMatrix[indList[[i]], ], newQuantiles)
    }
    
    if (length(X) != 0) {
        message("Normalization of the X-chromosome...")
        newQuantiles <- returnFitX(model.matrix = model.matrix, quantiles = cnList[[4]], 
            nPCs = nPCs, sex = predictedSex)
        cnMatrix[X, ] <- normalizeByType(cnMatrix[X, ], newQuantiles)
    }
    
    if (length(Y) != 0) {
        message("Quantile normalization of the Y-chromosome...")
        cnMatrix[Y, predictedSex == "Male"] <- preprocessCore::normalize.quantiles(cnMatrix[Y, 
            predictedSex == "Male"])
        cnMatrix[Y, predictedSex == "Female"] <- preprocessCore::normalize.quantiles(cnMatrix[Y, 
            predictedSex == "Female"])
    }
    
    message("Normalization done.")
    
    cnMatrix
}

################################################################################  
