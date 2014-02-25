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
        column[is.na(column)] <- mean(column, na.rm = TRUE)
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
    
    
    for (i in 1:nrow(quantiles)) {
        model <- lm(res[i, ] ~ model.matrix - 1)
        res[i, ] <- res[i, ] - model$fitted.values
    }
    
    newQuantiles <- meanFunction + res
    newQuantiles
}

################################################################################ 





################################################################################ 

# Normalize a matrix of intensities with the corrected quantile distributions
normalizeByType <- function(intMatrix, newQuantiles) {
    normMatrix <- matrix(NA, nrow(intMatrix), ncol(intMatrix))
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
        result <- preprocessCore::normalize.quantiles.use.target(matrix(crtColumn.reduced), 
            target)
        result
    })
}

################################################################################ 

# Main function call for normalization
functionalNormalization <- function(cnMatrix, extractedData, manifest, nPCs = 4) {
    
    probesI <- getProbeInfo(manifest, type = "I")
    probesII <- getProbeInfo(manifest, type = "II")
    probesIGrn <- probesI$Name[probesI$Color == "Grn"]
    probesIRed <- probesI$Name[probesI$Color == "Red"]
    
    uProbeNames <- rownames(cnMatrix)
    uProbesIGrn <- intersect(uProbeNames, probesIGrn)
    uProbesIRed <- intersect(uProbeNames, probesIRed)
    uProbesII <- intersect(uProbeNames, probesII$Name)
    
    II <- match(uProbesII, uProbeNames)
    IRed <- match(uProbesIRed, uProbeNames)
    IGreen <- match(uProbesIGrn, uProbeNames)
    
    
    types <- c("II", "IRed", "IGreen")
    indList <- list(II, IRed, IGreen)
    messages <- c("type II", "type I Green", "type I Red")
    
    
    model.matrix <- buildControlMatrix450k(extractedData)
    cnList <- list(extractedData$cnQuantiles$II, extractedData$cnQuantiles$IRed, 
        extractedData$cnQuantiles$IGrn)
    
    ### Normalization
    for (i in 1:3) {
        message(paste0("Normalization of the ", messages[i], " probes..."))
        
        newQuantiles <- returnFit(model.matrix = model.matrix, quantiles = cnList[[i]], 
            nPCs = nPCs)
        
        cnMatrix[indList[[i]], ] <- normalizeByType(cnMatrix[indList[[i]], ], newQuantiles)
    }
    
    message("Normalization done.")
    
    cnMatrix
}

################################################################################  
