######################################################
## Functional normalization of the 450k array
## Jean-Philippe Fortin 
## Aug 30 2013
## Code in development
#####################################################

library(fda)
library(refund)
library(preprocessCore) 
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)

### Plot the gender clusters
#################################################################
plotSex <- function(extractedData){
	plot(log2(extractedData$cnQuantiles$Y[250,]) - 
	                   log2(extractedData$cnQuantiles$X[250,]), rep(0,length(extractedData$cnQuantiles$Y[250,])))
}
#################################################################

### Predict the sex of the samples
#################################################################
predictSex <- function(extractedData, cutoff){
	diffs<- log2(extractedData$cnQuantiles$Y[250,]) - 
	                   log2(extractedData$cnQuantiles$X[250,])
	predictedSex <- rep(1, length(extractedData$cnQuantiles$Y[250,]))
	predictedSex[which(diffs < cutoff)] <- 2
	return(predictedSex)
}   
#################################################################

### Build the dye bias covariate with the internal normalization probes
#################################################################
buildDyebias <- function(extractedData) {
	
	greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
	controlNames <- names(greenControls)           
	 
	index <- match(c("NORM_A"), controlNames)
	normA <- colMeans(redControls[[index]])
	index <- match(c("NORM_T"), controlNames)
	normT <- colMeans(redControls[[index]])
	index <- match(c("NORM_C"), controlNames)
	normC <- colMeans(greenControls[[index]])
	index <- match(c("NORM_G"), controlNames)
	normG <- colMeans(greenControls[[index]])
    
    dyebias <- log((normC + normG)/(normA+normT))
          
    return(dyebias)
}
#################################################################

### Build the matrix of control probes intensities for the
### correction of the quantile distributions 
#################################################################
### Extraction of the Control matrix
buildControlMatrix450k <- function(extractedData) {
	 
    greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
	controlNames <- names(greenControls)                                 	
	 	
	# Dye bias: (only for Type II probes)
	index <- match("NEGATIVE",controlNames)
	greenControls.current <- greenControls[[index]]
	redControls.current <- redControls[[index]]
	dyebiasMatrix <- log2(greenControls.current / redControls.current)
	dyebias <- apply(dyebiasMatrix, 2, median)

	# Bisulfite conversion extraction for probe type II:
	index <- match("BISULFITE CONVERSION II",controlNames)
	redControls.current <- redControls[[ index ]]
	bisulfite2 <- colMeans(redControls.current)
	
	# Bisulfite conversion extraction for probe type I:
	index <- match("BISULFITE CONVERSION I",controlNames)
	redControls.current <- redControls[[ index ]][7:9,]
	greenControls.current <- redControls[[ index ]][1:3,]
	bisulfite1 <- colMeans( redControls.current + greenControls.current)
	
	# Staining
	index <- match("STAINING",controlNames)
	sg <- greenControls[[ index ]][3, ] 
	sr <- redControls[[ index ]][1, ] 

	# Extension
	index <- match("EXTENSION",controlNames)
	redControls.current <- redControls[[index]]
	greenControls.current <- greenControls[[index]]
	extr <- redControls.current[1:2,]
	extg <- greenControls.current[3:4,]

	# Hybridization should be monitored only in the green channel
	index <- match("HYBRIDIZATION",controlNames)
	h1 <- greenControls[[index]][1, ]
	h2 <- greenControls[[index]][2, ]
	h3 <- greenControls[[index]][3, ]

	# Target removal should be low compared to hybridization probes
	index <- match("TARGET REMOVAL",controlNames)
	tar <- greenControls[[index]]
	
	# Non-polymorphic probes
	index <- match("NON-POLYMORPHIC",controlNames)
	npr <- redControls[[index]][1:2, ]
	npg <- greenControls[[index]][3:4, ]
	
	# Specificity II
	index <- match("SPECIFICITY II",controlNames)
	greenControls.current <- greenControls[[index]]
	redControls.current <- redControls[[index]]
	spec2g <- colMeans(greenControls.current)
	spec2r <- colMeans(redControls.current)
	spec2ratio <- spec2g / spec2r
	spec2g <- greenControls.current
	spec2r <- redControls.current
	
	# Specificity I
	index <- match("SPECIFICITY I", controlNames)
	greenControls.current <- greenControls[[index]][1:3,]
	redControls.current <- redControls[[index]][7:9,]
	spec1g <- greenControls.current
	spec1r <- redControls.current
	greenControls.current <- greenControls[[index]][1:3,]
	redControls.current <- redControls[[index]][1:3,]
	ratio1 <- colMeans(redControls.current) / colMeans(greenControls.current)
	greenControls.current <- greenControls[[index]][7:9,]
	redControls.current <- redControls[[index]][7:9,]
	ratio2 <- colMeans(greenControls.current) / colMeans(redControls.current)
    spec1ratio <- (ratio1 + ratio2) / 2
	
	# Normalization probes:
	index <- match(c("NORM_A"), controlNames)
	normA <- colMeans(redControls[[index]])
	index <- match(c("NORM_T"), controlNames)
	normT <- colMeans(redControls[[index]])
	index <- match(c("NORM_C"), controlNames)
	normC <- colMeans(redControls[[index]])
	index <- match(c("NORM_G"), controlNames)
	normG <- colMeans(redControls[[index]])
    dyebias2 <- (normC + normG)/(normA+normT)
	
  	model.matrix <- cbind(bisulfite1, bisulfite2, t(extg), t(extr), h1, h2,h3, sg, sr, t(npg),
			t(npr), t(tar), t(spec1g), t(spec1r), t(spec2g), t(spec2r),ratio1,
   			spec1ratio, spec2ratio, ratio2, normA, normC, normT, normG, dyebias2)
   
    # Imputation
    for (colindex in 1:ncol(model.matrix)) {
    	column <- model.matrix[,colindex]
     	column[is.na(column)]<- mean(column, na.rm=T)
     	model.matrix[ ,colindex] <- column
      }          
           
    # Scaling   
    model.matrix <- scale(model.matrix)               
    return(model.matrix)
}
#################################################################


### Return the corrected quantile distributions 
#################################################################
returnFit <- function(model.matrix, quantiles, oobQuantiles, nPCs, dyebias = NULL) {
	library(fda)
	library(refund)
	newQuantiles<- quantiles
	model.matrix <- prcomp(model.matrix)$x[,1:nPCs]
	#model.matrix <- princomp(model.matrix)$scores[,1:nPCs]
	
	### If correcting for dyebias (e.g. when normalizing directly on the beta-values)
	if (!is.null(dyebias)){
		model.matrix <- cbind(model.matrix,dyebias)
	}	
	
	### The model is: Y(s) = u(s)  + UD(s)  + \int_t OA(s,t)dt + error(s)
	meanFunction <- apply(quantiles,1,mean)
	res <- quantiles - meanFunction
	
	### Estimation of UD(s)
	for (i in 1:nrow(quantiles)){
		model <- lm(res[i,] ~ model.matrix-1)
	    res[i,] <- res[i,] - model$fitted.values
	} 

	### Function-on-function regression on the residuals (Estimation of \int_t OA(s,t)dt )
	oobQuantiles <- log(oobQuantiles +1)
	oobQuantiles <- t(scale(t(oobQuantiles),scale=FALSE))
	coeff <- matrix(0,nrow(quantiles),nrow(oobQuantiles))
	for (r in 1:nrow(quantiles)){
		for (s in 1:nrow(oobQuantiles)){
			model <- lm(res[r,]~oobQuantiles[s,]-1)
			coeff[r,s] <- summary(model)$coeff[1,1]
		}
	}
	res <- res - coeff%*%oobQuantiles/nrow(oobQuantiles)
	newQuantiles <- meanFunction + res
	
	newQuantiles[newQuantiles < 0] = 0
	return(newQuantiles)
}
#################################################################

### Assumes that the predicted sex is 1 and 2. The sex prediction function used by default respects this
### Return the corrected quantile distributions for the X-chromosome intensities
#################################################################
returnFitX <- function(model.matrix, quantiles, oobQuantiles, nPCs, dyebias = NULL) {
	
	quantiles1 <- quantiles[, sex == 1]
	model.matrix1 <- model.matrix[sex == 1, ]
	dyebias1 <- dyebias[sex == 1]
	oobQuantiles1 <- oobQuantiles[,sex == 1]
	
	newQuantiles1 <- returnFit(model.matrix = model.matrix1, quantiles = quantiles1, 
			oobQuantiles = oobQuantiles1, nPCs = nPCs, dyebias = dyebias1)
	
	quantiles2 <- quantiles[, sex == 2]
	model.matrix2 <- model.matrix[sex == 2, ]
	dyebias2 <- dyebias[sex == 2]
	oobQuantiles2 <- oobQuantiles[,sex == 2]
	
	newQuantiles2 <- returnFit(model.matrix = model.matrix2, quantiles = quantiles2,
			oobQuantiles = oobQuantiles2, nPCs = nPCs,dyebias = dyebias2)
	
	newQuantiles <- quantiles
	newQuantiles[, sex==1] <- newQuantiles1
	newQuantiles[, sex==2] <- newQuantiles2
	
	return(newQuantiles)
}
#################################################################

### Normalize a matrix of intensities with 
### the corrected quantile distributions
#################################################################
normalizeByType <- function(intMatrix, newQuantiles) {
	normMatrix <- matrix(NA, nrow(intMatrix), ncol(intMatrix))
	library(preprocessCore) 
	n <- 500
	
	normMatrix <- sapply(1:ncol(intMatrix), function(i) {
				crtColumn <- intMatrix[ , i]
				crtColumn.reduced <- crtColumn[!is.na(crtColumn)]
				### Generation of the corrected intensities:
				target <- sapply(1:499, function(j) {
							start <- newQuantiles[j,i]
							end  <- newQuantiles[j+1,i]
							sequence <- seq(start, end,( end-start)/n)[-n]
							return(sequence)
				})

				target <- as.vector(target)
				result <- normalize.quantiles.use.target(matrix(crtColumn.reduced), target)
				print(paste("Sample ", i ," is normalized"))
				return(result)
			})
	return(normMatrix)
}
#################################################################

### Main function call for normalization
### Possibility to protect phenotype data from the normalization
#################################################################
normalizeFunNorm450kCN <- function(cnMatrix, extractedData, nPCs = 4, predictedSex) {
	
	oobQuantiles <- extractedData$oob$greenOOB
	
	probesI <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I")
	probesII <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "II")
    
    ### Chr probes:
	locations <- getLocations(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
	autosomal <- names(locations[seqnames(locations) %in% paste0("chr", 1:22)])
	chrY <- names(locations[seqnames(locations)=="chrY"])
	chrX <- names(locations[seqnames(locations)=="chrX"])
    
    probesIGrn <- intersect( probesI$Name[probesI$Color=="Grn"], autosomal )
	probesIRed <- intersect( probesI$Name[probesI$Color=="Red"], autosomal )
	probesII <- intersect( probesII$Name, autosomal )
    
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
	types <- c("II", "IRed", "IGreen","X","Y")
	indList <- list(II, IRed, IGreen, X, Y)

	cnList <-  list(extractedData$cnQuantiles$II, extractedData$cnQuantiles$IRed, 
			extractedData$cnQuantiles$IGrn, extractedData$cnQuantiles$X,extractedData$cnQuantiles$Y )
	                                            
	model.matrix <- buildControlMatrix450k(extractedData)
	dyebias <- buildDyebias(extractedData)
	dbList <- list(NULL,NULL,dyebias)

	messages <- c("type II","type I Green","type I Red")
	print("Normalization of the autosomal probes...")
	for (i in 1:3) {
		print(paste0("Normalization of the ",messages[i]," probes..."))
		print("Generating the adjusted quantile distributions...")
					 
		newQuantiles <- returnFit(model.matrix = model.matrix, quantiles = cnList[[i]], 
				oobQuantiles = oobQuantiles, nPCs = nPCs, dyebias = dbList[[i]])
				                                                
		print("Normalizing subset of probes...")                                          
		cnMatrix[indList[[i]] , ] <- normalizeByType(cnMatrix[indList[[i] ] ,], newQuantiles)
	}
			
	if (length(X)!=0) {
		print("Normalization of the X-chromosome...")
		newQuantiles <- returnFitX(model.matrix = model.matrix, quantiles = cnList[[4]],
				oobQuantiles = oobQuantiles, nPCs=nPCs, dyebias = dyebias)
					                                             
		cnMatrix[X,] <- normalizeByType(cnMatrix[X,], newQuantiles)
	}
	
	if (length(Y)!=0) {
		print("Quantile normalization of the Y-chromosome...")
		cnMatrix[Y, predictedSex == 1] <- preprocessCore::normalize.quantiles(cnMatrix[Y, predictedSex == 1])
		cnMatrix[Y, predictedSex == 2] <- preprocessCore::normalize.quantiles(cnMatrix[Y, predictedSex == 2])
	}

	print("Normalization done.")
	
	return(cnMatrix)

}
#################################################################
