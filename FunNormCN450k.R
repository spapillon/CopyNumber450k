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
buildDyebias <- function(extractedData){
	 greenControls <- extractedData$greenControls
     redControls <- extractedData$redControls
	 controlNames <- names(greenControls)           
	 
	index <-    match(c("NORM_A"), controlNames)
	normA <-	apply(redControls[[index]], 2, mean)
	index <-    match(c("NORM_T"), controlNames)
	normT <-   apply(redControls[[index]], 2, mean)
	index <-    match(c("NORM_C"), controlNames)
	normC <-apply(greenControls[[index]], 2, mean)
	index <-    match(c("NORM_G"), controlNames)
	normG <- 	apply(greenControls[[index]], 2, mean)
    
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
	 	
	# Dye bias:               (only for Type II probes)
	index     <-match("NEGATIVE",controlNames)
	greenControls.current   <- greenControls[[index]]
	redControls.current     <- redControls[[index]]
	dyebiasMatrix <- log2(greenControls.current / redControls.current)
	dyebias <- apply(dyebiasMatrix, 2, median)

	# Bisulfite conversion extraction for probe type II:
	index <-    match("BISULFITE CONVERSION II",controlNames)
	redControls.current     <- redControls[[ index ]]
	bisulfite2 <- apply(redControls.current, 2, mean)
	
	# Bisulfite conversion extraction for probe type I:
	index <-    match("BISULFITE CONVERSION I",controlNames)
	redControls.current     <- redControls[[ index ]][7:9,]
	greenControls.current     <- redControls[[ index ]][1:3,]
	bisulfite1 <- apply( redControls.current + greenControls.current, 2 ,mean)
               	
	# Staining
	index <-    match("STAINING",controlNames)
	sg  <- greenControls[[ index ]][3, ] 
	sr   <- redControls[[ index ]][1, ] 

	# Extension
	index <-    match("EXTENSION",controlNames)
	redControls.current     <- redControls[[index]]
	greenControls.current     <- greenControls[[index]]
	extr <- redControls.current[1:2,]
	extg <- greenControls.current[3:4,]

	# Hybridization should be monitored only in the green channel
	index <-    match("HYBRIDIZATION",controlNames)
	h1      <-    greenControls[[index]][1, ]
	h2     <-     greenControls[[index]][2, ]
	h3     <-     greenControls[[index]][3, ]

	# Target removal should be low compared to hybridization probes
	index <-    match("TARGET REMOVAL",controlNames)
	tar    <- greenControls[[index]]
	
	# Non-polymorphic probes
	index <-    match("NON-POLYMORPHIC",controlNames)
	npr <- redControls[[index]][1:2, ]
	npg <- greenControls[[index]][3:4, ]
	
	# Specificity II
	index <-    match("SPECIFICITY II",controlNames)
	greenControls.current     <- greenControls[[index]]
	redControls.current     <- redControls[[index]]
	spec2g      <- apply(greenControls.current, 2, mean)
	spec2r       <- apply(redControls.current, 2, mean)
	spec2ratio <- spec2g / spec2r
	spec2g      <- greenControls.current
	spec2r       <- redControls.current
	
	# Specificity I
	index <- match("SPECIFICITY I", controlNames)
	greenControls.current     <- greenControls[[index]][1:3,]
	redControls.current         <- redControls[[index]][7:9,]
	spec1g  <- greenControls.current
	spec1r   <- redControls.current
	greenControls.current     <- greenControls[[index]][1:3,]
	redControls.current         <- redControls[[index]][1:3,]
	ratio1 <- apply(redControls.current, 2, mean) /
	                   apply(greenControls.current, 2, mean)
	greenControls.current     <- greenControls[[index]][7:9,]
	redControls.current         <- redControls[[index]][7:9,]
	ratio2 <- apply(greenControls.current, 2, mean) / 
	                  apply(redControls.current, 2, mean)
    spec1ratio <- (ratio1 + ratio2) / 2
	
	# Normalization probes:
	index <-    match(c("NORM_A"), controlNames)
	normA <-	apply(redControls[[index]], 2, mean)
	index <-    match(c("NORM_T"), controlNames)
	normT <-   apply(redControls[[index]], 2, mean)
	index <-    match(c("NORM_C"), controlNames)
	normC <-apply(greenControls[[index]], 2, mean)
	index <-    match(c("NORM_G"), controlNames)
	normG <- 	apply(greenControls[[index]], 2, mean)
    
    dyebias2 <- (normC + normG)/(normA+normT)
	
  	 model.matrix <- cbind(
   		    bisulfite1, bisulfite2, t(extg), t(extr), h1, h2,h3, sg, sr, t(npg),
   			t(npr), t(tar), t(spec1g), t(spec1r), t(spec2g), t(spec2r),ratio1,
   			spec1ratio, spec2ratio, ratio2, normA, normC, normT, normG, dyebias2)
   
      # Imputation
      for (colindex in 1:ncol(model.matrix)){
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
## Assuming the phenoMatrix is a model.matrix
returnFit <- function(model.matrix, quantiles, oobQuantiles, nPCs, method, phenoMatrix = NULL, dyebias = NULL){
	library(fda)
	library(refund)
	newQuantiles<- quantiles
	model.matrix <- prcomp(model.matrix)$x[,1:nPCs]
	#model.matrix <- princomp(model.matrix)$scores[,1:nPCs]
	
	### If correcting for dyebias (e.g. when normalizing directly on the beta-values)
	if (!is.null(dyebias)){
		model.matrix <- cbind(model.matrix,dyebias)
	}
	
	# Pointwise estimation (using a separate linear 
	# regression model at each point)
	if (method == "pointwise"){
		
		
	### The model is: Y(s) = u(s) + XB(s) + UD(s)  + \int_t OA(s,t)dt + error(s)
		
	meanFunction <- apply(quantiles,1,mean)
	res <- quantiles - meanFunction
	
	### If a phenotype matrix X is provided (XB(s) estimation)
	if (!is.null(phenoMatrix)){
		phenoRes <- res
		for (i in 1:nrow(quantiles)){
			model <- lm(res[i,]~phenoMatrix-1)
			phenoRes[i,] <- model$fitted.values
		}
	res <- res - phenoRes	
	}
	
	### Estimation of UD(s)
	for (i in 1:nrow(quantiles)){
		model     <- lm(res[i,] ~ model.matrix-1)
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
	
	if (!is.null(phenoMatrix)){
		newQuantiles <- meanFunction + phenoRes + res
	} else {
		newQuantiles <- meanFunction + res
	}
	
	
    # Using function-on-scalar regression (can be 
    # very slow and memory consuming)
	} else if (method == "fosr") {
        breaks <- c(seq(0, 0.05, 0.01), 
                              seq(0.10, 0.90, 0.05), 
                                  seq(0.95, 1, 0.01))         
        basis <- create.bspline.basis(
                               rangeval = c(0,1), 
                                    norder = 10, 
                                        breaks = breaks)
                                        
        smooth.quantiles <- smooth.basis(
                                   seq(0, 1, 1/499), 
                                       quantiles, basis)$fd    
        fit <- fosr(
                 fdobj = smooth.quantiles,
                     X = cbind(1, model.matrix))    
        mean  <- eval.fd(seq(0, 1, 1/499), fit$fd[1])
        residuals <- eval.fd(seq(0, 1, 1/499), fit$resid)
        newQuantiles <- as.vector(mean) + residuals
        
        
    # Using 2-step function-on-scalar 
    # regression (for fast computation)
	} else {
	    fit <- fosr2s(
	           t(quantiles[-c(1, 500), ]), 
	               X = model.matrix, 
	                   argvals = seq(0, 1, 1/499)[ -c(1, 500) ],
	                        basistype = "bspline",
	                             nbasis = 30)
		meanQuantile <- fit$est.func[,1]
		residuals <- quantiles[-c(1, 500),] - t(fit$yhat)
		newQuantiles[-c(1, 500), ] <- as.vector(meanQuantile) + residuals
	}
	
	
	newQuantiles[newQuantiles < 0] = 0
	return(newQuantiles)
}
#################################################################













### Assumes that the predicted sex is 1 and 2. The sex prediction function used by default respects this
### Return the corrected quantile distributions for the X-chromosome intensities
#################################################################
returnFitX <- function(model.matrix, quantiles, oobQuantiles, nPCs, method, sex, phenoMatrix = NULL, dyebias = NULL){
	
	quantiles1 <- quantiles[, sex == 1]
	model.matrix1 <- model.matrix[sex == 1, ]
	phenoMatrix1 <- phenoMatrix[sex == 1,]
	dyebias1 <- dyebias[sex == 1]
	oobQuantiles1 <- oobQuantiles[,sex==1]
	
	newQuantiles1 <- returnFit(model.matrix = model.matrix1,
	                                   quantiles = quantiles1, 
	                                     oobQuantiles = oobQuantiles1, 
	                                       nPCs = nPCs,
	                                        method = method,
	                                          phenoMatrix = phenoMatrix1, 
	                                          dyebias = dyebias1)
	
	quantiles2 <- quantiles[, sex == 2]
	model.matrix2 <- model.matrix[sex == 2, ]
	phenoMatrix2 <- phenoMatrix[sex == 2,]
	dyebias2 <- dyebias[sex == 2]
	oobQuantiles2 <- oobQuantiles[,sex==2]
	
	newQuantiles2 <- returnFit(model.matrix = model.matrix2,
	                                   quantiles = quantiles2, 
	                                     oobQuantiles = oobQuantiles2, 
	                                       nPCs = nPCs,
	                                        method = method,
	                                          phenoMatrix = phenoMatrix2, 
	                                          dyebias = dyebias2)
	
	newQuantiles <- quantiles
	newQuantiles[, sex==1] <- newQuantiles1
	newQuantiles[, sex==2] <- newQuantiles2
	
	return(newQuantiles)
}
#################################################################

















### Normalize a matrix of intensities with 
### the corrected quantile distributions
#################################################################
normalizeByType <- function(intMatrix, newQuantiles){
	normMatrix <- matrix(NA, nrow(intMatrix), ncol(intMatrix))
	library(preprocessCore) 
	x <- seq(0,1,1/499)
	n <- 1000
	for (si in 1:ncol(intMatrix)){
		crtColumn <- intMatrix[ , si]
		crtColumn.reduced <- crtColumn[!is.na(crtColumn)]

		### Simulation of the normalized distribution: 
		
		target <- c()
		for (cell in 1:499){
			start <- newQuantiles[cell,si]
			end  <- newQuantiles[cell+1,si]
			sequence <- seq(start, end,( end-start)/n)[-n]
			target <- c(target,sequence)
		 }
			 
				
		### Quantile normalization with the simulated normalized distribution:
		normMatrix[,si][!is.na(crtColumn)] <- 
		               normalize.quantiles.use.target(
		                     matrix(crtColumn.reduced), 
		                         target)
		print(paste("Sample ", si ," is normalized"))
	}
	return(normMatrix)
}
#################################################################















### Main function call for normalization
### Possibility to protect phenotype data from the normalization
#################################################################
normalizeFunNorm450kCN <- function(cnMatrix, 
                                                                   extractedData, 
                                                                     nPCs = 4, 
                                                                         method = "pointwise",
                                                                             predictedSex, 
                                                                                 phenoMatrix = NULL ){
	
	
	oobQuantiles <- extractedData$oob$greenOOB
	
	probesI <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "I")
	        
    probesII <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "II")
    
    ### Chr probes:
	locations <- getLocations(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
	autosomal <- names(locations[seqnames(locations) %in% paste0("chr", 1:22)])
	chrY <- names(locations[seqnames(locations)=="chrY"])
	chrX <- names(locations[seqnames(locations)=="chrX"])
    
    probesIGrn <- intersect( probesI$Name[probesI$Color=="Grn"], autosomal )
	probesIRed   <- intersect( probesI$Name[probesI$Color=="Red"], autosomal )
	probesII     <- intersect( probesII$Name, autosomal )
    
    uProbeNames <- rownames(cnMatrix)
    uProbesIGrn <- intersect(uProbeNames, probesIGrn)
	uProbesIRed <- intersect(uProbeNames, probesIRed)
	uProbesII   <- intersect(uProbeNames, probesII)
    
	II <- match(uProbesII, uProbeNames)
	IRed <- match(uProbesIRed, uProbeNames)
	IGreen <- match(uProbesIGrn, uProbeNames)
	X <- match(chrX, uProbeNames)
	Y <- match(chrY, uProbeNames)
	types <- c("II", "IRed", "IGreen","X","Y")
	indList <- list(II, IRed, IGreen, X, Y)

   
	
	cnList <-  list(extractedData$cnQuantiles$II,
	                                extractedData$cnQuantiles$IRed,
	                                    extractedData$cnQuantiles$IGrn,
	                                        extractedData$cnQuantiles$X,
	                                            extractedData$cnQuantiles$Y )
	                                            
	model.matrix <- buildControlMatrix450k(extractedData)
	dyebias <- buildDyebias(extractedData)
	dbList <- list(NULL,NULL,dyebias)

	
		
				messages <- c("type II","type I Green","type I Red")
		        print("Normalization of the autosomal probes...")
				for (i in 1:3){
					 print(paste0("Normalization of the ",messages[i]," probes..."))
					 print("Generating the adjusted quantile distributions...")
					 
				     newQuantiles <- returnFit(model.matrix = model.matrix,
				                                                quantiles = cnList[[i]], 
				                                                oobQuantiles = oobQuantiles, 
				                                                nPCs = nPCs,
				                                                method = method,
				                                                phenoMatrix = phenoMatrix,
				                                                dyebias = dbList[[i]])
				                                                
				     print("Normalizing subset of probes...")                                          
				     cnMatrix[indList[[i]] , ] <- normalizeByType(cnMatrix[indList[[i] ] ,], newQuantiles)
				}
			
				print("Normalization of the X-chromosome...")
				newQuantiles <- returnFitX(model.matrix = model.matrix,
				                                             quantiles = cnList[[4]],
				                                             oobQuantiles = oobQuantiles,
				                                             nPCs=nPCs, 
				                                             method = method, 
				                                             sex = predictedSex,
				                                             phenoMatrix = phenoMatrix,
				                                             dyebias = dyebias)
				                                             
				cnMatrix[X,] <- normalizeByType(cnMatrix[X,], newQuantiles)
				
				print("Quantile normalization of the Y-chromosome...")
				cnMatrix[Y, predictedSex == 1] <- preprocessCore::normalize.quantiles(cnMatrix[Y, predictedSex == 1])
				cnMatrix[Y, predictedSex == 2] <- preprocessCore::normalize.quantiles(cnMatrix[Y, predictedSex == 2])
	
	 			print("Normalization done.")
	 			
	 			return(cnMatrix)
	 			
	 					

	}
#################################################################









