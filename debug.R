################################################################################ 
library(CopyNumber450k)
load("~/git/CopyNumber450k/data/control_RGset.RData")

pData(control_RGset)$Sample_Group[1] <- "case"

control_RGset <- updateObject(control_RGset)

mcds <- CNV450kSet(control_RGset)



pData(mcds)$Sample_Sex <- predictSampleSexes(mcds)

mcds.f <- normalize(mcds, "functional")
mcds.q <- normalize(mcds, "quantile")
par(mfrow=c(1,3))
plotDensity(mcds)
plotDensity(mcds.f)
plotDensity(mcds.q)

################################################################################  
