################################################################################  
# Using minfiData
################################################################################  

library(CopyNumber450k)
load("~/git/CopyNumber450k/data/control_RGset.RData")

baseDir <- system.file("extdata", package = "minfiData")
targets <- read.450k.sheet(baseDir)
case_RGset <- read.450k.exp(base = baseDir, targets = targets)
control_RGset <- updateObject(control_RGset)
RGset <- combine(control_RGset, case_RGset)
mcds <- CNV450kSet(RGset)
mcds <- dropSNPprobes(mcds)
mcds.q <- normalize(mcds, "quantile")
mcds.q.seg <- segmentize(mcds.q)

par(mfrow=c(1,3))
plotDensity(mcds)
plotDensity(mcds.f)
plotDensity(mcds.q)