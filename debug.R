################################################################################ 
library(CopyNumber450k)
load("~/git/CopyNumber450k/data/control_RGset.RData")

#pData(control_RGset)$Sample_Group[1] <- "case"
path <- "~/Documents/iChange/data_ETMR/"
case_RGset <- read.450k.exp(base = path, targets = read.450k.sheet(path))
#pData(case_RGset)$Sample_Group[-c(1,2)]  <- "control"
control_RGset <- updateObject(control_RGset)
RGset <- combine(control_RGset, case_RGset)
mcds <- CNV450kSet(RGset)
mcds <- dropSNPprobes(RGset)
#mcds.f <- normalize(mcds, "functional")
mcds.q <- normalize(mcds, "quantile")
mcds.q.seg <- segmentize(mcds.q)

par(mfrow=c(1,3))
plotDensity(mcds)
plotDensity(mcds.f)
plotDensity(mcds.q)

################################################################################  
