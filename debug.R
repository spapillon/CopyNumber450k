################################################################################  
# Using minfiData
################################################################################  

library(CopyNumber450k)
library(CopyNumber450kData)
library(minfiData)
data(RGcontrolSetEx)
data(RGsetEx)

RGset <- combine(RGcontrolSetEx, RGsetEx)

mcds <- CNV450kSet(RGset)
mcds <- dropSNPprobes(mcds)
mcds.q <- normalize(mcds, "quantile")
mcds.q <- segmentize(mcds.q)
