\name{CNV450kSet-class}
\Rdversion{1.0}
\docType{class}
\alias{CopyNumber450k}
\alias{CNV450kSet-class}
\alias{CNV450kSet}
% Methods
\alias{getManifest,CNV450kSet-method}
\alias{getSegments,CNV450kSet-method}
\alias{getSummary,CNV450kSet-method}
\alias{dropSNPprobes,CNV450kSet-method}
\alias{normalize,CNV450kSet-method}
\alias{segmentize,CNV450kSet-method}
\alias{findCNV,CNV450kSet-method}
\alias{intersectCNV,CNV450kSet-method}
\alias{subgroupDifference,CNV450kSet-method}
\alias{plotSample,CNV450kSet-method}
\alias{plotPCA,CNV450kSet-method}
\alias{plotDensity,CNV450kSet-method}
\alias{computeSignificance,CNV450kSet-method}
\alias{write.csv,CNV450kSet-method}

\alias{getManifest}
\alias{getSegments}
\alias{getSummary}
\alias{dropSNPprobes}
\alias{normalize}
\alias{segmentize}
\alias{findCNV}
\alias{intersectCNV}
\alias{subgroupDifference}
\alias{plotSample}
\alias{plotPCA}
\alias{plotDensity}
\alias{computeSignificance}

\title{CNV450kSet instances}

\description{
  This class holds Illumina 450k methylation microarray data and annotations for CNV calling.
}
\usage{
  ## Constructor
  CNV450kSet(RGChannelSet)
}
\arguments{
  \item{RGChannelSet}{A class defined in the \code{minfi} package that represents raw (unprocessed) data from a two color microarray; specifically an Illumina methylation array.}
}
\details{
  This class inherits from \code{eSet}.  The class is a representation of an intensity matrix that is computed from the \code{MethylSet} class produced by \code{RGChannelSet} in the \code{minfi} package.
}
\section{Constructor}{
Instances are constructed using the \code{CNV450kSet} function with the arguments outlined above.
}
\section{Accessors}{
This class has a few accessors defined on top of those provided in \code{eSet}.  In the following code, \code{object} is a \code{CNV450kSet}.

  \describe{
    \item{\code{getManifest(object)}}{Gets the manifest associated with the object.}
    \item{\code{getSegments(object)}}{Gets list of segments produced by the \code{segmentize} method.}
    \item{\code{getSummary(object)}}{Gets statistics extracted from raw data (\code{RGChannelSet}) and preprocessed methylation data (\code{MethylSet}) that are used internally.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{dropSNPprobes(object, maf_threshold = 0)}}{Removes probes mapping to or targeting known SNPs from consideration.  Returns a new CNV450kSet object in which SNP-containing probes are discarded.}
    \item{\code{normalize(object, type = c("functional", "quantile"))}}{Normalizes the intensity matrix and returns normalized object.  Refer to the vignette for details on the types of normalization.}
    \item{\code{segmentize(object, verbose = TRUE, p.adjust.method = "bonferroni")}}{Arranges probes into segments related by signal intensity and proximity by circular binary segmentation.  Certain class methods require that segments be created.  Refer to the vignette for more details on the algorithm.}
    \item{\code{computeSignificance(object, p.value.threshold = 0.01, num.mark.threshold = 10)}}{Computes the significance for each segment in each sample.  Requires that \code{segmentize} be called on the object beforehand.}
  }
}
\section{Comparison Methods}{
  The following methods require that \code{segmentize} be called on the object before their calls.

  \describe{
    \item{\code{findCNV(object, gene_names, type = c("gain", "loss", "both"))}}{Returns a matrix in which \code{i, j} denotes the presence or absence of a CNV event for gene \code{i} and sample \code{j}. \code{gene_names} is a vector containing the gene symbols of interest.}
    \item{\code{intersectCNV(object, sample_indices, type = c("gain", "loss", "both"))}}{Returns a vector of gene symbols corresponding to gain or loss of genes within a group, sorted by CNV abundance; \code{sample_indices} is a vector containing the indices of the samples belonging to the group.}
    \item{\code{subgroupDifference(object, group1_indices, group2_indices)}}{Returns two vectors (gains and losses) containing Fisher's exact test p-values on gene-based CNV counts between two sample groups (or conditions); \code{group1_indices} and \code{group2_indices} are vectors containing the indices of the samples belonging to each respective group.}
  }
}
\section{Plotting Methods}{
  \describe{
    \item{\code{plotSample(object, index, chr, start, end)}}{Produces a plot of the genomic segments and relative values for sample at \code{index}. \code{chr, start ,end} are optional parameters to be used to zoom in a specific genomic location.}
    \item{\code{plotDensity(object, color.by = c("array.row", "array.col", "sample.group", "slide", "origin"),\\
    			color.function = rainbow, legend.position = "topright")}}{Plots the density distribution of the intensity matrix of the object.}
    \item{\code{plotPCA(object, color.by = c("array.row", "array.col", "sample.group", "slide", "origin"),\\
    			color.function = rainbow, legend.position = "topright")}}{Plots the PCA scatter plot of the intensity matrix of the object.}
    \item{\code{write.csv(object, ...)}}{Writes the segment output for each sample in csv format.}
  }
}
\examples{
library(CopyNumber450kData)
library(minfiData)

# Load the CopyNumber450kData control set
data(RGcontrolSetEx)

# Load example data (n=6) from minfiData
data(RGsetEx)
# In order to reduce example time, let's use only one sample
# and 30 controls (instead of 52). In real life situations, it is advised to
# use all the available controls
RGsetEx <- RGsetEx[, 5]
RGcontrolSetEx <- RGcontrolSetEx[, sample(1:ncol(RGcontrolSetEx), 30)]

# Combine both RGsets in a single RGset
RGset <- combine(RGcontrolSetEx, RGsetEx)

# Create the object
mcds <- CNV450kSet(RGset)

# In order to speed up example computation, we will randomly subset the
# probes used by CopyNumber450k. THIS SHOULD NEVER BE DONE AS IT SERVES
# ONLY FOR SPEEDING UP THE EXAMPLE.
mcds <- mcds[sample(1:nrow(mcds), 10000), ]

# Drop SNP probes
mcds <- dropSNPprobes(mcds, maf_threshold=0.01)

# Normalization
mcds <- normalize(mcds, "quantile")

# Some plots
plotDensity(mcds, main="Density plot of functional normalized data")
plotPCA(mcds, main="PCA plot of functional normalized data")

# Segmentation
mcds <- segmentize(mcds)

# Plotting the results 
plotSample(mcds, 1, main="Genomic view of Sample 1")
plotSample(mcds, 1, chr="chr1", ylim=c(-.25,.25))
        
# Saving the results in csv format
write.csv(mcds, file="segments.csv")
}
\author{Simon Papillon-Cavanagh, Jean-Philippe Fortin, Nicolas De Jay}
\seealso{
  \code{\linkS4class{eSet}} for the basic class structure.  Objects of this class are typically created from a \code{\link{RGChannelSet}} using the constructor.
}
