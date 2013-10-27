\name{CNV450kSet-class}
\Rdversion{1.0}
\docType{class}
\alias{CNV450kSet-Class}
\alias{CNV450kSet}
\alias{getManifest}
\alias{getSegments}
\alias{getSummary}
\alias{normalize}
\alias{segmentize}
\alias{findCNV}
\alias{intersectCNV}
\alias{subgroupDifference}
\title{CNV450kSet instances}
\description{
  This class holds Illumina 450k methylation micro array data and annotations for CNV calling.
}
\usage{
  ## Constructor

  cs <- CNV450kSet(RGChannelSet)
}
\arguments{
  \item{RGChannelSet}{A class defined in the \code{minfi} package that represents raw (unprocessed) data from a two color micro array; specifically an Illumina methylation array.}
}
\details{
  This class inherits from \code{eSet}.  The class is a representation of an intensity matrix that is computed from the \code{MethylSet} class produced by \code{RGChannelSet} in the \code{minfi} package.
}
\section{Constructor}{
Instances are constructed using the \code{CNV450kSet} function with the arguments outlined above.
}
\section{Accessors}{
This class has a few accessors defined on top of those provided in \code{eSet}.  In the following code, \code{object} is a \code{MethylCNVDataSet}.

  \describe{
    \item{\code{getManifest(object)}}{Gets the manifest associated with the object.}
    \item{\code{getSegments(object)}}{Gets list of segments produced by the \code{segmentize} method.}
    \item{\code{getSummary(object)}}{Gets statistics extracted from raw data (\code{RGChannelSet}) and preprocessed methylation data (\code{MethylSet}) that are used internally.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{dropSNPprobes(object, maf_threshold = 0)}}{Removes probes mapping to known SNPs from consideration.  Returns a new CNV450kSet object in which SNP-containing probes are discarded.}
    \item{\code{predictSampleSexes(object)}}{Predicts sex for each sample based on methylation data and returns the sex-sample mapping.}
    \item{\code{normalize(object, type = c("functional", "quantile"))}}{Normalizes the intensity matrix and returns normalized object.  Refer to the vignette for details on the types of normalization.}
    \item{\code{segmentize(object, verbose = TRUE, p.adjust.method = "bonferroni", plotting = FALSE)}}{Arranges probes into segments related by signal intensity and proximity.  Certain class methods require that segments be created.  Can optionally produce a plot ???.  Refer to the vignette for more details on the algorithm.}
  }
}
\section{Comparison Methods}{
  The following methods require that \code{segmentize} be called on the object before their calls.

  \describe{
    \item{\code{findCNV(object, CNVs, type = c("gain", "loss", "both"))}}{Returns ??? corresponding to CNV gain, loss or both for a given set of \code{CNVs}.}
    \item{\code{intersectCNV(object, sample_indices, type = c("gain", "loss", "both"))}}{Returns a list of gene symbols corresponding to gain or loss of CNV within a group; \code{sample_indices} is a vector containing the indices of the samples belonging to the group.}
    \item{\code{subgroupDifference(object, group1_indices, group2_indices)}}{Infers (Differential) CNV gains and losses between two sample groups (or conditions); \code{group1_indices} and \code{group2_indices} are vectors containing the indices of the samples belonging to each respective group.}
  }
}
\section{Plotting Methods}{
  \describe{
    \item{\code{plotSample(object, index, chr, start, end)}}{}
    \item{\code{plotSex(object, color.by = c("array.row", "array.col", "sample.group", "slide", "origin"), color.function = rainbow, legend.position = "topright")}}{}
    \item{\code{plotDensity(object, color.by = c("array.row", "array.col", "sample.group", "slide", "origin"), color.function = rainbow, legend.position = "topright")}}{}
    \item{\code{plotPCA(object, color.by = c("array.row", "array.col", "sample.group", "slide", "origin"), color.function = rainbow, legend.position = "topright")}}{}
    \item{\code{getColoring(object, color.by = c("array.row", "array.col", "sample.group", "slide", "origin"), color.function = rainbow)}}{}
  }
}
\author{Simon Papillon-Cavanagh, Jean-Philippe Fortin, Nicolas De Jay}
\seealso{
  \code{\linkS4class{eSet}} for the basic class structure.  Objects of this class are typically created from a \code{\link{RGChannelSet}} using the constructor.
}