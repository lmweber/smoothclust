% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothclust.R
\name{smoothclust}
\alias{smoothclust}
\title{smoothclust}
\usage{
smoothclust(
  input,
  method = c("uniform", "kernel"),
  bandwidth = 0.05,
  truncate = 0.05
)
}
\arguments{
\item{input}{Input data, assumed to be provided as \code{SpatialExperiment}
object containing matrix of spatial coordinates in \code{spatialCoords}
slot and log-transformed normalized counts in assay named \code{logcounts}.}

\item{method}{Method used for smoothing. The \code{uniform} method calculates
average logcounts (unweighted) across all measurement locations within a
circular window with radius \code{bandwidth} at each measurement location,
which smooths out spatial variability as well as sparsity due to sampling
variability. The \code{kernel} method calculates a weighted average using a
truncated exponential kernel applied to Euclidean distances with a length
scale parameter equal to \code{bandwidth}, which provides a more
sophisticated approach to smoothing out spatial variability but may be
affected by sparsity due to sampling variability, especially sparsity at
the index point.}

\item{bandwidth}{Bandwidth parameter for smoothing, expressed as proportion
of width or height (whichever is greater) of tissue area. Smoothing is
performed across neighboring values of logcounts at each point. For
\code{method = "uniform"}, the bandwidth represents the radius of a circle,
and unweighted average logcounts are calculated across points within this
circle. For \code{method = "kernel"}, the averaging is weighted by
distances scaled using a truncated exponential kernel applied to Euclidean
distances. For example, a bandwidth of 0.05 will smooth logcounts across
neighbors weighted by distances scaled using a truncated exponential kernel
with length scale equal to 5% of the width or height (whichever is greater)
of the tissue area. Weights for \code{method = "kernel"} are truncated at
small values for computational efficiency.}

\item{truncate}{Truncation threshold parameter if \code{method = "kernel"}.
Kernel weights below this value are set to zero for computational
efficiency. Ignored if \code{method == "uniform"}.}
}
\value{
Returns the \code{SpatialExperiment} object with a new assay named
  \code{logcounts_smooth} containing spatially smoothed logcounts values,
  which can be used as the input for further downstream analyses such as
  clustering.
}
\description{
Spatial clustering algorithm for spatial transcriptomics data.
}
\details{
Spatial clustering algorithm for spatial transcriptomics data based on the
principle of smoothing expression measurements across neighboring spatial
locations. The algorithm can be used to define spatial domains consisting of
a single cell type or a consistent mixture of cell types, with clearly
defined spatial boundaries.
}
\examples{
library(SpatialExperiment)
library(STexampleData)
library(scran)

# download data object
spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

# calculate highly variable genes (HVGs)
spe <- logNormCounts(spe)
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
spe <- spe[!is_mito, ]
# keep full object for plotting
spe_full <- spe

dec <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
spe <- spe[top_hvgs, ]
dim(spe)

# run smoothclust
spe <- smoothclust(spe)

# check
assayNames(spe)

}
