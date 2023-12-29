#' smoothclust
#' 
#' Function to run 'smoothclust' spatial clustering algorithm.
#' 
#' Function to run 'smoothclust' spatial clustering algorithm.
#' 
#' 
#' @param input Input data, assumed to be provided as \code{SpatialExperiment}
#'   object containing matrix of spatial coordinates in \code{spatialCoords}
#'   slot and log-transformed normalized counts in assay named \code{logcounts}.
#' 
#' @param method Method used for smoothing. The \code{uniform} method calculates
#'   average logcounts (unweighted) across all measurement locations within a
#'   circle with radius \code{bandwidth} at each measurement location, which
#'   smooths out spatial variability as well as sparsity due to sampling
#'   variability. The \code{kernel} method calculates a weighted average using a
#'   truncated exponential kernel applied to Euclidean distances with a length
#'   scale parameter equal to \code{bandwidth}, which provides a more
#'   sophisticated approach to smoothing out spatial variability but may be
#'   affected by sparsity due to sampling variability, especially at the index
#'   point.
#' 
#' @param bandwidth Parameter defining the bandwidth for smoothing, expressed as
#'   the proportion of the width or height (whichever is greater) of the tissue
#'   area. Smoothing is performed across neighboring values of logcounts at each
#'   point. For \code{method = "uniform"}, the bandwidth represents the radius
#'   of a circle, and unweighted average logcounts are calculated across points
#'   within this circle. For \code{method = "kernel"}, the averaging is weighted
#'   by distances scaled using a truncated exponential kernel applied to
#'   Euclidean distances. For example, a bandwidth of 0.05 will smooth logcounts
#'   across neighbors weighted by distances scaled using a truncated exponential
#'   kernel with length scale equal to 5% of the width or height (whichever is
#'   greater) of the tissue area. Weights for \code{method = "kernel"} are
#'   truncated for computational efficiency.
#' 
#' @param truncate Truncation threshold parameter if \code{method = "kernel"}.
#'   Kernel weights below this value are set to zero for computational
#'   efficiency. Ignored if \code{method == "uniform"}.
#' 
#' 
#' @return Returns the \code{SpatialExperiment} object with a new assay named
#'   \code{logcounts_smooth} containing spatially smoothed logcounts values.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assays 'assays<-'
#' @importFrom spdep dnearneigh nbdists
#' @importFrom methods is as
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
#' library(scran)
#' 
#' # download data object
#' spe <- Visium_humanDLPFC()
#' spe <- spe[, colData(spe)$in_tissue == 1]
#' 
#' # calculate highly variable genes (HVGs)
#' spe <- logNormCounts(spe)
#' is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
#' spe <- spe[!is_mito, ]
#' # keep full object for plotting
#' spe_full <- spe
#' 
#' dec <- modelGeneVar(spe)
#' top_hvgs <- getTopHVGs(dec, prop = 0.1)
#' spe <- spe[top_hvgs, ]
#' dim(spe)
#' 
#' # run smoothclust
#' spe <- smoothclust(spe)
#' 
smoothclust <- function(input, method = c("uniform", "kernel"), 
                        bandwidth = 0.05, truncate = 0.05) {
  
  method <- match.arg(method, c("uniform", "kernel"))
  
  spe <- input
  
  stopifnot(is(spe, "SpatialExperiment"))
  
  spatialcoords <- spatialCoords(spe)
  
  # convert bandwidth to same units as distances
  range_x <- abs(diff(range(spatialcoords[, 1])))
  range_y <- abs(diff(range(spatialcoords[, 2])))
  range_max <- max(range_x, range_y)
  bandwidth_scaled <- bandwidth * range_max
  
  if (method == "uniform") {
    # calculate neighbors (note self is excluded)
    neigh <- dnearneigh(spatialcoords, d1 = 0, d2 = bandwidth_scaled)
  }
  
  if (method == "kernel") {
    # calculate neighbors (note self is excluded)
    neigh <- dnearneigh(spatialcoords, d1 = 0, d2 = Inf)
    # calculate distances
    dists <- nbdists(neigh, coords = spatialcoords)
  }
  
  # put back self within set of neighbors for each point
  stopifnot(length(neigh) == ncol(spe))
  # include index of self point
  neigh <- mapply(c, as.list(seq_len(ncol(spe))), neigh, SIMPLIFY = FALSE)
  
  if (method == "kernel") {
    stopifnot(length(dists) == ncol(spe))
    # include distance (zero) to self point
    dists <- mapply(c, 0, dists, SIMPLIFY = FALSE)
  }
  
  # calculate weights for kernel method
  if (method == "kernel") {
    # calculate exponential kernel weights
    exp_kernel <- function(d) {exp(-d / bandwidth_scaled)}  ## d = Euclidean distance
    weights <- lapply(dists, exp_kernel)
    
    # truncate kernel weights below threshold
    keep <- lapply(weights, function(w) {w >= truncate})
    
    stopifnot(length(weights) == length(keep))
    stopifnot(length(neigh) == length(keep))
    
    weights <- mapply(function(w, k) {w[k]}, weights, keep)
    neigh <- mapply(function(n, k) {n[k]}, neigh, keep)
  }
  
  # calculate smoothed logcounts
  # note: using dense matrix to ensure zeros are included in averaging
  lc <- as.matrix(logcounts(spe))
  lc_smooth <- matrix(as.numeric(NA), nrow = nrow(spe), ncol = ncol(spe))
  
  pb <- txtProgressBar(0, ncol(lc_smooth), style = 3)
  
  if (method == "uniform") {
    for (i in seq_len(ncol(lc_smooth))) {
      setTxtProgressBar(pb, i)
      # extract values
      lc_sub <- lc[, neigh[[i]]]
      # calculate average
      lc_smooth[, i] <- rowMeans(lc_sub)
    }
  }
  
  if (method == "kernel") {
    for (i in seq_len(ncol(lc_smooth))) {
      setTxtProgressBar(pb, i)
      # extract values and weights
      lc_sub <- lc[, neigh[[i]]]
      weights_rep <- t(replicate(nrow(lc_sub), weights[[i]]))
      stopifnot(all(dim(lc_sub) == dim(weights_rep)))
      # calculate weighted average
      vals <- rowSums(lc_sub * weights_rep) / sum(weights[[i]])
      lc_smooth[, i] <- vals
    }
  }
  
  close(pb)
  
  stopifnot(nrow(lc_smooth) == nrow(spe))
  stopifnot(ncol(lc_smooth) == ncol(spe))
  rownames(lc_smooth) <- rownames(spe)
  colnames(lc_smooth) <- colnames(spe)
  
  # store in object
  assays(spe)[["logcounts_smooth"]] <- lc_smooth
  
  spe
}
