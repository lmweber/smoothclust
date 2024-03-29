#' smoothclust
#' 
#' Method for segmentation of spatial domains and spatially-aware clustering.
#' 
#' Method for segmentation of spatial domains and spatially-aware clustering in
#' spatial transcriptomics data.
#' 
#' Method for segmentation of spatial domains and spatially-aware clustering in
#' spatial transcriptomics data. The method generates spatial domains with
#' smooth boundaries by smoothing gene expression profiles across neighboring
#' spatial locations, followed by unsupervised clustering. Spatial domains
#' consisting of consistent mixtures of cell types may then be further
#' investigated by applying cell type compositional analyses or differential
#' analyses.
#' 
#' 
#' @param input Input data, which can be provided as either a
#'   \code{SpatialExperiment} object or a numeric matrix. If this is a
#'   \code{SpatialExperiment} object, it is assumed to contain either raw
#'   expression counts or logcounts in the \code{assay} slots and spatial
#'   coordinates in the \code{spatialCoords} slot. If this is a numeric matrix,
#'   it is assumed to contain either raw expression counts or logcounts, and
#'   spatial coordinates need to be provided separately with the
#'   \code{spatial_coords} argument.
#' 
#' @param assay_name For a \code{SpatialExperiment} input object, this argument
#'   specifies the name of the \code{assay} containing the expression values to
#'   be smoothed. In most cases, this will be \code{counts}, which contains raw
#'   expression counts. Alternatively, \code{logcounts} may also be used. Note
#'   that if \code{logcounts} are used, the smoothed values represent geometric
#'   averages, which are more difficult to interpret. We recommend using raw
#'   counts if possible. This argument is only used if the input is a
#'   \code{SpatialExperiment} object. Default = \code{counts}.
#' 
#' @param spatial_coords Numeric matrix of spatial coordinates, assumed to
#'   contain x coordinates in first column and y coordinates in second column.
#'   This argument is only used if the input is a numeric matrix.
#' 
#' @param method Method used for smoothing. Options are \code{uniform},
#'   \code{kernel}, and \code{knn}. The \code{uniform} method calculates
#'   unweighted averages across spatial locations within a circular window with
#'   radius \code{bandwidth} at each spatial location, which smooths out spatial
#'   variability as well as sparsity due to sampling variability. The
#'   \code{kernel} method calculates a weighted average using a truncated
#'   exponential kernel applied to Euclidean distances with a length scale
#'   parameter equal to \code{bandwidth}, which provides a more sophisticated
#'   approach to smoothing out spatial variability but may be affected by
#'   sparsity due to sampling variability (especially sparsity at the index
#'   point), and is computationally slower. The \code{knn} method calculates an
#'   unweighted average across the index point and its k nearest neighbors, and
#'   is the fastest method. Default = \code{uniform}.
#' 
#' @param bandwidth Bandwidth parameter for smoothing, expressed as proportion
#'   of width or height (whichever is greater) of tissue area. Only used for
#'   \code{method = "uniform"} or \code{method = "kernel"}. For \code{method =
#'   "uniform"}, the bandwidth represents the radius of a circle, and unweighted
#'   averages are calculated across neighboring points within this circle. For
#'   \code{method = "kernel"}, the averaging is weighted by distances scaled
#'   using a truncated exponential kernel applied to Euclidean distances. For
#'   example, a bandwidth of 0.05 will smooth values across neighbors weighted
#'   by distances scaled using a truncated exponential kernel with length scale
#'   equal to 5% of the width or height (whichever is greater) of the tissue
#'   area. Weights for \code{method = "kernel"} are truncated at small values
#'   for computational efficiency. Default = 0.05.
#' 
#' @param k Number of nearest neighbors parameter for \code{method = "knn"}.
#'   Only used for \code{method == "knn"}. Unweighted averages are calculated
#'   across the index point and its k nearest neighbors. Default = 18 (based on
#'   two layers in honeycomb pattern for 10x Genomics Visium platform).
#' 
#' @param truncate Truncation threshold parameter if \code{method = "kernel"}.
#'   Kernel weights below this value are set to zero for computational
#'   efficiency. Only used for \code{method = "kernel"}. Default = 0.05.
#' 
#' 
#' @return Returns spatially smoothed expression values, which can then be used
#'   as the input for further downstream analyses. Results are returned either
#'   as a \code{SpatialExperiment} object containing a new \code{assay} named
#'   \code{<assay_name>_smooth} (e.g. \code{counts_smooth} or
#'   \code{logcounts_smooth}), or as a numeric matrix, depending on the input
#'   type.
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assays 'assays<-' assayNames
#' @importFrom sparseMatrixStats rowMeans2 rowWeightedMeans
#' @importFrom spdep dnearneigh nbdists knearneigh
#' @importFrom methods is as
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export
#' 
#' @examples
#' library(STexampleData)
#' 
#' # load data
#' spe <- Visium_humanDLPFC()
#' # keep spots over tissue
#' spe <- spe[, colData(spe)$in_tissue == 1]
#' 
#' # run smoothclust
#' # using "knn" method for faster runtime in this example
#' spe <- smoothclust(spe, method = "knn", k = 6)
#' 
#' # see vignette for extended example using default method and including
#' # downstream analysis steps
#' 
smoothclust <- function(input, assay_name = "counts", spatial_coords = NULL, 
                        method = c("uniform", "kernel", "knn"), 
                        bandwidth = 0.05, truncate = 0.05, k = 18) {
  
  method <- match.arg(method, c("uniform", "kernel", "knn"))
  
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
    vals <- assays(spe)[[assay_name]]
    spatial_coords <- spatialCoords(spe)
  } else {
    stopifnot(is.numeric(input))
    vals <- input
    stopifnot(is.numeric(spatial_coords))
  }
  
  if (method %in% c("uniform", "kernel")) {
    # convert bandwidth to same units as distances
    range_x <- abs(diff(range(spatial_coords[, 1])))
    range_y <- abs(diff(range(spatial_coords[, 2])))
    range_max <- max(range_x, range_y)
    bandwidth_scaled <- bandwidth * range_max
  }
  
  if (method == "uniform") {
    # calculate neighbors (note self is excluded)
    neigh <- dnearneigh(spatial_coords, d1 = 0, d2 = bandwidth_scaled)
  }
  
  if (method == "kernel") {
    # calculate neighbors (note self is excluded)
    neigh <- dnearneigh(spatial_coords, d1 = 0, d2 = Inf)
    # calculate distances
    dists <- nbdists(neigh, coords = spatial_coords)
  }
  
  if (method %in% c("uniform", "kernel")) {
    # put back self within set of neighbors for each point
    # note: self point is first element in vector
    stopifnot(length(neigh) == ncol(vals))
    # include index of self point
    neigh <- mapply(c, as.list(seq_len(ncol(vals))), neigh, SIMPLIFY = FALSE)
    if (method == "kernel") {
      stopifnot(length(dists) == ncol(vals))
      # include distance (zero) to self point
      dists <- mapply(c, 0, dists, SIMPLIFY = FALSE)
    }
  }
  
  # calculate weights for kernel method
  if (method == "kernel") {
    # calculate exponential kernel weights
    exp_kernel <- function(d) {exp(-d / bandwidth_scaled)}  ## d = Euclidean distance
    weights <- lapply(dists, exp_kernel)
    
    # truncate kernel weights below threshold
    keep <- lapply(weights, function(w) {w >= truncate})
    
    stopifnot(length(weights) == length(neigh))
    stopifnot(length(weights) == length(dists))
    stopifnot(length(weights) == length(keep))
    stopifnot(all(sapply(weights, length) == length(weights)))
    stopifnot(all(sapply(keep, length) == length(keep)))
    
    # truncate weights and neighbors
    # weights <- mapply(function(w, k) {w[k]}, weights, keep)
    # neigh <- mapply(function(n, k) {n[k]}, neigh, keep)
    
    # truncate weights and fill vector with zeros
    # note rowWeightedMeans() requires full-length weights vectors so cannot subset
    weights <- mapply(function(w, k) {
      w_trunc <- rep(0, length(w))
      w_trunc[k] <- w[k]
      w_trunc
    }, weights, keep, SIMPLIFY = FALSE)
    
    # order weights to match original order of points
    weights <- mapply(function(w, n) {
      w[order(n)]
    }, weights, neigh, SIMPLIFY = FALSE)
  }
  
  if (method == "knn") {
    neigh <- knearneigh(spatial_coords, k = k)$nn
    # include index point
    stopifnot(nrow(neigh) == ncol(vals))
    neigh <- cbind(seq_len(nrow(neigh)), neigh)
  }
  
  # calculate smoothed values
  # note: using sparse matrices
  
  vals_smooth <- matrix(as.numeric(NA), nrow = nrow(vals), ncol = ncol(vals))
  
  # sparse matrix class for sparseMatrixStats
  vals <- as(vals, "CsparseMatrix")
  stopifnot(all(dim(vals) == dim(vals_smooth)))
  
  pb <- txtProgressBar(0, ncol(vals_smooth), style = 3)
  
  if (method == "uniform") {
    for (i in seq_len(ncol(vals_smooth))) {
      setTxtProgressBar(pb, i)
      # calculate average over subset of neighbors
      vals_smooth[, i] <- rowMeans2(vals, cols = neigh[[i]])
    }
  }
  
  if (method == "kernel") {
    for (i in seq_len(ncol(vals_smooth))) {
      setTxtProgressBar(pb, i)
      # calculate weighted average over subset of neighbors
      vals_smooth[, i] <- rowWeightedMeans(vals, w = weights[[i]])
    }
  }
  
  if (method == "knn") {
    stopifnot(nrow(neigh) == ncol(vals_smooth))
    for (i in seq_len(ncol(vals_smooth))) {
      setTxtProgressBar(pb, i)
      # calculate average over subset of neighbors
      vals_smooth[, i] <- rowMeans2(vals, cols = neigh[i, ])
    }
  }
  
  close(pb)
  
  stopifnot(nrow(vals_smooth) == nrow(input))
  stopifnot(ncol(vals_smooth) == ncol(input))
  rownames(vals_smooth) <- rownames(input)
  colnames(vals_smooth) <- colnames(input)
  
  # return results (smoothed values)
  if (is(input, "SpatialExperiment")) {
    assay_name_smooth <- paste0(assay_name, "_smooth")
    assays(spe)[[assay_name_smooth]] <- vals_smooth
    spe
  } else {
    vals_smooth
  }
}
