#' smoothclust
#' 
#' Spatial clustering algorithm for spatial transcriptomics data.
#' 
#' Spatial clustering algorithm for spatial transcriptomics data based on the
#' principle of smoothing expression measurements across neighboring spatial
#' locations. The algorithm can be used to define spatial domains consisting of
#' a single cell type or a consistent mixture of cell types with smooth spatial
#' boundaries.
#' 
#' 
#' @param input Input data, assumed to be provided as \code{SpatialExperiment}
#'   object containing spatial coordinates in \code{spatialCoords} slot and
#'   expression counts (raw or transformed) in \code{assays} slot.
#' 
#' @param assay_name Name of \code{assay} in input object containing expression
#'   count values to be smoothed. Usually this will be either \code{counts} (raw
#'   counts) or \code{logcounts} (log-transformed normalized counts). We
#'   recommend using raw counts (\code{counts}) if available, in which case the
#'   methodology will apply smoothing to the raw counts and then calculate
#'   logcounts. Alternatively, logcounts (\code{logcounts}) may also be used as
#'   the input, in which case the smoothed values represent geometric averages,
#'   which are more difficult to interpret. Default = \code{counts}.
#' 
#' @param method Method used for smoothing. The \code{uniform} method calculates
#'   unweighted averages across spatial locations within a circular window with
#'   radius \code{bandwidth} at each spatial location, which smooths out spatial
#'   variability as well as sparsity due to sampling variability. The
#'   \code{kernel} method calculates a weighted average using a truncated
#'   exponential kernel applied to Euclidean distances with a length scale
#'   parameter equal to \code{bandwidth}, which provides a more sophisticated
#'   approach to smoothing out spatial variability but may be affected by
#'   sparsity due to sampling variability (especially sparsity at the index
#'   point), and is computationally slower. Default = \code{uniform}.
#' 
#' @param bandwidth Bandwidth parameter for smoothing, expressed as proportion
#'   of width or height (whichever is greater) of tissue area. Smoothing is
#'   performed across neighboring values at each point. For \code{method =
#'   "uniform"}, the bandwidth represents the radius of a circle, and unweighted
#'   averages are calculated across points within this circle. For \code{method
#'   = "kernel"}, the averaging is weighted by distances scaled using a
#'   truncated exponential kernel applied to Euclidean distances. For example, a
#'   bandwidth of 0.05 will smooth values across neighbors weighted by distances
#'   scaled using a truncated exponential kernel with length scale equal to 5%
#'   of the width or height (whichever is greater) of the tissue area. Weights
#'   for \code{method = "kernel"} are truncated at small values for
#'   computational efficiency. Default = 0.05.
#' 
#' @param truncate Truncation threshold parameter if \code{method = "kernel"}.
#'   Kernel weights below this value are set to zero for computational
#'   efficiency. Only used for \code{method = "kernel"}. Default = 0.05.
#' 
#' @param assay_name_smooth Name of assay to store smoothed values. Default =
#'   \code{counts_smooth}. Alternatively, set to \code{logcounts_smooth} if
#'   logcounts were provided instead of raw counts for the input values.
#' 
#' @param calc_logcounts Whether to calculate logcounts using smoothed values of
#'   raw counts and store these in new assay named \code{logcounts_smooth} in
#'   output object. Uses methods from \code{scran} package to calculate
#'   normalization and log-transformation. Default = TRUE.
#' 
#' 
#' @return Returns the \code{SpatialExperiment} object with new assays (default
#'   names \code{counts_smooth} and \code{logcounts_smooth}) containing
#'   spatially smoothed expression values, which can be used as the input for
#'   further downstream analyses such as clustering.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment assays 'assays<-' assayNames
#' @importFrom spdep dnearneigh nbdists
#' @importFrom methods is as
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom scuttle normalizeCounts
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
#' # check
#' assayNames(spe)
#' 
smoothclust <- function(input, method = c("uniform", "kernel"), 
                        assay_name = "counts", 
                        bandwidth = 0.1, truncate = 0.05, 
                        assay_name_smooth = "counts_smooth", 
                        calc_logcounts = TRUE) {
  
  method <- match.arg(method, c("uniform", "kernel"))
  
  spe <- input
  
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(assay_name %in% assayNames(spe))
  
  vals <- assays(spe)[[assay_name]]
  
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
  
  # calculate smoothed values
  # note: using dense matrix to ensure zeros are included in averaging
  vals <- as.matrix(vals)
  vals_smooth <- matrix(as.numeric(NA), nrow = nrow(spe), ncol = ncol(spe))
  
  pb <- txtProgressBar(0, ncol(vals_smooth), style = 3)
  
  if (method == "uniform") {
    for (i in seq_len(ncol(vals_smooth))) {
      setTxtProgressBar(pb, i)
      # extract values
      vals_sub <- vals[, neigh[[i]]]
      # calculate average
      vals_smooth[, i] <- rowMeans(vals_sub)
    }
  }
  
  if (method == "kernel") {
    for (i in seq_len(ncol(vals_smooth))) {
      setTxtProgressBar(pb, i)
      # extract values and weights
      vals_sub <- vals[, neigh[[i]]]
      weights_rep <- t(replicate(nrow(vals_sub), weights[[i]]))
      stopifnot(all(dim(vals_sub) == dim(weights_rep)))
      # calculate weighted average
      out <- rowSums(vals_sub * weights_rep) / sum(weights[[i]])
      vals_smooth[, i] <- out
    }
  }
  
  close(pb)
  
  stopifnot(nrow(vals_smooth) == nrow(spe))
  stopifnot(ncol(vals_smooth) == ncol(spe))
  rownames(vals_smooth) <- rownames(spe)
  colnames(vals_smooth) <- colnames(spe)
  
  # store in object
  assays(spe)[[assay_name_smooth]] <- vals_smooth
  
  # add logcounts
  if (calc_logcounts) {
    lc_smooth <- normalizeCounts(vals_smooth)
    assays(spe)[["logcounts_smooth"]] <- lc_smooth
  }
  
  spe
}
