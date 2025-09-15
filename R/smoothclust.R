#' smoothclust
#' 
#' Method for identification of spatial domains and spatially-aware clustering.
#' 
#' Method for identification of spatial domains and spatially-aware clustering
#' in spatial transcriptomics data.
#' 
#' Method for identification of spatial domains and spatially-aware clustering
#' in spatial transcriptomics data. The method generates spatial domains with
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
#'   averages. This argument is only used if the input is a
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
#' @importFrom BiocNeighbors findNeighbors findKNN
#' @importFrom Matrix sparseMatrix
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
#' # run smoothclust using default parameters
#' spe <- smoothclust(spe)
#' 
#' # see vignette for extended example including downstream analyses
#' 
smoothclust <- function(input, assay_name = "counts", spatial_coords = NULL, 
                        method = c("uniform", "kernel", "knn"), 
                        bandwidth = 0.05, k = 18, truncate = 0.05) {
  
  method <- match.arg(method, c("uniform", "kernel", "knn"))
  
  stopifnot(is.character(assay_name) && length(assay_name) == 1)
  stopifnot(is.numeric(bandwidth))
  stopifnot(is.numeric(k))
  stopifnot(is.numeric(truncate))
  
  if (is(input, "SpatialExperiment")) {
    spe <- input
    stopifnot(assay_name %in% assayNames(spe))
    vals <- assays(spe)[[assay_name]]
    spatial_coords <- spatialCoords(spe)
  } else {
    stopifnot(is.numeric(input))
    vals <- input
  }
  
  stopifnot(!is.null(spatial_coords), 
            is.numeric(spatial_coords), 
            is.matrix(spatial_coords), 
            ncol(spatial_coords) == 2)
  
  # convert vals to CsparseMatrix for efficient multiplication
  vals <- as(vals, "CsparseMatrix")
  N <- ncol(vals)
  
  if (method %in% c("uniform", "kernel")) {
    # convert bandwidth to same units as distances
    range_x <- abs(diff(range(spatial_coords[, 1])))
    range_y <- abs(diff(range(spatial_coords[, 2])))
    range_max <- max(range_x, range_y)
    bandwidth_scaled <- bandwidth * range_max
  }
  
  # --- method-specific steps to construct weights matrix (W) ---
  
  if (method == "uniform") {
    # 1. fast neighbor search
    nn_data <- findNeighbors(spatial_coords, threshold = bandwidth_scaled, 
                             get.index = TRUE, get.distance = FALSE)
    nn_list <- nn_data$index
    
    # 2. get number of neighbors for each spot
    n_neighbors <- lengths(nn_list)
    
    # 3. prepare indices and values for sparse weights matrix W
    # row indices: the neighbors themselves
    i_idx <- unlist(nn_list, use.names = FALSE)
    # column indices: the spot being considered (repeated for each of its neighbors)
    j_idx <- rep(seq_along(nn_list), n_neighbors)
    # values: 1 / (number of neighbors for that column)
    x_val <- rep(1 / n_neighbors, n_neighbors)
    
  } else if (method == "kernel") {
    # 1. determine max search distance based on truncation threshold
    # to avoid slow all-vs-all distance calculation
    max_dist <- -bandwidth_scaled * log(truncate)
    
    # 2. fast neighbor search within this radius
    nn_data <- findNeighbors(spatial_coords, threshold = max_dist, 
                             get.index = TRUE, get.distance = TRUE)
    
    # 3. calculate raw exponential kernel weights
    i_idx <- unlist(nn_data$index, use.names = FALSE)
    j_idx <- rep(seq_along(nn_data$index), lengths(nn_data$index))
    dists <- unlist(nn_data$distance, use.names = FALSE)
    
    raw_weights <- exp(-dists / bandwidth_scaled)
    
    # 4. normalize weights so each column in W sums to 1
    col_sums <- as.vector(tapply(raw_weights, j_idx, sum))
    x_val <- raw_weights / col_sums[j_idx]
    
  } else if (method == "knn") {
    # 1. fast k-nearest neighbor search (k+1 to include self)
    nn_data <- findKNN(spatial_coords, k = k + 1, 
                       get.index = TRUE, get.distance = FALSE)
    
    # 2. prepare indices and values for W; weight is uniform 1/(k+1)
    i_idx <- as.vector(t(nn_data$index))
    j_idx <- rep(seq_len(N), each = k + 1)
    x_val <- rep(1 / (k + 1), length(i_idx))
  }
  
  # --- construct weights matrix W and perform single matrix multiplication ---
  
  # 1. construct weights matrix
  W <- sparseMatrix(i = i_idx, j = j_idx, x = x_val, dims = c(N, N), repr = "C")
  
  # 2. perform smoothing operation in one matrix multiplication step
  vals_smooth <- vals %*% W
  
  # --- return results ---
  
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
