#' Function for evaluation metric
#' 
#' Function for clustering smoothness evaluation metric
#' 
#' Function to calculate smoothness of cluster boundaries, defined as average
#' number of nearest neighbors per point that are from a different cluster.
#' 
#' 
#' @param input Input object, assumed to be a \code{SpatialExperiment} object
#'   containing a column of cluster labels in `colData`.
#' 
#' @param label Name of column in `colData` containing cluster labels. Default =
#'   "label".
#' 
#' @param k Number of k nearest neighbors to use in calculation. Default = 6
#'   (from 10x Genomics Visium platform).
#' 
#' 
#' @return Returns value of evaluation metric.
#' 
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom spdep knearneigh
#' 
#' @export
#' 
#' @examples
#' # Example shown in vignette.
#' 
smoothness_metric <- function(input, label = "label", k = 6) {
  
  spe <- input
  stopifnot(is(spe, "SpatialExperiment"))
  
  # extract spatial coordinates for each point
  spatialcoords <- spatialCoords(spe)
  
  # calculate k nearest neighbors for each point
  neigh <- knearneigh(spatialcoords, k = k)$nn
  
  # calculate ordered columns of cluster labels
  labels <- as.numeric(colData(spe)$label)
  neigh_labels <- matrix(NA, nrow = nrow(neigh), ncol = ncol(neigh))
  for (i in seq_len(ncol(neigh_labels))) {
    neigh_labels[, i] <- labels[neigh[, i]]
  }
  
  # calculate number of non-matching labels
  stopifnot(length(labels) == nrow(neigh_labels))
  res <- rep(0, length(labels))
  for (i in seq_len(ncol(neigh_labels))) {
    res <- res + as.numeric(labels != neigh_labels[, i])
  }
  
  # average value
  mean(res)
}
