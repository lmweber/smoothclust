#' Function for smoothness metric
#' 
#' Function for clustering smoothness evaluation metric
#' 
#' Function to calculate smoothness of cluster boundaries, defined as average
#' number of nearest neighbors per point that are from a different cluster.
#' 
#' 
#' @param spatialcoords Numeric matrix containing spatial coordinates of points,
#'   formatted as nrow = number of points, ncol = 2 (assuming x and y
#'   dimensions). For example, `spatialcoords = spatialCoords(spe)`.
#' 
#' @param labels Numeric vector of cluster labels for each point. For example,
#'   `labels <- as.numeric(colData(spe)$label)`.
#' 
#' @param k Number of k nearest neighbors to use in calculation. Default = 6
#'   (from 10x Genomics Visium platform).
#' 
#' 
#' @return Returns value of evaluation metric.
#' 
#' 
#' @importFrom spdep knearneigh
#' 
#' @export
#' 
#' @examples
#' # Example to be shown in vignette.
#' 
smoothness_metric <- function(spatialcoords, labels, k = 6) {
  
  stopifnot(length(labels) == nrow(spatialcoords))
  stopifnot(ncol(spatialcoords) != 2)
  
  # calculate k nearest neighbors for each point
  neigh <- knearneigh(spatialcoords, k = k)$nn
  
  # calculate ordered columns of cluster labels
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
