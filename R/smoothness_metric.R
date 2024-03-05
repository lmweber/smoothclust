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
#' library(STexampleData)
#' library(scran)
#' library(scater)
#' 
#' # load data
#' spe <- Visium_humanDLPFC()
#' # keep spots over tissue
#' spe <- spe[, colData(spe)$in_tissue == 1]
#' 
#' # run smoothclust
#' # using "knn" method for faster runtime in this example
#' # see vignette for example using default method and downstream analyses
#' spe <- smoothclust(spe, method = "knn", k = 6)
#' 
#' # calculate logcounts
#' spe <- logNormCounts(spe)
#' 
#' # preprocessing steps for downstream clustering
#' # remove mitochondrial genes
#' is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
#' spe <- spe[!is_mito, ]
#' # select top highly variable genes (HVGs)
#' dec <- modelGeneVar(spe)
#' top_hvgs <- getTopHVGs(dec, prop = 0.1)
#' spe <- spe[top_hvgs, ]
#' 
#' # dimensionality reduction
#' set.seed(123)
#' spe <- runPCA(spe)
#' 
#' # run k-means clustering
#' set.seed(123)
#' k <- 5
#' clust <- kmeans(reducedDim(spe, "PCA"), centers = k)$cluster
#' colLabels(spe) <- factor(clust)
#' 
#' # calculate smoothness metric
#' smoothness_metric(spatialCoords(spe), as.numeric(colData(spe)$label))
#' 
smoothness_metric <- function(spatialcoords, labels, k = 6) {
  
  stopifnot(length(labels) == nrow(spatialcoords))
  stopifnot(ncol(spatialcoords) == 2)
  
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
