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
#' @param bandwidth Parameter defining the bandwidth for smoothing, expressed as
#'   the proportion of the width or height (whichever is greater) of the tissue
#'   area. For example, a bandwidth of 0.1 will smooth logcounts across values
#'   measured within a circle of radius equal to 10% of the width or height
#'   (whichever is greater) of the tissue area.
#' 
#' 
#' @return Returns the \code{SpatialExperiment} object with a new assay named
#'   \code{logcounts_smooth} containing spatially smoothed logcounts values.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assays 'assays<-'
#' @importFrom spdep dnearneigh
#' @importFrom methods is as
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
smoothclust <- function(input, bandwidth = 0.1) {
  
  spe <- input
  
  stopifnot(is(spe, "SpatialExperiment"))
  
  spatialcoords <- spatialCoords(spe)
  
  # convert bandwidth argument to distance units
  range_x <- abs(diff(range(spatialcoords[, 1])))
  range_y <- abs(diff(range(spatialcoords[, 2])))
  range_max <- max(range_x, range_y)
  bandwidth_dist <- bandwidth * range_max
  
  # calculate neighbors
  neigh <- dnearneigh(spatialcoords, d1 = 0, d2 = bandwidth_dist)
  # include self within set of neighbors for each point
  stopifnot(length(neigh) == ncol(spe))
  self <- as.list(seq_along(neigh))
  neigh <- mapply(c, self, neigh, SIMPLIFY=FALSE)
  # format as matrix (with NAs to fill)
  # using utility functions
  neigh_mx <- do.call(.cbind.fill, c(neigh, fill = NA))
  rownames(neigh_mx) <- NULL
  colnames(neigh_mx) <- NULL
  
  # to do: alternative using kernel smoothing
  
  # calculate average logcounts across neighbors
  # note: missing values (non-detected genes) are interpreted as zeros when averaging
  # note: formatting as dense 3-D array / tensor
  logcounts_neigh <- array(NA, dim = c(nrow(spe), ncol(spe), nrow(neigh_mx)))
  lc <- as.matrix(logcounts(spe))
  for (i in seq_len(ncol(spe))) {
    logcounts_neigh[, i, ] <- lc[, neigh_mx[, i]]
  }
  
  # slightly slow step (runtime: seconds)
  logcounts_smooth <- apply(logcounts_neigh, c(1, 2), mean, na.rm = TRUE)
  
  stopifnot(nrow(logcounts_smooth) == nrow(spe))
  stopifnot(ncol(logcounts_smooth) == ncol(spe))
  rownames(logcounts_smooth) <- rownames(spe)
  colnames(logcounts_smooth) <- colnames(spe)
  
  # store in object
  assays(spe)[["logcounts_smooth"]] <- logcounts_smooth
  
  spe
}
