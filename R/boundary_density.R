#' Boundary density metric
#' 
#' Function for boundary density metric
#' 
#' Function to calculate the raw boundary density metric, defined as the
#' average fraction of nearest neighbors per point that are from a different
#' cluster. This metric can be used to quantify and compare the relative density
#' of the boundaries of clusters or spatial domains. Optional composition
#' adjustments can be used to compare the raw boundary density to an expected
#' value under random relabeling.
#' 
#' 
#' @param spatial_coords Numeric matrix containing spatial coordinates of
#'   points, formatted as nrow = number of points, ncol = 2 (assuming x and y
#'   dimensions). For example, `spatial_coords = spatialCoords(spe)` if using a
#'   \code{SpatialExperiment} object.
#' 
#' @param labels Atomic vector or factor containing cluster labels for each
#'   point. Missing values are not allowed. For example, `labels <-
#'   colData(spe)$label` if using a \code{SpatialExperiment} object.
#' 
#' @param k Number of k nearest neighbors to use in calculation. Default = 6
#'   (e.g. for hexagonal arrangement in 10x Genomics Visium platform).
#' 
#' @param n_threads Number of threads to use for nearest-neighbor searches.
#'   Default = 1.
#' 
#' @param adjust Composition adjustment to calculate. Options are \code{none},
#'   \code{analytic}, and \code{permutation}. With \code{none}, only the raw
#'   boundary density is returned. With \code{analytic}, the raw boundary
#'   density is compared to a deterministic expectation based on the observed
#'   label counts. With \code{permutation}, labels are randomly permuted across
#'   spatial locations \code{n_permutations} times. Default = \code{none}.
#' 
#' @param n_permutations Number of random label permutations to use when
#'   \code{adjust = "permutation"}. Default = 1000.
#' 
#' @param seed Optional random seed to use when \code{adjust = "permutation"}.
#'   Default = \code{NULL}.
#' 
#' 
#' @return Returns a list containing values at each point (i.e. the number of
#'   nearest neighbors that are from a different cluster, the number of nearest
#'   neighbors, and the local raw boundary density) as well as the mean
#'   discordant neighbor count and the sample-level raw boundary density. If
#'   \code{adjust = "analytic"}, the list also contains the expected boundary
#'   density, relative boundary density, and excess boundary density. If
#'   \code{adjust = "permutation"}, the list also contains the permutation
#'   mean, standard deviation, relative boundary density, standardized score,
#'   and number of permutations.
#' 
#' 
#' @importFrom BiocNeighbors findKNN
#' 
#' @export
#' 
#' @examples
#' spatial_coords <- cbind(c(0, 1, 3, 6), 0)
#' labels <- c("A", "B", "A", "A")
#' 
#' # calculate raw boundary density metric
#' res <- boundary_density(spatial_coords, labels, k = 2)
#' str(res)
#' res$n_discordant
#' res$mean_discordant
#' res$boundary_density
#' 
#' # calculate analytic composition adjustment
#' res_adj <- boundary_density(spatial_coords, labels, k = 2, 
#'                             adjust = "analytic")
#' res_adj$expected_boundary_density
#' res_adj$relative_boundary_density
#' 
boundary_density <- function(spatial_coords, labels, k = 6, n_threads = 1, 
                             adjust = c("none", "analytic", "permutation"), 
                             n_permutations = 1000, seed = NULL) {
  
  adjust <- match.arg(adjust)
  
  stopifnot(!is.null(spatial_coords), 
            is.numeric(spatial_coords), 
            is.matrix(spatial_coords), 
            ncol(spatial_coords) == 2, 
            all(is.finite(spatial_coords)))
  stopifnot((is.atomic(labels) || is.factor(labels)) && 
              is.null(dim(labels)))
  stopifnot(length(labels) == nrow(spatial_coords))
  stopifnot(!anyNA(labels))
  stopifnot(is.numeric(k) && length(k) == 1 && is.finite(k) && 
              k > 0 && k == floor(k) && k < nrow(spatial_coords))
  stopifnot(is.numeric(n_threads) && length(n_threads) == 1 && 
              is.finite(n_threads) && n_threads > 0 && 
              n_threads == floor(n_threads))
  if (adjust == "permutation") {
    stopifnot(is.numeric(n_permutations) && length(n_permutations) == 1 && 
                is.finite(n_permutations) && n_permutations > 1 && 
                n_permutations == floor(n_permutations))
    stopifnot(is.null(seed) || 
                (is.numeric(seed) && length(seed) == 1 && is.finite(seed) && 
                   seed >= 0 && seed == floor(seed)))
  }
  
  # --- fast k-nearest neighbor search ---
  
  # search for k nearest neighbors
  nn_data <- findKNN(spatial_coords, k = k, 
                     get.index = TRUE, get.distance = FALSE, 
                     num.threads = n_threads)
  neigh <- nn_data$index
  
  # --- vectorized label lookup and calculation ---
  
  # create matrix of neighbor labels using a single matrix-indexing operation
  neigh_labels <- matrix(labels[neigh], ncol = ncol(neigh))
  
  # compare 'labels' vector against each column of 'neigh_labels'
  n_discordant <- rowSums(labels != neigh_labels)
  n_neighbors <- rep.int(ncol(neigh), nrow(neigh))
  local_boundary_density <- n_discordant / n_neighbors
  boundary_density <- mean(local_boundary_density)
  
  # --- return results ---
  
  # return local and sample-level values
  out <- list(n_discordant = n_discordant, 
              n_neighbors = n_neighbors, 
              local_boundary_density = local_boundary_density, 
              mean_discordant = mean(n_discordant), 
              boundary_density = boundary_density)
  
  if (adjust == "analytic") {
    label_counts <- as.numeric(table(labels))
    n <- length(labels)
    expected_boundary_density <- 
      1 - sum(label_counts * (label_counts - 1)) / (n * (n - 1))
    relative_boundary_density <- 
      if (expected_boundary_density > 0) {
        boundary_density / expected_boundary_density
      } else {
        NA_real_
      }
    excess_boundary_density <- 
      boundary_density - expected_boundary_density
    
    out <- c(out, 
             list(expected_boundary_density = expected_boundary_density, 
                  relative_boundary_density = relative_boundary_density, 
                  excess_boundary_density = excess_boundary_density))
  } else if (adjust == "permutation") {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), 
                add = TRUE)
      } else {
        on.exit(if (exists(".Random.seed", envir = .GlobalEnv, 
                           inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }, add = TRUE)
      }
      set.seed(seed)
    }
    
    permutation_boundary_density <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      permuted_labels <- sample(labels)
      perm_neigh_labels <- matrix(permuted_labels[neigh], ncol = ncol(neigh))
      permutation_boundary_density[i] <- 
        mean(rowSums(permuted_labels != perm_neigh_labels) / n_neighbors)
    }
    
    permutation_mean <- mean(permutation_boundary_density)
    permutation_sd <- sqrt(
      sum((permutation_boundary_density - permutation_mean)^2) / 
        (n_permutations - 1))
    permutation_relative_boundary_density <- 
      if (permutation_mean > 0) {
        boundary_density / permutation_mean
      } else {
        NA_real_
      }
    permutation_z <- 
      if (permutation_sd > 0) {
        (boundary_density - permutation_mean) / permutation_sd
      } else {
        NA_real_
      }
    
    out <- c(out, 
             list(permutation_mean = permutation_mean, 
                  permutation_sd = permutation_sd, 
                  permutation_relative_boundary_density = 
                    permutation_relative_boundary_density, 
                  permutation_z = permutation_z, 
                  n_permutations = n_permutations))
  }
  
  out
}
