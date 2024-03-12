# smoothclust function versions using dense matrices (for debugging)

# calculate smoothed values
# note: using dense matrices

vals <- as.matrix(vals)
vals_smooth <- matrix(as.numeric(NA), nrow = nrow(vals), ncol = ncol(vals))

pb <- txtProgressBar(0, ncol(vals_smooth), style = 3)

if (method == "uniform") {
  for (i in seq_len(ncol(vals_smooth))) {
    setTxtProgressBar(pb, i)
    # extract values
    vals_sub <- vals[, neigh[[i]], drop = FALSE]
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

if (method == "knn") {
  stopifnot(nrow(neigh) == ncol(vals_smooth))
  for (i in seq_len(ncol(vals_smooth))) {
    setTxtProgressBar(pb, i)
    # extract values
    vals_sub <- vals[, neigh[i, ], drop = FALSE]
    # calculate average
    vals_smooth[, i] <- rowMeans(vals_sub)
  }
}

close(pb)
