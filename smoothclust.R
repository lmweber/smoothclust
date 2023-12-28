# set up object and calculate HVGs

library(STexampleData)
spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

library(scran)
spe <- logNormCounts(spe)

is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
spe <- spe[!is_mito, ]

# keep full object for plotting
spe_full <- spe

dec <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
spe <- spe[top_hvgs, ]

dim(spe)


# calculate neighbors (including self)

colData(spe)$neighbors <- as.list(rep(NA, ncol(spe)))

arrayrow <- colData(spe)$array_row
arraycol <- colData(spe)$array_col

# slow step (can improve)
for (i in 1:ncol(spe)) {
  neighbors <- i
  for (j in setdiff(1:ncol(spe), i)) {
    # defining neighbors as within row and column distance threshold
    if (abs(arrayrow[i] - arrayrow[j]) <= 5 & (abs(arraycol[i] - arraycol[j])) <= 5) {
      neighbors <- c(neighbors, j)
    }
  }
  # neighbors stored as indices
  colData(spe)$neighbors[[i]] <- neighbors
  print(i)
}

head(colData(spe))
head(colData(spe)$neighbors)

# check
colData(spe)$neighbors[[1]]


# calculate average logcounts across neighbors (vectorized for faster runtime)
# (alternatively: use kernels for more sophisticated approach)

smoothed_logcounts <- matrix(NA, nrow = nrow(spe), ncol = ncol(spe))

for (i in 1:ncol(smoothed_logcounts)) {
  smoothed_logcounts[, i] <- rowMeans(logcounts(spe)[, colData(spe)$neighbors[[i]], drop = FALSE])
  print(i)
}

rownames(smoothed_logcounts) <- rownames(spe)
colnames(smoothed_logcounts) <- colnames(spe)

assays(spe)[["smoothed_logcounts"]] <- smoothed_logcounts

assayNames(spe)


### too slow (nested loop)

# smoothed_logcounts <- matrix(NA, nrow = nrow(spe), ncol = ncol(spe))
# 
# for (g in 1:nrow(smoothed_logcounts)) {
#   for (j in 1:ncol(smoothed_logcounts)) {
#     smoothed_logcounts[g, j] <- mean(logcounts(spe)[g, colData(spe)$neighbors[[j]]])
#     print(j)
#   }
#   print(g)
# }


# plots for checking (PCP4 gene)

library(ggplot2)

df <- cbind(as.data.frame(spatialCoords(spe)), 
            gene = smoothed_logcounts["ENSG00000183036", ], 
            gene_original = logcounts(spe_full)[33335, ])

colnames(df) <- c("x", "y", "gene", "gene_original")

ggplot(df, aes(x = x, y = y, color = gene)) + geom_point(size = 0.5) + scale_y_reverse() + coord_fixed()
ggplot(df, aes(x = x, y = y, color = gene_original)) + geom_point(size = 0.5) + scale_y_reverse() + coord_fixed()


# dimensionality reduction and clustering

library(scater)
library(scran)

# compute PCA
set.seed(123)
spe <- runPCA(spe, exprs_values = "smoothed_logcounts")
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)


# graph-based clustering
set.seed(123)
#k <- 10
k <- 100
g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe) <- factor(clus)


# alternatively: k-means clustering (works better)
k <- 6
set.seed(100)
clust <- kmeans(reducedDim(spe, "PCA"), centers = k)$cluster
table(clust)
colLabels(spe) <- factor(clust)


# plots
library(ggspavis)

plotSpots(spe, annotate = "label", palette = "libd_layer_colors")
plotSpots(spe, annotate = "label", palette = unname(palette.colors(36, "Polychrome 36")))

plotSpots(spe, annotate = "ground_truth", palette = "libd_layer_colors")


# ARI

library(mclust)
adjustedRandIndex(colData(spe)$label, colData(spe_full)$ground_truth)

