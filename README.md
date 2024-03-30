# smoothclust

[![R build status](https://github.com/lmweber/smoothclust/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/smoothclust/actions)


## Overview

`smoothclust` is a method for segmentation of spatial domains and spatially-aware clustering in spatial transcriptomics data. The method generates spatial domains with smooth boundaries by smoothing gene expression profiles across neighboring spatial locations, followed by unsupervised clustering. Spatial domains consisting of consistent mixtures of cell types may then be further investigated by applying cell type compositional analyses or differential analyses.


## Installation

The `smoothclust` package has been submitted to Bioconductor. The release version of the package will be installable from Bioconductor as follows (using R version 4.4 and Bioconductor version 3.19, available from approximately 2024-05-01 onwards). Additional details will be shown on the [Bioconductor](https://bioconductor.org/packages/smoothclust) package landing page.

```
install.packages("BiocManager")
BiocManager::install("smoothclust")
```

The latest development version of the package will also be installable from the [development version of Bioconductor](https://contributions.bioconductor.org/use-devel.html) in the future.

Alternatively, the current development version of the package can be installed from GitHub as follows. (Note the argument `ref = R_release`, which is required if you are using R version 4.3, which is the current stable release version of R.)

```
remotes::install_github("lmweber/smoothclust", ref = "R_release")
```


## Tutorial

For a tutorial and example workflow, see the package vignette.


## Citation

In preparation.
