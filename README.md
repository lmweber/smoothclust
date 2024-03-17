# smoothclust

[![R build status](https://github.com/lmweber/smoothclust/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/smoothclust/actions)


## Overview

`smoothclust` is a method for segmentation of spatial domains and spatially-aware clustering in spatial transcriptomics data. The method generates spatial domains with smooth boundaries by smoothing gene expression profiles across neighboring spatial locations, followed by unsupervised clustering. Spatial domains consisting of consistent mixtures of cell types may then be further investigated by applying cell type compositional analyses or differential analyses.


## Installation

The package has been submitted to Bioconductor, and the release version of the package will be installable from Bioconductor in the future.

The latest development version of the package can be installed from GitHub using the code below. (Note the argument `ref = R_release`, which is required if you are using a recent release version of R (e.g. 4.3). If you try to install without this argument, you will also need the latest development version of R (4.4), which is required for the current development version of Bioconductor.)

```
remotes::install_github("lmweber/smoothclust", ref = "R_release")
```


## Tutorial

For a tutorial and example workflow, see the package vignette.


## Citation

In preparation.
