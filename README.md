# smoothclust

[![R build status](https://github.com/lmweber/smoothclust/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/smoothclust/actions)


## Overview

`smoothclust` is a method for segmentation of spatial domains and spatially-aware clustering in spatial transcriptomics data. The method generates spatial domains with smooth boundaries by smoothing gene expression profiles across neighboring spatial locations, followed by unsupervised clustering. Spatial domains consisting of consistent mixtures of cell types may then be further investigated by applying cell type compositional analyses or differential analyses.


## Installation

The development version of the package can be installed from GitHub:

```
remotes::install_github("lmweber/smoothclust")
```

The release version will be made available from Bioconductor.


## Tutorial

For a tutorial and example workflow, see the package vignette.


## Citation

In preparation.
