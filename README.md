# smoothclust

[![R build status](https://github.com/lmweber/smoothclust/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/smoothclust/actions)


## Overview

Tools for spatial domain analysis in spatial transcriptomics data. Includes spatial smoothing of gene expression profiles for downstream clustering and a boundary density (BD) metric for summarizing spatial domain organization. The smoothing function returns processed expression values for subsequent clustering. BD requires only spatial coordinates and categorical labels, and provides local and sample-level summaries of label transitions.


## Installation

The `smoothclust` package can be installed from Bioconductor as follows. Additional details are shown on the [Bioconductor](https://bioconductor.org/packages/smoothclust) package landing page.

```
install.packages("BiocManager")
BiocManager::install("smoothclust")
```

The latest development version can also be installed from the [development version of Bioconductor](https://contributions.bioconductor.org/use-devel.html) or from GitHub as follows.

```
remotes::install_github("lmweber/smoothclust")
```

or

```
devtools::install_github("lmweber/smoothclust")
```


## Tutorial

For a tutorial and example workflow, see the package vignette.


## Citation

Rusnak and Weber (2026), *smoothclust: Spatial smoothing and boundary-density summaries for spatial transcriptomics*, bioRxiv (preprint). [Available here](https://www.biorxiv.org/content/10.64898/2026.10.26.761912v1).
