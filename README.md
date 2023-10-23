# GloScope Package

This repository implements the `GloScope` methodology for sample-level analysis of scRNA-Seq data. The details of this methodology can be found in the following pre-print:

[H. Wang, W. Torous, B. Gong, E. Purdom (2023).
Visualizing scRNA-Seq Data at Population Scale with GloScope. bioRxiv.](https://doi.org/10.1101/2023.05.29.542786)

A vignette demonstrating the basic usage of this package can be found in the [vignettes folder](https://github.com/epurdom/GloScope/tree/main/vignettes) of this repository.

## Installation

This package is available on Bioconductor, and the following R code can be used to install it.

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GloScope")
```
We expect this package to be available on BioConductor shortly.
