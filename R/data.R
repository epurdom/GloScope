#' example of data frame as the input for gloscope
#' @rdname data
#' @name example_data
#' @aliases example_small_data
#' @description `example_data` is a list containing PCA embeddings and metadata
#'   of cells from 20 COVID-infected and healthy PBMC samples. Each original
#'   sample is reduced to a subset of 500 cells, for a total of 10,000 cells.
#'   The list element `metadata` is a data.frame containing the associated
#'   sample ID and phenotype for each cell, and the element `pca_embeddings`
#'   contains the first 50 PCs. These embeddings are provided by the authors of
#'   "Single-cell multi-omics analysis of the immune response in COVID-19"
#'   (Stephenson et al., 2021; doi: 10.1038/s41591-021-01329-2).
#' @description `example_small_data` is a list of the same structure as
#'   `example_data`, only consisting of the first five samples. This is a
#'   smaller set for examples.
#' @format A list containing two data.frame: `metadata` (10000 x 2) and
#' `pca_embeddings` (10000 x 50)
#' @examples
#' # Code to make the small set of data from `example_data`
#' # Reduce to only 5 samples for small illustration
#' data(example_data)
#' sample_ids <- example_data$metadata$sample_id
#' whKeep<-which(sample_ids %in% unique(sample_ids)[seq_len(5)])
#' example_small_data<-list(
#'  metadata=example_data$metadata[whKeep,],
#'  pca_embeddings=example_data$pca_embeddings[whKeep,])
#'
#'

NULL
