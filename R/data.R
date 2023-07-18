#' SingleCellExperiment containing example inputs to GloScope
#' @rdname data
#' @name example_SCE
#' @aliases example_SCE_small
#' @description `example_SCE` is a SingleCellExperiment object which
#'   contains PCA embeddings and metadata for PBMCs from 20 COVID-infected
#'   and healthy control patients. Each sample is reduced
#'   to a random subset of 500 cells, for a total of 10,000 cells. The 
#'   `colData` slot of the object contains the metadata for each cell, 
#'   its sample ID and phenotype. The dimensionality reductions slot 
#'   contains the first 50 PCs, and these embeddings are provided by
#'   the authors of "Single-cell multi-omics analysis of the immune
#'   response in COVID-19" (Stephenson et al., 2021; doi: 10.1038/s41591-021-01329-2).
#' @description `example_SCE_small` is a SingleCellExperiment with the same structure as
#'   `example_SCE`, but only containing data from the first five samples. This is a
#'   smaller set for examples.
#' @format A SingleCellExperiment object with metadata and PCA embeddings 
#' @return A SingleCellExperiment object
#' @examples
#' # Code to create the small SCE from the full sample
#' # Reduction to 5 samples demonstrates data extraction from SCE objects
#' data(example_SCE)
#' sample_ids <- SingleCellExperiment::colData(example_SCE)$sample_id
#' whKeep <- which(sample_ids %in% unique(sample_ids)[seq_len(5)])
#' example_SCE_small <- SingleCellExperiment::SingleCellExperiment(
#'  assays=list(counts=matrix(rep(0,2500),ncol=2500)),
#'  colData=SingleCellExperiment::colData(example_SCE)[whKeep,],
#'  reducedDims=list("PCA"=SingleCellExperiment::reducedDim(example_SCE,"PCA")[whKeep,]))
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SingleCellExperiment colData

NULL
