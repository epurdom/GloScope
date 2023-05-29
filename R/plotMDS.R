#' @title Plot the multidimensional scaling of the GloScope represenation
#'
#' @description This function creates a multidimensional scaling plot for
#'  a set of samples using their GloScope divergence. Each sample's scatter will
#' be color-coded based on their phenotype.
#'
#' @param dist_mat The divergence matrix output of `gloscope()`
#' @param metadata_df A data frame contains each sample's metadata. Note this is NOT at the cell-level.
#' @param sample_id The column name or index in metadata_df that contains the sample ID
#' @param group_id The column name or index in metadata_df that contains the patient condition
#' @param n Number of MDS dimension to generate, default = 10
#' @return A list containing the MDS embedding and plot of the distance matrix
#'
#' @examples
#' \donttest{
#' data(example_data)
#' sample_ids <- example_data$metadata$sample_id
#' pca_embeddings <- example_data$pca_embeddings
#' pca_embeddings_subset <- pca_embeddings[,1:10] # select the first 10 PCs
#' dist_result <- gloscope(pca_embeddings_subset, sample_ids,
#'                     BPPARAM = BiocParallel::SerialParam(RNGseed=2))
#' mds_result <- plotMDS(dist_mat = dist_result, metadata_df =  unique(example_data$metadata),
#' "sample_id", "phenotype",n=2)
#' mds_result$plot
#' mds_result$mds
#' }
#' @import ggplot2
#' @importFrom MASS isoMDS
#' @export
#'
plotMDS <- function(dist_mat, metadata_df, sample_id, group_id, n=10){
  if(nrow(dist_mat)!= nrow(metadata_df)){
    stop(paste("Not consistent patient number. Make sure your
               distance matrix and meta info have the same patient
               number."))
  }
  if(length(intersect(rownames(dist_mat), metadata_df[, sample_id]))!= nrow(dist_mat)){
    stop(paste("Not consistent patient IDs. Make sure your
               distance matrix and meta info have the same patients."))
  }
  if(!(identical(rownames(dist_mat), metadata_df[, sample_id]))){
    metadata_df <- metadata_df[match(rownames(dist_mat), metadata_df[,sample_id]),]
  }

  fit_df <- MASS::isoMDS(dist_mat, k = n, trace = FALSE)
  colnames(fit_df$points) <- paste0("Coordinate",seq_len(n))
  mds_df <- cbind(metadata_df, fit_df$points)

  colour_palette <- .paletteBig()
  mds_plot <- ggplot2::ggplot(mds_df, aes(x = .data$Coordinate1, y = .data$Coordinate2,
    color = .data[[group_id]])) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values=colour_palette) +
    ggplot2::theme_bw()

  return(list(mds = mds_df, plot = mds_plot))
}
