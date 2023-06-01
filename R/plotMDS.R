#' @title Plot the multidimensional scaling of the GloScope represenation
#'
#' @description This function creates a multidimensional scaling plot for
#'  a set of samples using their GloScope divergence. Each sample's scatter will
#' be color-coded based on their phenotype. The function calls the `isoMDS` function
#' from the `MASS` package.
#'
#' @param dist_mat The divergence matrix output of `gloscope()`
#' @param metadata_df A data frame contains each sample's metadata. Note this is NOT at the cell-level.
#' @param sample_id The column name or index in metadata_df that contains the sample ID
#' @param group_id The column name or index in metadata_df that contains the patient condition
#' @param k Number of MDS dimension to generate, default = 10
#' @return A list containing the MDS embedding and plot of the distance matrix
#' \itemize{
#'   \item mds - A data.frame containing the MDS embedding, with the number of rows equal to the number of samples.
#'   \item plot - A ggplot object containing the plot object. `print` of the object will create a plot.
#' }
#'
#' @examples
#' data(example_small_data)
#' sample_ids <- example_small_data$metadata$sample_id
#' # Run gloscope on first 10 PCA embeddings
#' # We use 'KNN' option for speed ('GMM' is slightly slower)
#' pca_embeddings <- example_small_data$pca_embeddings
#' pca_embeddings_subset <- pca_embeddings[,seq_len(10)] # select the first 10 PCs
#' dist_result <- gloscope(pca_embeddings_subset, sample_ids,
#'    dens="KNN",
#'    BPPARAM = BiocParallel::SerialParam(RNGseed=2))
#' # make a per-sample metadata
#' sample_metadata <- unique(example_small_data$metadata)
#' mds_result <- plotMDS(dist_mat = dist_result, metadata_df = sample_metadata ,
#' "sample_id", "phenotype",k=2)
#' mds_result$plot
#' head(mds_result$mds)
#' @import ggplot2
#' @importFrom MASS isoMDS
#' @export
#'
plotMDS <- function(dist_mat, metadata_df, sample_id, group_id, k=10){
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

  fit_df <- MASS::isoMDS(dist_mat, k = k, trace = FALSE)
  colnames(fit_df$points) <- paste0("Coordinate",seq_len(k))
  mds_df <- cbind(metadata_df, fit_df$points)

  colour_palette <- paletteBig()
  mds_plot <- ggplot2::ggplot(mds_df, aes(x = .data$Coordinate1, y = .data$Coordinate2,
    color = .data[[group_id]])) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values=colour_palette) +
    ggplot2::theme_bw()

  return(list(mds = mds_df, plot = mds_plot))
}
