#' @title Plot a heatmap of the GloScope represenation
#'
#' @description This function creates a heatmap of the given GloScope divergence
#'   matrix.
#'
#' @param dist_mat The divergence matrix output of `gloscope()`. Should be a
#'   symmetric, square matrix.
#' @param metadata_df A data frame contains each sample's metadata. Note this is
#'   NOT at the cell-level, and should have the same number of rows as dist_mat.
#' @param sample_id The column name or index in metadata_df that contains the
#'   sample ID. This is for ensuring alignment between the dist_mat and the
#'   metadata_df. The rownames of dist_mat are expected to match the sample_id
#'   values.
#' @param color_by A vector of column names or indices in metadata_df that should be used
#'   to color/annotate the samples. 
#' @param which_side One of "columns","rows", or "both", indicating whether the
#'   annotation of the samples in `color_by` should be on the rows, columns, or
#'   on both.
#' @param ... parameters passed to \code{\link[pheatmap]{pheatmap}}.
#' @return Invisibly returns the output of \code{\link[pheatmap]{pheatmap}}
#'
#' @details The function is a wrapper function to \code{\link[pheatmap]{pheatmap}}.
#'   `color_by` is used to create subset of the `metadata_df` to pass to
#'   `annotation_col` (if `which_side="columns"`) or `annotation_row` (if
#'   `which_side="rows"`). If `which_side="both"`, then it is passed to both,
#'   and `annotation_names_row` argument is set to `FALSE`, suppressing labeling
#'   both the columns and rows (which user can thus not override). All other
#'   arguments to \code{\link[pheatmap]{pheatmap}} can be passed directly by the user
#' @seealso \code{\link[pheatmap]{pheatmap}}
#' @examples
#' data(example_SCE_small)
#' sample_ids <- SingleCellExperiment::colData(example_SCE_small)$sample_id
#' # Run gloscope on first 10 PCA embeddings
#' # We use 'KNN' option for speed ('GMM' is slightly slower)
#' pca_embeddings <- SingleCellExperiment::reducedDim(example_SCE_small,"PCA")
#' pca_embeddings_subset <- pca_embeddings[,seq_len(10)] # select the first 10 PCs
#' dist_result <- gloscope(pca_embeddings_subset, sample_ids,
#'    dens="KNN",
#'    BPPARAM = BiocParallel::SerialParam(RNGseed=2))
#' # make a per-sample metadata
#' sample_metadata <- as.data.frame(unique(SingleCellExperiment::colData(example_SCE_small)[,c(1,2)]))
#' plotHeatmap(dist_mat = dist_result, metadata_df = sample_metadata ,
#' sample_id="sample_id", color_by="phenotype")
#' # Pass additional options to pheatmap to control colors of groups
#' library(RColorBrewer)
#' plotHeatmap(dist_mat = dist_result, metadata_df = sample_metadata ,
#' sample_id="sample_id", color_by="phenotype", which_side="both", 
#' annotation_colors=list(phenotype = c(Covid = "magenta", Healthy = "white")), 
#' color = colorRampPalette(brewer.pal(9, "PuBuGn"))(100))

#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dist
#' @export
#'
plotHeatmap <- function(dist_mat, metadata_df, sample_id, color_by, 
                        which_side=c("columns","rows","both"),
                        annotation_colors = NA,...){
    which_side<-match.arg(which_side)
    if(nrow(dist_mat)!= nrow(metadata_df)){
    stop("Not consistent patient number. Make sure your
            distance matrix and meta info have the same patient number.")
  }
  if(length(intersect(rownames(dist_mat), metadata_df[, sample_id]))!= nrow(dist_mat)){
    stop("Not consistent patient IDs. Make sure your
            distance matrix and meta info have the same patients.")
  }
  if(!(identical(rownames(dist_mat), metadata_df[, sample_id]))){
    metadata_df <- metadata_df[match(rownames(dist_mat), metadata_df[,sample_id]),,drop=FALSE]
  }
  # Create a data frame for the annotations
  if(!missing(color_by)){
    rownames(metadata_df)<-metadata_df[,sample_id]
    annotcol <- metadata_df[,color_by,drop=FALSE]
    if(which_side=="columns"){
      invisible(pheatmap(
        mat = as.dist(dist_mat),
        annotation_col = annotcol,
        ...)
      )
    }
    else if(which_side=="rows"){
      invisible(pheatmap(
        mat = as.dist(dist_mat),
        annotation_row = annotcol,
        ...)
      )
    }
    else if(which_side=="both"){
      invisible(pheatmap(
        mat = as.dist(dist_mat),
        annotation_row = annotcol,
        annotation_col = annotcol,
        annotation_names_row=FALSE,
        ...)
      )
    }
  }
  else{
    invisible(pheatmap(mat = as.dist(dist_mat),...))
  }

}