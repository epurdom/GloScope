#' @title Plot the multidimensional scaling of the GloScope represenation
#'
#' @description This function calculates the multidimensional scaling for a
#'   GloScope divergence matrix and returns a ggplot object that plots it.
#'
#' @param dist_mat The divergence matrix output of `gloscope()`. Should be a
#'   symmetric, square matrix.
#' @param metadata_df A data frame contains each sample's metadata. Note this is
#'   NOT at the cell-level, and should have the same number of rows as dist_mat.
#' @param sample_id The column name or index in metadata_df that contains the
#'   sample ID. This is for ensuring alignment between the dist_mat and the
#'   metadata_df. The rownames of dist_mat are expected to match the sample_id
#'   values.
#' @param color_by The column name or index in metadata_df that should be used
#'   to color the points by. If missing all points will be the same color.
#' @param shape_by The column name or index in metadata_df that should be used
#'   to determine the shape of the points. If missing all points will be the
#'   same shape.
#' @param k Number of MDS dimension to generate, default = 10
#' @return A list containing the MDS embedding and plot of the distance matrix
#' \itemize{
#'   \item mds - A data.frame containing the MDS embedding, with the number of rows equal to the number of samples.
#'   \item plot - A ggplot object containing the plot object. `print` of the object will create a plot.
#' }
#'
#' @details The function calls \code{\link[MASS]{isoMDS}} from the MASS package,
#'   calculates the requested k coordinates of the MDS plot. It also creates a
#'   ggplot object that will plot the first two dimensions color or shape coded
#'   by the given variables in the metadata data frame.
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
#' mds_result <- plotMDS(dist_mat = dist_result, metadata_df = sample_metadata ,
#' sample_id="sample_id", color_by="phenotype",k=2)
#' head(mds_result$mds)
#' require(ggplot2)
#' mds_result$plot
#' # Add additional ggplot2 components to adapt figure
#' mds_result$plot + theme_bw()  + scale_color_manual(values=alpha(c("red","blue"),0.5))

#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_bw
#' @importFrom MASS isoMDS
#' @importFrom rlang .data
#' @export
#'
plotMDS <- function(dist_mat, metadata_df, sample_id, k=10, color_by, shape_by){
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

    fit_df <- MASS::isoMDS(dist_mat, k = k, trace = FALSE)
    colnames(fit_df$points) <- paste0("Coordinate",seq_len(k))
    mds_df <- cbind(metadata_df, fit_df$points)
    if(missing(color_by) & missing(shape_by)){
      mds_plot <- ggplot2::ggplot(mds_df, ggplot2::aes(x = .data$Coordinate1, y = .data$Coordinate2))
    }
    else if(missing(shape_by) & !missing(color_by)){
      mds_plot <- ggplot2::ggplot(mds_df, ggplot2::aes(x = .data$Coordinate1, y = .data$Coordinate2,
                                                       color = .data[[color_by]]))
      
    }
    else if(!missing(shape_by) & missing(color_by)){
      mds_plot <- ggplot2::ggplot(mds_df, ggplot2::aes(x = .data$Coordinate1, y = .data$Coordinate2,
                                                       shape = .data[[shape_by]]))
    } 
    else{
      mds_plot <- ggplot2::ggplot(mds_df, ggplot2::aes(x = .data$Coordinate1, y = .data$Coordinate2,
                                                       color = .data[[color_by]],shape = .data[[shape_by]]))
    }
    mds_plot <- mds_plot+ ggplot2::geom_point() 
    return(list(mds = mds_df, plot = mds_plot))
}
