#' @title plot MDS representation of the distance matrix
#'
#' @description This function loads the distance matrix,
#' and the info for each patients, including patient id
#' and patient condition group. It will create the MDS plot
#' for the distance matrix and color-coded each patient's dot
#' based on the condition.
#'
#' @param dist_mat the distance matrix of comparing each pairs
#' of patients' distributions
#' @param x a data frame contains each patients' info
#' @param sample_id column name in x that contains the sample ID
#' @param group_id column name in x that contains the patient condition
#' @return MDS plot of the distance matrix, color-coded by patient conditions
#'
#' @examples
#' data(example_data)
#' set.seed(1)
#' dist_mat <- distMat(example_data, "donor_label", "PC", 1:10)
#' pat_info <- unique(example_data[, c("donor_label", "joint_region_label"])])
#'
#' plotMDS(dist_mat = dist_mat, n = 10,
#'  x =  pat_info, "donor_label", "joint_region_label")
#'
#'
#' @import ggplot2
#' @importFrom MASS isoMDS
#' @importFrom clusterExperiment bigPalette
#' @export
#'
plotMDS = function(dist_mat, n=10, x, sample_id, group_id){
  if(nrow(dist_mat)!= nrow(meta)){
    stop(paste("Not consistent patient number. Make sure your
               distance matrix and meta info have the same patient
               number."))
  }
  if(length(intersect(rownames(dist_mat), x[, sample_id]))!= nrow(dist_mat)){
    stop(paste("Not consistent patient IDs. Make sure your
               distance matrix and meta info have the same patients."))
  }
  if(!(identical(rownames(dist_mat), x[, sample_id]))){
   x = x[match(rownames(dist_mat), x[,sample_id]),]
   }
    fit_df <- isoMDS(dist_mat, k = n, trace = FALSE)
    colnames(fit_df$points) <- paste0("Coordinate",1:n)
    mds_df <- cbind(x, fit_df$points)

  mds_plot = ggplot(mds_df, aes_string(x = "Coordinate1", y = "Coordinate2",
                                    color = group_id)) + geom_point() +
    scale_color_manual(values=bigPalette) +
    theme_bw()

  return(mds_plot)
}
