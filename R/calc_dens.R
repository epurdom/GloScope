#' @title calculate density for each samples based on dimension reduction embedding
#'
#' @description This function loads a vector of the sample names, with length of the
#'  number of samples whose density is to be estimated, a data frame x contains the info
#'  of cells, with row number equals the cell number, a data matrix or frame dim_redu
#'  contains dimension reduction embedding for each cell, which is used to calculate the
#'  density for each sample.
#'
#' @param x a data frame with row corresponding to single cell, contains the info of
#' each cell
#' @param sample_id column names in x that contains the sample ID
#' @param dim_redu_data a data matrix or frame contains the dimension reduction embedding
#' for each cell
#' @param dens method used to estimate density, options are GMM (Gaussian mixture model)
#' and KNN (K-nearest Neighbor)
#' @param k number of nearest neighbor, used when dens = "KNN".
#' @return A list of length number of samples, contains the estimated density for each
#' sample
#'
#'
#' @examples
#' data(example_data)
#' set.seed(1)
#' sample_name = as.character(unique(example_data[, "donor_label"]))
#' dim_redu_data = example_data[,str_detect(colnames(example_data), "PC")]
#' mod_list = calc_dens(sample_name, sample_id = "donor_label",
#'                      x = example_data, dim_redu_data, dens = "GMM")
#'
#'
#' @importFrom mclust densityMclust
#' @importFrom stringr str_detect
#' @importFrom RANN nn2
#' @export



calc_dens = function(sample_name, sample_id, x, dim_redu_data, dens = c("GMM", "KNN"),k = NULL){
  mod_list = list()
  for (s in sample_name){
    subset_idx <- x[, sample_id] == s
    subset_dt <- as.matrix(dim_redu[subset_idx, 1:ndim])

    if(den == "GMM"){
      mod_list[[s]] <- densityMclust(subset_dt, G=1:9)
    }else if(den == "KNN"){
      if(is.null(k)){
        stop(paste("Did not specify k for KNN density estimation."))
      }
      mod_list[[s]] <- nn2(dim_redu[subset_idx, ], query = dim_redu, k = k)
    }
  }
  return(mod_list)
}
