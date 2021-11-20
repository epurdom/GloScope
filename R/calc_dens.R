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
#' BiocParallel::register(BiocParallel::SerialParam())
#' sample_name = as.character(unique(example_data[, "donor_label"]))
#'
#' #working with large data set, use BiocParallel
#' mod_list = calc_dens(sample_id = "donor_label",
#'                      x = example_data, dim_redu = "PC", dens = "GMM",
#'                      BPPARAM=BiocParallel::SerialParam())
#'
#'
#' @importFrom mclust densityMclust
#' @importFrom stringr str_detect
#' @importFrom RANN nn2
#' @import BiocParallel
#' @export



calc_dens = function(df_list, dens, k,
                     BPPARAM=BiocParallel::bpparam()){


  if(dens == "GMM"){
    mod_list <- bplapply(df_list, function(z) densityMclust(z, G = 1:9),
                              BPPARAM=BPPARAM)
  }else if(dens == "KNN"){
    if(is.null(k)){
      stop(paste("Did not specify k for KNN density estimation."))
    }
    mod_list <- bplapply(df_list, function(z) nn2(z, query = z, k = k)$nn.dists[,k],
                         BPPARAM=BPPARAM)
  }

  return(mod_list)
}
