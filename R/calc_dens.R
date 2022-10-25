#' @title calculate density for each samples based on dimension reduction embedding
#'
#' @description This function loads a vector of the sample names, with length of the
#'  number of samples whose density is to be estimated, a data frame x contains the info
#'  of cells, with row number equals the cell number, a data matrix or frame dim_redu
#'  contains dimension reduction embedding for each cell, which is used to calculate the
#'  density for each sample.
#'
#' @param df_list: the list contains each samples' dimension reduction embedding
#' @param dens: method used to estimate density, options are GMM (Gaussian mixture model)
#' and KNN (K-nearest Neighbor)
#' @param k: number of k nearest negibhour for KNN density estimation, default k = 50.
#' @param num_components: a vector of integers for the number of components to fit GMMS to, default is 1:9
#' @param BPPARAM: BiocParallel parameters
#' @return mod_list: a list of length number of samples, contains the estimated density for each
#' sample
#'
#'
#' @examples
#' data("example_data")
#' library(stringr)
#' set.seed(1)
#' sample_name <- as.character(unique(example_data[, "patient_id"]))
#' example_data[,"patient_id"] <- as.character(example_data[,"patient_id"])
#' df_list <- split(example_data, example_data[,"patient_id"])
#' df_list <- lapply(df_list, function(y) y[,str_detect(colnames(y), "PC")])
#' df_list <- lapply(df_list, function(y) as.matrix(y[,1:10]))
#' #working with large data set, use BiocParallel
#' mod_list <- calc_dens(df_list, dens = "GMM", BPPARAM = BiocParallel::SerialParam())
#'
#'
#' @importFrom mclust densityMclust
#' @importFrom stringr str_detect
#' @importFrom RANN nn2
#' @import BiocParallel
#' @export



calc_dens = function(df_list, dens = "GMM", k = 50, num_components = c(1:9),
                     BPPARAM = BiocParallel::bpparam()){


  if(dens == "GMM"){
    mod_list <- BiocParallel::bplapply(df_list, function(z) densityMclust(z, G = num_components, verbose = F, plot = F),
                              BPPARAM=BPPARAM)
  }else if(dens == "KNN"){

    mod_list <- df_list
  }

  return(mod_list)
}
