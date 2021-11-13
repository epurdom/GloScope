#' @title calculate symmetrised  KL for each pair of samples
#'
#' @description This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param x a data frame with row corresponding to single cell, contains a column of sample ID, and columns of dimension
#' reduction embedding
#' @param sample_id column names in x that contains the sample ID
#' @param dim_reduc dimension reduction index in x file (e.g. dim_reduc = "PC" for "PC_1").
#' @return A distance matrix contains the symmetrised KL divergence value calculated for each pair of samples.
#'
#' @examples
#' data(example_data)
#' set.seed(1)
#' dist_mat <- distMat(example_data, "donor_label", "PC", 1:10)
#'
#' #print out the distance matrix using PCA embedding.
#' dist_mat
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom stringr str_detect
#' @export







distMat = function(x, sample_id, dim_redu, ndim){
  sample_names = as.character(unique(x[, sample_id]))
  dim_redu_data = x[,str_detect(colnames(x), dim_redu)]

  mod_list = list()

  for (s in sample_names){
    subset_idx <- x[, sample_id] == s
    subset_dt <- as.matrix(dim_redu_data[subset_idx, 1:ndim])


    mod_list[[s]] <- densityMclust(subset_dt, G=1:9)
  }


  all_combn <- t(combn(sample_names, 2))
  dist_vec <- c()

  # calculate the distance
  for (i in 1:nrow(all_combn)){
    s1 <- all_combn[i, 1]
    s2 <- all_combn[i, 2]
    dist_vec <- c(dist_vec, calc_dist(mod_list[[s1]], mod_list[[s2]]))
  }

  dist_mat <- matrix(0, ncol = length(sample_names), nrow = length(sample_names))
  rownames(dist_mat) <- sample_names
  colnames(dist_mat) <- sample_names

  for (i in 1:nrow(all_combn)){
    dist_mat[all_combn[i, 1], all_combn[i, 2]] <- dist_vec[i]
    dist_mat[all_combn[i, 2], all_combn[i, 1]] <- dist_vec[i]
  }

  return(dist_mat)
}


