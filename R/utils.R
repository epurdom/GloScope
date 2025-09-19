#' @title Draw samples from a fit GMM for Monte Carlo approximation
#'
#' @description The helper function `.sample_mclust` draws `r` samples from a
#'   GMM fit to a sample's latent embedding. This is used in `.calc_KL()` to
#'   compute a Monte Carlo approximation of the symmetric KL divergence between
#'   sample pairs.
#'
#' @param mclust_mod An `mclustDensity` object of a sample's estimated density
#' @param r Number of Monte Carolosimulations to generate, default = 10000.
#' @return A matrix contains the symmetrised KL divergence value calculated for
#'   each pair of samples.
#'
#' @importFrom stats rmultinom
#' @importFrom mvnfast rmvn
#' @noRd

# Helper function to sample from a fit `mclust` GMM
.sample_mclust <- function(mclust_mod, r){
    p <- table(mclust_mod$classification)/mclust_mod$n
    z <- stats::rmultinom(1, size = r, prob = p)[,1]
    non_zero <- which(z!=0) # clusters with zero samples raise a sampling error
    samples <- lapply(non_zero, function(i){
        mvnfast::rmvn(z[i], mu = mclust_mod$parameters$mean[, i],
        sigma = mclust_mod$parameters$variance$sigma[, , i])
    })
    samples <- do.call("rbind", samples)
    return(samples)
}

#' @title Helper function to query nearest neighbour distances
#'
#' @description This helper function `.knn_query` gets the distance from each
#'   cell in a query matrix to its k-th nearest neighbour in a target matrix.
#'   This is used to approximate various statistical divergences in
#'   `R/calc_dist.R`.
#'
#' @param df_list A named list with each sample's reduced dimension embedding
#' @param query The name or index of the sample whose query matrix is used to
#'   compute distances from
#' @param input The name or index of the sample whose matrix is used to find the
#'   k-th NN
#' @param KNN_params a list of arguments for the RANN::nn2 funciton, including `k`
#' @return A matrix contains the symmetrised KL divergence value calculated for
#'   each pair of samples.
#'
#' @importFrom RANN nn2
#' @noRd
.knn_query <- function(df_list, query, input, KNN_params){
    KNN_params$data <- df_list[[input]]
    KNN_params$query <- df_list[[query]]
    knnq_dist <- do.call(RANN::nn2,KNN_params)$nn.dists[,KNN_params$k]
    return(knnq_dist)
}



#' @title Helper function to pick number of GMM components
#'
#' @description `mclust::densityMclust` will raise an error if the 
#'   number of components vector `G` contains values at least as 
#'   large as the number of samples. If the user sets `num_components`
#'   with invalid entries, this function replaces them with n-1 and
#'   raises a warning.
#'
#' @param num_obs The number of observations to fit
#' @param num_components The user's requested vector of components
#' @return A valid vector of GMM number of components to fit
#'                                                                           
#' @noRd  
get_gmm_num_components_vec <- function(num_obs, num_components){
    if(any(num_components > (num_obs - 1))){
    warning(paste0("Unable to fit a GMM with ",
            paste(num_components[num_components > (num_obs - 1)], collapse = ", "),
        " components to a sample with ", num_obs, " cells.",
        " Replacing those parameters with ", num_obs - 1, " components instead."))
    }
    num_components[num_components > (num_obs - 1)] <- (num_obs - 1)
    return(unique(num_components)) # avoid duplicates
}

.as_full_matrix <- function(x) {
  m <- as.matrix(x)
  # ensure square, named, and symmetric indexability
  if (is.null(rownames(m)) && !is.null(colnames(m))) rownames(m) <- colnames(m)
  if (is.null(colnames(m)) && !is.null(rownames(m))) colnames(m) <- rownames(m)
  m
}

.testDistMeta<-function(dist_mat,metadata_df,sample_id){
  if(!c(sample_id) %in% names(metadata_df)) stop("sample_id does not define a variable in metadata_df")
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
  return(metadata_df)
  
}
