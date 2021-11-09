#' calculate symmetrised  KL for each pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param meta metadata file
#' @param sample_id column names in metadata file that contains the sample ID
#' @param disease column names in metadata file that contains the disease group annotation
#' @param dim_reduc dimension reduction index in metadata file (e.g. dim_reduc = "PC_" for "PC_1").
#' @return A distance matrix contains the symmetrised KL divergence value calculated for each pair of samples.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom mvrnorm
#' @export



# monte-carlo simualtion
sample_mclust <- function(mclust_mod, n = 10000){
  p <- table(mclust_mod$classification)/mclust_mod$n
  z <- rmultinom(1, size = n, prob = p)[,1]
  samples <- list()
  for ( i in 1:length(z)){
    samples[[i]] <- mvrnorm(z[i],
                            mu = mclust_mod$parameters$mean[, i],
                            Sigma = mclust_mod$parameters$variance$sigma[, , i])
  }
  samples <- do.call("rbind", samples)
  return(samples)
}


# KL divergence
calc_kl <- function(mclust_mod1, mclust_mod2, n = 10000, ep = 1e-64){
  s <- sample_mclust(mclust_mod1, n)
  dens1 <- predict(mclust_mod1, s, what = "dens")
  dens2 <- predict(mclust_mod2, s, what = "dens")
  kl <- sum(log(dens1/(dens2 + ep))) / n
  return(kl)
}


# make symmatrised KL
calc_dist <- function(mclust_mod1, mclust_mod2, n = 10000){
  calc_kl(mclust_mod1, mclust_mod2, n) + calc_kl(mclust_mod2, mclust_mod1, n)
}




SKL = function(meta, sample_id, disease,dim_reduc){
  meta_sub <- unique(meta[, c(sample_id,disease)])
  meta_sub[, sample_id] = as.character(meta_sub[,sample_id])
  sample_names = meta_sub[, sample_id]

  mod_list = list()

  for (s in sample_names){
    subset_idx <- meta[, sample_id] == s
    subset_dt <- as.matrix(meta[subset_idx, paste0(dim_red, 1:10)])


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

  dist_mat <- matrix(0, ncol = nrow(meta_sub), nrow = nrow(meta_sub))
  rownames(dist_mat) <- sample_names
  colnames(dist_mat) <- sample_names

  for (i in 1:nrow(all_combn)){
    dist_mat[all_combn[i, 1], all_combn[i, 2]] <- dist_pca[i]
    dist_mat[all_combn[i, 2], all_combn[i, 1]] <- dist_pca[i]
  }

  return(dist_mat)
}


