

#' calculate  KL for a pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param n number of monte-carlo simulations to generate for \hat{p}
#' @param ep error term added to the KL divergence calculation
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param s1 sample 1 name
#' @param s2 sample 2 name
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param ndim number of dimension reduction to keep
#' @return a numeric value of distance between sample1 and sample2's distribution.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @rdname CalcDist
#' @export

# put into one rd file

# KL divergence
calc_kl <- function(mod_list, sample1, sample2, df_list, n,
                    dens, k,
                    ep, ndim){

  if(dens == "GMM"){
    mclust_mod1 <- mod_list[[sample1]]
    mclust_mod2 <- mod_list[[sample2]]
    s <- .sample_mclust(mclust_mod1, n)
    dens1 <- predict(mclust_mod1, s, what = "dens")
    dens2 <- predict(mclust_mod2, s, what = "dens")
    kl <- sum(log(dens1/(dens2 + ep))) / n
  }else if(dens == "KNN"){

    knn2 <- .knn_query(df_list, sample2, sample1, k = k)
    knn1 <- mod_list[[sample1]]
    kl <- mean(log(knn2/knn1))*ndim +
      log(nrow(df_list[[sample2]])/(nrow(df_list[[sample1]])-1))
  }

  return(kl)
}

#' calculate symmetrised  KL for a pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param n number of monte-carlo simulations to generate for \hat{p}
#' @param ep error term added to the KL divergence calculation
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param s1 sample 1 name
#' @param s2 sample 2 name
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param ndim number of dimension reduction to keep
#' @return a numeric value of estimatedsymmatrised KL divergence between
#' sample1 and sample2's distribution.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @rdname CalcDist

#' @export
#'
# make symmatrised KL
calc_dist <- function(mod_list, s1, s2, df_list, n = 10000,
                      dens,k, ep = 1e-64, ndim ){

  sym_kl = calc_kl(mod_list = mod_list, sample1 = s1, sample2 = s2, df_list = df_list,
                   dens = dens, k = k, ndim = ndim, n = n, ep = ep) +
    calc_kl(mod_list, sample1 = s2, sample2 = s1, df_list = df_list,
            dens = dens, k = k, ndim = ndim, n = n, ep = ep)
  return(sym_kl)
}
