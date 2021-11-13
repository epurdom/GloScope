

#' calculate  KL for a pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param mclust_mod1 mclust object for sample1, contains estimated density \hat{p_1}
#' @param mclust_mod2 mclust object for sample2, contains estimated density \hat{p_2}
#' @param n number of monte-carlo simulations to generate for \hat{p}
#' @param ep error term added to the KL divergence calculation
#' @return a numeric value of distance between sample1 and sample2's distribution.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @export



# KL divergence
calc_kl <- function(mclust_mod1, mclust_mod2, n = 10000, ep = 1e-64){
  s <- .sample_mclust(mclust_mod1, n)
  dens1 <- predict(mclust_mod1, s, what = "dens")
  dens2 <- predict(mclust_mod2, s, what = "dens")
  kl <- sum(log(dens1/(dens2 + ep))) / n
  return(kl)
}

#' calculate symmetrised  KL for a pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param mclust_mod1 mclust object for sample1, contains estimated density \hat{p_1}
#' @param mclust_mod2 mclust object for sample2, contains estimated density \hat{p_2}
#' @param n number of monte-carlo simulations to generate for \hat{p}
#' @param ep error term added to the KL divergence calculation
#' @return a numeric value of distance between sample1 and sample2's distribution.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @export
#'
# make symmatrised KL
calc_dist <- function(mclust_mod1, mclust_mod2, n = 10000, ep = 1e-64){
  sym_kl = calc_kl(mclust_mod1, mclust_mod2, n, ep) + calc_kl(mclust_mod2, mclust_mod1, n, ep)
  return(sym_kl)
}
