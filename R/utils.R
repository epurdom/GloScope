#' generate monte-carlo simulation from GMM estimated density for integrating expectation for KL divergence
#'
#' This function uses a object from mclust package as input, which contains the estimated density \hat{p} and
#' performs monte-carlo simulation based on \hat{p}
#'
#' @param mclust_mod mclust object which contains the estimated density \hat{p}
#' @param n number of simulations to generate. Default = 10,000.
#' @return A matrix contains the symmetrised KL divergence value calculated for each pair of samples.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#'

# ==============================================================================
# generate monte-carlo simulation from GMM estimated density \hat{p} for integrating
# expectation in KL divergence
# ------------------------------------------------------------------------------


.sample_mclust <- function(mclust_mod, n = 10000){
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
