#' generate monte-carlo simulation from GMM estimated density for integrating expectation for KL divergence
#'
#' This function uses a object from mclust package as input, which contains the estimated density \hat{p} and
#' performs monte-carlo simulation based on p hat.
#'
#' @param mclust_mod mclust object which contains the estimated density p hat.
#' @param n number of simulations to generate. Default = 10,000.
#' @return A matrix contains the symmetrised KL divergence value calculated for each pair of samples.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @noRd

# ==============================================================================
# generate monte-carlo simulation from GMM estimated density \hat{p} for integrating
# expectation in KL divergence
# ------------------------------------------------------------------------------


.sample_mclust <- function(mclust_mod, n){
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

.knn_query = function(df_list, query, input, k = k){
  knnq_dist = nn2(df_list[[input]], df_list[[query]], k = k)$nn.dists[,k]
  return(knnq_dist)
}

.KLvar <- function(pi_1, pi_2, mu_1,mu_2, cov_1, cov_2){
  Dvar <- 0
  for(i in 1:length(pi_1)){

    num <- don <- 0
    for(j in 1:length(pi_1)){
      num <- num + pi_1[j]*exp(-KLdiv(mu_1[,i], mu_1[,j], cov_1[,,i], cov_1[,,j]))
    }
    for(k in 1:length(pi_2)){
      don <- don + pi_2[k]*exp(-KLdiv(mu_1[,i], mu_2[,k], cov_1[,,i], cov_2[,,k]))
    }
    Dvar <- Dvar + pi_1[i]*(log(num) - log(don))
  }
  return(Dvar)
}


#' prepare for plotting colors

#' @rdname plottingColors
#' @export
bigPalette<-c(
	'#E31A1C',
	'#1F78B4',
	'#33A02C',
	'#FF7F00',
	'#6A3D9A',
	'#B15928',
	'#A6CEE3',
	'#bd18ea',
	'cyan',
	'#B2DF8A',
	'#FB9A99',
	"deeppink4",
	'#00B3FFFF',
	'#CAB2D6',
	'#FFFF99',
	'#05188a',
	'#CCFF00FF',
	'cornflowerblue',
	'#f4cc03',
	'black',
	'blueviolet',
	'#4d0776',
	'maroon3',
	'blue',
#	'grey',
	'#E5D8BD',
	'cadetblue4',
	'#e5a25a',
	"lightblue1",
	'#F781BF',
	'#FC8D62',
	'#8DA0CB',
	'#E78AC3',
	'green3',
	'#E7298A',
	'burlywood3',
	'#A6D854',
	"firebrick",
	'#FFFFCC',
	"mediumpurple",
	'#1B9E77',
	'#FFD92F',
	'deepskyblue4',
	"yellow3",
	'#00FFB2FF',
	'#FDBF6F',
	'#FDCDAC',
	"gold3",
	'#F4CAE4',
	'#E6F5C9',
	'#FF00E6FF',
	'#7570B3',
	"goldenrod",
	'#85848f',
	"lightpink3",
	"olivedrab",
#	"plum",
#	"lightskyblue3",
#	"mediumturquoise",
	'cadetblue3'
)

#' @importFrom grDevices colors
.rcolors<-function(){
	#so sure that setting seed doesn't mess up when installing package
	set.seed(23589)
	return(sample(colors()[-c(152:361)]))
}

#' @rdname plottingColors
#' @export
massivePalette<-unique(c(bigPalette,.rcolors()))
