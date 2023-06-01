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
#' @importFrom MASS mvrnorm
#' @noRd

# Helper function to sample from a fit `mclust` GMM
.sample_mclust <- function(mclust_mod, r){
  p <- table(mclust_mod$classification)/mclust_mod$n
  z <- stats::rmultinom(1, size = r, prob = p)[,1]
  z <- z[z!=0] # clusters with zero samples raise an error in mvrnorm, so we drop them
  samples <- list()
  for ( i in 1:length(z)){
    samples[[i]] <- mvnfast::rmvn(z[i],
                            mu = mclust_mod$parameters$mean[, i],
                            sigma = mclust_mod$parameters$variance$sigma[, , i])
  }
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
#' @param k Number of k nearest neighbours for KNN density estimation
#' @return A matrix contains the symmetrised KL divergence value calculated for
#'   each pair of samples.
#'
#' @importFrom RANN nn2
#' @noRd
.knn_query = function(df_list, query, input, k){
  knnq_dist <- RANN::nn2(df_list[[input]], df_list[[query]], k = k)$nn.dists[,k]
  return(knnq_dist)
}

#' @description `paletteBig` is a small helper function to create a large color
#'   palette for plotting
#'
#' @rdname plotMDS
paletteBig <- function(){
	bigPalette <- c(
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
		'cadetblue3')
	return(bigPalette)
}
