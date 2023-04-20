#' calculate  KL for a pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param r number of monte-carlo simulations to generate
#' @param ep error term added to the KL divergence calculation
#' @param epapp whether to apply the error term
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param sample1 sample 1 index
#' @param sample2 sample 2 index
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param ndim number of dimension reduction to keep
#' @param varapp logic variable for using variational approximation or not
#' @return a numeric value of distance between sample1 and sample2's distribution.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @importFrom FNN KL.dist
#' @importFrom rags2ridges KLdiv
#' @rdname CalcDist
#' @export

# put into one rd file

# KL divergence
calc_kl <- function(mod_list, sample1, sample2, df_list, r,
			dens, k,
			ndim, varapp,epapp, ep){
	if(dens == "GMM"){
		mclust_mod1 <- mod_list[[sample1]]
		mclust_mod2 <- mod_list[[sample2]]
		s <- .sample_mclust(mclust_mod1, r=r)
		pi_1 <- mclust_mod1$parameters$pro
		pi_2 <- mclust_mod2$parameters$pro
		cov_1 <- mclust_mod1$parameters$variance$sigma
		cov_2 <- mclust_mod2$parameters$variance$sigma

		mu_1 <- mclust_mod1$parameters$mean
		mu_2 <- mclust_mod2$parameters$mean

		if(varapp) {
			kl <- .KLvar(pi_1, pi_2,mu_1, mu_2, cov_1, cov_2)
		} else {
			dens1 <- predict(mclust_mod1, s, what = "dens", logarithm = TRUE)
			dens2 <- predict(mclust_mod2, s, what = "dens", logarithm = TRUE)
			if(epapp) {
				kl <- sum(dens1 - (dens2+ep))/r
			} else {
				kl <- sum(dens1 - dens2) / r
			}
		}
	}else if(dens == "KNN"){
		kl <- KL.dist(as.matrix(mod_list[[sample1]]), as.matrix(mod_list[[sample2]]), k = k)[k]
	}

	return(kl)
}

#' calculate  EMD for a pair of samples
#'
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param sample1 sample 1 index
#' @param sample2 sample 2 index
#' @param ndim number of dimension reduction to keep
#' @return a numeric value of distance between sample1 and sample2's distribution.
#'
#' @importFrom transport transport
#' @importFrom psych tr
#' @rdname CalcDist
#' @export

calc_EMD <- function(mod_list, sample1, sample2, dens, ndim){
	if(dens == "KNN"){
		stop("kNN density estimation is incompatible with the EMD metric.")
	}
	x1 <- mod_list[[sample1]]$parameters
	x2 <- mod_list[[sample2]]$parameters
	dist <- matrix(NA, length(x1$pro), length(x2$pro))
	for(i in 1:length(x1$pro)){
		for(j in 1:length(x2$pro)){
			dist[i,j] <- 1/2*(log(max(colSums(x2$variance$sigma[,,j]))/max(colSums(x1$variance$sigma[,,i]))) -ndim +
					tr(solve(x2$variance$sigma[,,j])%*%x1$variance$sigma[,,i]) +
					t(x2$mean[,j] - x1$mean[,i])%*%solve(x2$variance$sigma[,,j])%*%(x2$mean[,j] - x1$mean[,i])) +
					1/2*(log(max(colSums(x1$variance$sigma[,,i]))/max(colSums(x2$variance$sigma[,,j]))) - ndim +
					tr(solve(x1$variance$sigma[,,i])%*%x2$variance$sigma[,,j]) +
					t(x1$mean[,i] - x2$mean[,j])%*%solve(x1$variance$sigma[,,i])%*%(x1$mean[,i] - x2$mean[,j]))
		}
	}
	colnames(dist) <- 1:ncol(dist)
	w1 <- x1$pro
	w2 <- x2$pro
	flow <- transport(w1, w2, dist, method='primaldual')

	EMD <- 0

	for(i in 1:nrow(flow)){
		EMD <- EMD +  dist[flow[i,1], flow[i,2]]*flow[i,3]
	}
	return(EMD)
}

#' Calculate JS divergence for a pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param r number of monte-carlo simulations to generate
#' @param ep error term added to the KL divergence calculation
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param sample1 sample 1 index
#' @param sample2 sample 2 index
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param ndim number of dimension reduction to keep
#' @return a numeric value of distance between sample1 and sample2's distribution.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @importFrom FNN KL.dist
#' @rdname CalcDist
#' @export

# TODO: Put into one rd file?

calc_JS <- function(mod_list, sample1, sample2, df_list, r, dens, k, ep, ndim){
	if(dens == "GMM"){
		mclust_mod1 <- mod_list[[sample1]]
		mclust_mod2 <- mod_list[[sample2]]
		s1 <- .sample_mclust(mclust_mod1, r=r)
		s2 <- .sample_mclust(mclust_mod2, r=r)
		dens1_1 <- predict(mclust_mod1, s1, what = "dens", logarithm = TRUE)
		dens2_1 <- predict(mclust_mod2, s1, what = "dens", logarithm = TRUE)
		mixture_1 <- log(1/2*exp(dens1_1) + 1/2*exp(dens2_1))
		dens1_2 <- predict(mclust_mod1, s2, what = "dens", logarithm = TRUE)
		dens2_2 <- predict(mclust_mod2, s2, what = "dens", logarithm = TRUE)
		mixture_2 <- log(1/2*exp(dens1_2) + 1/2*exp(dens2_2))

		js <- sum(dens1_1 - mixture_1)/(2*r) + sum(dens2_2 - mixture_2)/(2*r)
	}else if(dens == "KNN"){
		knn2_1 <- .knn_query(df_list, input = sample2, query = sample1, k = k)
		knn1 <- mod_list[[sample1]]
		knn1_2 <- .knn_query(df_list, input = sample1, query = sample2, k = k)
		knn2 <- mod_list[[sample2]]
		js <- 1/(2*length(knn1))*sum(log(2*length(knn2) * knn2_1^ndim/(length(knn2) * knn2_1^ndim +
			(length(knn1)-1) * knn1^ndim))) +
			1/(2*length(knn2))*sum(log(2*length(knn1) * knn1_2^ndim/(length(knn1) * knn1_2^ndim +
			(length(knn2)-1) * knn2^ndim)))
	}

	return(js)
}

#' Main function for computing distance from samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param r number of monte-carlo simulations to generate
#' @param ep error term added to the KL divergence calculation
#' @param epapp whether to apply the error term
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param s1 sample 1 name
#' @param s2 sample 2 name
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param ndim number of dimension reduction to keep
#' @param dist_mat which distance metric to use
#' @param varapp logic variable for using variational approximation or not
#' @return a numeric value of estimated symmatrised KL divergence between
#' sample1 and sample2's distribution.
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @importFrom transport transport
#' @importFrom psych tr
#' @importFrom rags2ridges KLdiv
#' @rdname CalcDist
#' @export

calc_dist <- function(mod_list, s1, s2, df_list, r,
		      dens, k, ep, ndim, dist_mat, varapp, epapp){
	if(dist_mat == "KL"){
		if(dens == "KNN"){
			mydist <- calc_kl(mod_list = mod_list, sample1 = s1, sample2 = s2, df_list = df_list,
					 dens = dens, k = k, ndim = ndim, r = r, varapp = varapp,
					 epapp = epapp, ep = ep)
		} else{
			mydist <- calc_kl(mod_list = mod_list, sample1 = s1, sample2 = s2, df_list = df_list,
						dens = dens, k = k, ndim = ndim, r = r, varapp = varapp, epapp = epapp, ep = ep) +
					calc_kl(mod_list, sample1 = s2, sample2 = s1, df_list = df_list,
						dens = dens, k = k, ndim = ndim, r = r, varapp = varapp, epapp = epapp, ep = ep)
		}
	} else if(dist_mat == "EMD"){
		mydist <- calc_EMD(mod_list = mod_list, sample1 = s1, sample2 = s2,
				  dens = dens, ndim = ndim)
	} else if(dist_mat == "JS"){
		mydist <- calc_JS(mod_list = mod_list, sample1 = s1, sample2 = s2, df_list = df_list,
				 dens = dens, k = k, ndim = ndim, r = r, ep = ep)
	}

	return(mydist)
}
