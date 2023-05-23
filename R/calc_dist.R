#' @title Main function for computing distance from samples
#'
#' @description This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param mod_list a list contains each samples' estimated density
#' @param s1 sample 1 name
#' @param s2 sample 2 name
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param dist_mat which distance metric to use
#' @param dens type of density to estimate for.
#' @param r number of monte-carlo simulations to generate
#' @param k number of k nearest negibhour for KNN density estimation, default k = 50.
#' @param ndim number of dimension reduction to keep; only applicable for JS distance, default = 10
#' @param varapp logic variable for using variational approximation or not, default = FALSE
#' @param epapp whether to apply the error term, default = FALSE
#' @param ep error term added to the KL divergence calculation, default = NA
#' @return a numeric value of estimated symmatrised KL divergence between
#' sample1 and sample2's distribution.
#' @examples
#' data(example_data)
#' library(stringr)
#' sample_name <- as.character(unique(example_data[, "patient_id"]))
#' example_data[,"patient_id"] <- as.character(example_data[,"patient_id"])
#' df_list <- split(example_data, example_data[,"patient_id"])
#' df_list <- lapply(df_list, function(y) y[,str_detect(colnames(y), "PC")])
#' df_list <- lapply(df_list, function(y) as.matrix(y[,1:10]))
#' mod_list <- .calc_dens(df_list, dens = "KNN", BPPARAM = BiocParallel::SerialParam())
#' all_combn <- combn(sample_name, 2)
#' patient_pair_list <- lapply(seq_len(ncol(all_combn)), function(i) all_combn[,i])
#' distance_list <- BiocParallel::bplapply(patient_pair_list, function(w){
#'                                         .calc_dist(mod_list = mod_list, df_list = df_list, k = 50,
#'                                         s1 = w[1], s2 = w[2], dens = "KNN", ndim = 10,
#'                                         r=10000, ep = 1e-64, dist_mat = "KL", varapp = FALSE,
#'                                         epapp = FALSE)},BPPARAM=BiocParallel::SerialParam())
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @rdname CalcDist

.calc_dist <- function(mod_list, s1, s2, df_list, dist_mat, dens, r, k,
                      ndim = 10, varapp = FALSE, epapp = FALSE, ep = NA){
  if(dist_mat == "KL"){
    if(dens == "KNN"){
      mydist <- .calc_kl(mod_list = mod_list, sample1 = s1, sample2 = s2, df_list = df_list,
                        dens = dens, k = k, r = r, varapp = varapp,
                        epapp = epapp, ep = ep)
    } else{
      mydist <- .calc_kl(mod_list = mod_list, sample1 = s1, sample2 = s2, df_list = df_list,
                        dens = dens, k = k, r = r, varapp = varapp, epapp = epapp, ep = ep) +
        .calc_kl(mod_list, sample1 = s2, sample2 = s1, df_list = df_list,
                dens = dens, k = k, r = r, varapp = varapp, epapp = epapp, ep = ep)
    }
  }  else if(dist_mat == "JS"){
    mydist <- .calc_JS(mod_list = mod_list, sample1 = s1, sample2 = s2, df_list = df_list,
                      dens = dens, k = k, ndim = ndim, r = r, ep = ep)
  }

  return(mydist)
}

#' @title calculate KL distance for a pair of samples
#'
#' @description This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param r number of monte-carlo simulations to generate
#' @param ep error term added to the KL divergence calculation
#' @param epapp whether to apply the error term
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param k number of k nearest negibhour for KNN density estimation, default k = 50.
#' @param sample1 sample 1 index
#' @param sample2 sample 2 index
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param varapp logic variable for using variational approximation or not
#' @return a numeric value of distance between sample1 and sample2's distribution.
#' @examples
#' data("example_data")
#' library(stringr)
#' sample_name <- as.character(unique(example_data[, "patient_id"]))
#' example_data[,"patient_id"] <- as.character(example_data[,"patient_id"])
#' df_list <- split(example_data, example_data[,"patient_id"])
#' df_list <- lapply(df_list, function(y) y[,str_detect(colnames(y), "PC")])
#' df_list <- lapply(df_list, function(y) as.matrix(y[,1:10]))
#' mod_list <- .calc_dens(df_list, dens = "KNN", BPPARAM = BiocParallel::SerialParam())
#' dist_test <- .calc_kl(mod_list, sample_name[1], sample_name[2], df_list, dens = "KNN",
#'                      r = 10000,k=50,varapp=FALSE,epapp=FALSE, ep = 1e-64)
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @importFrom FNN KL.dist
#' @rdname CalcDist

# KL divergence
.calc_kl <- function(mod_list, sample1, sample2, df_list, r,
			dens, k, epapp, ep){
	if(dens == "GMM"){
		mclust_mod1 <- mod_list[[sample1]]
		mclust_mod2 <- mod_list[[sample2]]
		s <- .sample_mclust(mclust_mod2, r=r)
		pi_1 <- mclust_mod1$parameters$pro
		pi_2 <- mclust_mod2$parameters$pro
		cov_1 <- mclust_mod1$parameters$variance$sigma
		cov_2 <- mclust_mod2$parameters$variance$sigma

		mu_1 <- mclust_mod1$parameters$mean
		mu_2 <- mclust_mod2$parameters$mean

    ## Old option for approximating KL; have disabled it.
    varapp<-FALSE
		if(varapp) {
#			kl <- .KLvar(pi_1, pi_2,mu_1, mu_2, cov_1, cov_2)
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



#' Calculate JS divergence for a pair of samples
#'
#' This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param r number of monte-carlo simulations to generate
#' @param ep error term added to the KL divergence calculation
#' @param mod_list a list contains each samples' estimated density
#' @param dens type of density to estimate for.
#' @param k number of k nearest negibhour for KNN density estimation, default k = 50.
#' @param sample1 sample 1 index
#' @param sample2 sample 2 index
#' @param df_list a list contain each samples' dimension reduction embedding
#' @param ndim number of dimension reduction to keep
#' @return a numeric value of distance between sample1 and sample2's distribution.
#' @examples
#' data(example_data)
#' library(stringr)
#' sample_name <- as.character(unique(example_data[, "patient_id"]))
#' example_data[,"patient_id"] <- as.character(example_data[,"patient_id"])
#' df_list <- split(example_data, example_data[,"patient_id"])
#' df_list <- lapply(df_list, function(y) y[,str_detect(colnames(y), "PC")])
#' df_list <- lapply(df_list, function(y) as.matrix(y[,1:10]))
#' #working with large data set, use BiocParallel
#' mod_list <- .calc_dens(df_list, dens = "KNN", BPPARAM = BiocParallel::SerialParam())
#' dist_test <- .calc_JS(mod_list, sample_name[1], sample_name[2], df_list, dens = "KNN",
#'                      r = 10000,k=50,ep = 1e-64,ndim = 10)
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom
#' @importFrom MASS mvrnorm
#' @importFrom RANN nn2
#' @importFrom FNN KL.dist
#' @rdname CalcDist

.calc_JS <- function(mod_list, sample1, sample2, df_list, r, dens, k, ep, ndim){
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


