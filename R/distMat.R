#' @title calculate symmetrised  KL for each pair of samples
#'
#' @description This function loads the metadata, which contains the sample id, disease group, dimension reduction embedding which has names
#' in the form "dim_i". Dimension reduction embedding is used to calculate the density for each sample, denoted by their sample id.
#'
#' @param x a data frame with row corresponding to single cell, contains a column of sample ID, and columns of dimension
#' reduction embedding
#' @param sample_id column names in x that contains the sample ID
#' @param dim_redu dimension reduction index in x file (e.g. dim_reduc = "PC" for "PC_1").
#' @param n number of monte-carlo simulations to generate
#' @param ep error term added to the KL divergence calculation
#' @param epapp whether to apply the error term
#' @param dens type of density to estimate for.
#' @param num_components: a vector of integers for the number of components to fit GMMS to, default is 1:9
#' @param ndim number of dimension reduction to keep
#' @param BPPARAM BiocParallel parameters; NULL to let system pick
#' @param requested_cores if NULL BPPARAM, the number of requested cores
#' @param k number of k nearest negibhour for KNN density estimation, default k = 50.
#' @param dist_mat distance metric to calculate the distance
#' @param varapp logic variable for using variational approximation or not
#' @param fit_density a named list containing the fit density for each patient
#' @param returndens return the GMM parameter list or not
#' @return A distance matrix contains the symmetrised KL divergence value calculated for each pair of samples.
#'
#' @examples
#' data("example_data")
#' set.seed(1)
#' dist_result <- distMat(example_data, sample_id = "patient_id", dim_redu = "PC",
#'                     ndim = 10, dens = "KNN", n=10000, ep = 1e-64, dist_mat = "KL",
#'                     BPPARAM = BiocParallel::SerialParam(), varapp = FALSE,
#'                     returndens = FALSE, epapp = FALSE)
#'
#' #print out the distance matrix using PCA embedding.
#' dist_result
#'
#' @importFrom mclust densityMclust
#' @importFrom stats rmultinom predict
#' @importFrom utils combn
#' @importFrom MASS mvrnorm
#' @importFrom stringr str_detect
#' @importFrom RANN nn2
#' @importFrom FNN KL.dist
#' @importFrom rags2ridges KLdiv
#' @import BiocParallel
#' @rdname CalcDist
#' @export

distMat = function(x, sample_id, dim_redu, ndim, k=50 , dens = "GMM",
                   n = 10000,ep = 1e-64, dist_mat = "KL", num_components = c(1:9),
		   BPPARAM=NULL, requested_cores=1,
                   varapp = FALSE, returndens = FALSE, epapp = FALSE,
		   fit_density=NULL){
	# check available cores for parallelzation unless user has specified BPPARAM
	BPPARAM <- set_BPPARAM(BPPARAM,request_cores)

	if(is.null(fit_density)){
		sample_names <- as.character(unique(x[, sample_id]))
	} else{
		sample_names <- as.character(names(fit_density))
	}
	x[,sample_id] = as.character(x[,sample_id])

	df_list = split(x, x[,sample_id])
	df_list = lapply(df_list, function(y) y[,str_detect(colnames(y), dim_redu)])
	df_list = lapply(df_list, function(y) as.matrix(y[,1:ndim]))

	if(is.null(fit_density)){
		mod_list <- calc_dens(df_list, dens = dens, k = k, BPPARAM = BPPARAM)
	} else {
		mod_list <- fit_density
	}

	all_combn <- t(combn(sample_names, 2))

	dist_vec <- c()

  # calculate the distance
  for (i in 1:nrow(all_combn)){
    s1 <- all_combn[i, 1]
    s2 <- all_combn[i, 2]
    dist_vec <- c(dist_vec, calc_dist(mod_list = mod_list, df_list = df_list, k = k,
                                      s1 = s1, s2 = s2, dens = dens, ndim = ndim,
                                      n=n, ep = ep, dist_mat = dist_mat, varapp = varapp,
                                      epapp = epapp))
  }

  dist_mat <- matrix(0, ncol = length(sample_names), nrow = length(sample_names))
  rownames(dist_mat) <- sample_names
  colnames(dist_mat) <- sample_names

  for (i in 1:nrow(all_combn)){
    dist_mat[all_combn[i, 1], all_combn[i, 2]] <- dist_vec[i]
    dist_mat[all_combn[i, 2], all_combn[i, 1]] <- dist_vec[i]
  }
  if(dens == "GMM"){
    mod_list = lapply(mod_list, function(x) x[c("data", "classification", "uncertainty", "density")] = NULL)
  }
  if(returndens){
    return(list(dist = dist_mat,
                modlist = mod_list))
  }else{
    return(dist_mat)
    }
}


