#' @title calculate statistical divergence for each pair of samples
#'
#' @description This function calculates a matrisx of pairwise divergences between input GloScope representations.
#'
#' @param embedding_matrix: a matrix of latent embeddings with rows corresponding to cells and columns to dimensions
#' @param cell_sample_ids: a list of the samples IDs each cell comes from. Length must match the number of rows in `embedding_matrix`
#' @param dens: the density estimation. One of c("GMM","KNN")
#' @param dist_mat: distance metric to calculate the distance. One of c("KL","JS") 
#' @param r: number of monte-carlo simulations to generate
#' @param num_components: a vector of integers for the number of components to fit GMMS to, default is 1:9
#' @param k: number of k nearest negibhour for KNN density estimation, default k = 50.
#' @param BPPARAM: BiocParallel parameters, default is running in serial
#' @param prefit_density: a named list of `densityMclust` objects for each sample, default is NULL
#' @param return_density: return the GMM parameter list or not (if applicable), default is FALSE
#' @return A matrix containing the paitwise divergence or distance between all pairs of samples
#'
#' @examples
#' \dontrun{
#' data(example_data)
#' sample_ids <- seurat_object[[]]$sample
#' pca_embeddings <- seurat_object@reductions$pca@cell.embeddings
#' pca_embeddings_subset <- pca_embeddings[,1:10]
#' dist_result <- gloscope(pca_embeddings_subset, sample_ids,
#'                     BPPARAM = BiocParallel::SerialParam(RNGseed=2))
#' dist_result
#' }
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
#' @rdname gloscope
#' @export

gloscope <- function(embedding_matrix, cell_sample_ids, dens = "GMM", dist_mat = "KL",
		r = 10000, num_components = c(1:9), k=50,
		BPPARAM = BiocParallel::SerialParam(),
		fit_density = NULL, return_density = FALSE){

	# Input safety check
	if(length(cell_sample_ids)!=nrow(embedding_matrix)){
		stop("The number of cells in the embedding matrix does not match the number of sample labels.")
	}

	# Extract the unique sample IDs
	cell_sample_ids <- as.character(cell_sample_ids)
	unique_sample_ids <- unique(cell_sample_ids)
	names(unique_sample_ids) <- unique_sample_ids
	
	# Results may be unreliable for samples with less than 500 cells
	# We raise a warning denoting samples for which this is the case
	MIN_CELLS <- 500
	cells_per_sample <- sapply(unique_sample_ids, function(x){nrow(embedding_matrix[cell_sample_ids==x,])})
	if(sum(cells_per_sample < MIN_CELLS) > 0){
		small_samples <- names(cells_per_sample)[cells_per_sample < MIN_CELLS]
		warning(paste0("The following samples have fewer than ", MIN_CELLS,
			       " cells. This may lead to unreliable results.\n"),
				paste(shQuote(small_samples, type="cmd"), collapse=", "))
	}

	# If a sample has fewer cells than the `k` argument and `dens=="KNN"`,
	# its cells are ignored and its pairwise divergences are NA-valued.
	knn_na_values <- FALSE # Flag to add NA divergences to output matrix
	if(dens == "KNN" && (sum(cells_per_sample < k) > 0)){
		knn_na_values <- TRUE
		knn_withheld_samples <- names(cells_per_sample)[cells_per_sample < k]
		warning(paste0("The following samples have fewer than the minimum of ", k,
			       " cells required to run k-Nearest Neighbors. 
			       Divergence pairs including these samples have value NA in the divergence matrix.\n",
			       paste(shQuote(knn_withheld_samples, type="cmd"), collapse=", ")))
		# Remove cells from the affected units
		embedding_matrix <- embedding_matrix[!(cell_sample_ids %in% knn_withheld_samples),]
		cell_sample_ids <- cell_sample_ids[!(cell_sample_ids %in% knn_withheld_samples)]
		# Remove the associated sample IDs to avoid unnecessary computation
		unique_sample_ids <- unique(cell_sample_ids)
		names(unique_sample_ids) <- unique_sample_ids
	}

	# We create a list indexed by sample ID containing its embedding matrix
	sample_matrix_list <- lapply(unique_sample_ids,function(x){embedding_matrix[(cell_sample_ids==x),]})
	# Saved `mclust` densities can be used instead of running the package by setting the
	# `fit_density` optional argument to a list indexed by sample ID containing a fit `densityMclust` object
	if(is.null(fit_density)){
		mod_list <- calc_dens(sample_matrix_list, dens = dens, k = k, BPPARAM = BPPARAM, num_components = num_components)
	} else {
		mod_list <- fit_density
	}

	sample_pairs <- utils::combn(unique_sample_ids, 2)
	# Convert patient pairs to a list for BiocParallel::bplapply
	patient_pair_list <- lapply(seq_len(ncol(sample_pairs)), function(i) sample_pairs[,i])
	# IMPORTANT: There are additional algorithms for divergence estimation implemented in this package
	# which are not accessible from the `gloscope` function. The optional arguments
	# `ndim`, `varapp`, `epapp`, and `ep` must be manually set below. See `R/calc_dist.R` for their details.
	divergence_list <- BiocParallel::bplapply(patient_pair_list,
		function(w){ calc_dist(mod_list = mod_list, s1 = w[1], s2 = w[2],
				       df_list = df_list, dist_mat = dist_mat, dens = dens,
				       r = r, k = k,
				       ndim = 10, varapp = FALSE, epapp = FALSE, ep = NA)},BPPARAM=BPPARAM)

	divergence_vec <- unlist(divergence_list)
	# Convert pair-wise distances to a symmetric distance matrix
	divergence_matrix <- matrix(0, ncol = length(unique_sample_ids), nrow = length(unique_sample_ids))
	rownames(divergence_matrix) <- unique_sample_ids
	colnames(divergence_matrix) <- unique_sample_ids

	for (i in seq_len(ncol(sample_pairs))){
	  divergence_matrix[sample_pairs[1, i], sample_pairs[2, i]] <- divergence_vec[i]
	  divergence_matrix[sample_pairs[2, i], sample_pairs[1, i]] <- divergence_vec[i]
	}
	if(knn_na_values){
		full_samples <- c(unique_sample_ids ,knn_withheld_samples)
		num_total <- length(full_samples)	
		num_included = length(unique_sample_ids)
		num_withheld <- length(knn_withheld_samples)
		na_pad_col <- matrix(, nrow = num_included, ncol = num_withheld)
		na_pad_row <- matrix(, nrow = num_withheld, ncol = num_total)
		divergence_matrix <- cbind(divergence_matrix,na_pad_col)
		divergence_matrix <- rbind(divergence_matrix,na_pad_row)
		rownames(divergence_matrix) <- full_samples
		rownames(divergence_matrix) <- full_samples
	}

	if(dens == "GMM"){
		mod_list <- lapply(mod_list, function(x) x[c("data", "classification", "uncertainty", "density")] = NULL)
	} 
	if(return_density && dens == "GMM"){
		return(list(dist = divergence_matrix,modlist = mod_list))
	} else{
		return(divergence_matrix)
	}
}
