#' @title Calculate statistical divergence between all sample pairs
#'
#' @description This function calculates a matrix of pairwise divergences
#'   between input samples of single cell data.
#'
#' @param embedding_matrix a matrix of latent embeddings with rows corresponding
#'   to cells and columns to dimensions
#' @param cell_sample_ids a list of the samples IDs each cell comes from. Length
#'   must match the number of rows in `embedding_matrix`
#' @param dens the density estimation. One of c("GMM","KNN", "Prop)
#' @param dist_mat distance metric to calculate the distance. One of
#'   c("KL","JS")
#' @param r number of Monte Carlo simulations to generate
#' @param num_components a vector of integers for the number of components to
#'   fit GMMS to, default is seq_len(9)
#' @param k number of nearest neighbours for KNN density estimation, default k =
#'   50.
#' @param BPPARAM BiocParallel parameters, default is running in serial. Set
#'   random seed with `RNGseed` argument
#' @param prefit_density a named list of pre-fit `densityMclust` objects for
#'   each sample, default is NULL
#' @param return_density return the GMM parameter list or not (if applicable),
#'   default is FALSE
#' @return A matrix containing the pairwise divergence or distance between all
#'   pairs of samples
#'
#' @examples
#' # Bring in small example data of single cell embeddings
#' data(example_SCE_small)
#' sample_ids <- SingleCellExperiment::colData(example_SCE_small)$sample_id 
#' pca_embeddings <- SingleCellExperiment::reducedDim(example_SCE_small,"PCA")
#' # Run gloscope on first 10 PCA embeddings
#' # We use 'KNN' option for speed ('GMM' is slightly slower)
#' pca_embeddings_subset <- pca_embeddings[,seq_len(10)] # select the first 10 PCs
#' dist_result <- gloscope(pca_embeddings_subset, sample_ids,
#'          dens="KNN", BPPARAM = BiocParallel::SerialParam(RNGseed=2))
#' dist_result
#'
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel SerialParam
#' @importFrom utils combn
#' @rdname gloscope
#' @export

gloscope <- function(embedding_matrix, cell_sample_ids,
                dens = c("GMM","KNN"), dist_mat = c("KL","JS"),
                r = 10000, num_components = seq_len(9), k = 50,
                BPPARAM = BiocParallel::SerialParam(),
                prefit_density = NULL, return_density = FALSE){
    dens<-match.arg(dens)
    dist_mat<-match.arg(dist_mat)
    # Input safety check
    if(length(cell_sample_ids)!=nrow(embedding_matrix)){
        stop("The number of cells in the embedding matrix does not match the number of sample labels.")
    }

    # Extract the unique sample IDs as characters
    cell_sample_ids <- as.character(cell_sample_ids)
    unique_sample_ids <- unique(cell_sample_ids)
    names(unique_sample_ids) <- unique_sample_ids

    # Results may be unreliable for samples with less than 500 cells
    # We raise a warning denoting samples for which this is the case
    MIN_CELLS <- 500
    cells_per_sample <- vapply(unique_sample_ids,
        function(x){nrow(embedding_matrix[cell_sample_ids==x,,drop=FALSE])},integer(1))
    if(sum(cells_per_sample < MIN_CELLS) > 0){
        small_samples <- names(cells_per_sample)[cells_per_sample < MIN_CELLS]
        samples_to_warn <- paste(shQuote(small_samples, type="cmd"), collapse=", ")
        warning("The following samples have fewer than ", MIN_CELLS,
                " cells. This may lead to unreliable results.\n",
                samples_to_warn)
    }

    # If a sample has fewer cells than the `k` argument and `dens=="KNN"`,
    # its cells are ignored and its pairwise divergences are NA-valued.
    knn_na_values <- FALSE # Flag to add NA divergences to output matrix
    if(dens == "KNN" && (sum(cells_per_sample < k) > 0)){
        knn_na_values <- TRUE
        knn_withheld_samples <- names(cells_per_sample)[cells_per_sample < k]
        samples_to_warn <- paste(shQuote(knn_withheld_samples, type="cmd"), collapse=", ")
        warning("The following samples have fewer than the minimum of ", k,
                " cells required to run k-Nearest Neighbors. ",
                "Divergence pairs including these samples have value NA in the divergence matrix.\n",
                samples_to_warn)
        # Remove cells from the affected units
        embedding_matrix <- embedding_matrix[!(cell_sample_ids %in% knn_withheld_samples),,drop=FALSE]
        cell_sample_ids <- cell_sample_ids[!(cell_sample_ids %in% knn_withheld_samples)]
        # Remove the associated sample IDs to avoid an error when running kNN
        unique_sample_ids <- unique(cell_sample_ids)
        names(unique_sample_ids) <- unique_sample_ids
    }

    # We create a list indexed by sample ID containing each sample's embedding
    # matrix
    sample_matrix_list <- lapply(unique_sample_ids,
        function(x){embedding_matrix[(cell_sample_ids==x),,drop=FALSE]})
    # Saved `mclust` densities can be used instead of running the package by
    # setting the `prefit_density` optional argument to a list indexed by sample ID
    # containing a fit `densityMclust` object for each
    if(is.null(prefit_density)){
        mod_list <- .calc_dens(sample_matrix_list, dens = dens, k = k, BPPARAM = BPPARAM, num_components = num_components)
    } else {
        mod_list <- prefit_density
    }

    sample_pairs <- utils::combn(unique_sample_ids, 2)
    # Convert patient pairs to a list for BiocParallel::bplapply
    patient_pair_list <- lapply(seq_len(ncol(sample_pairs)), function(i) sample_pairs[,i])
    # IMPORTANT: There are additional algorithms for divergence estimation
    # implemented in this package which are not accessible from the `gloscope`
    # function. The optional arguments `varapp`, `epapp`, and `ep` must be manually
    # set below. See `R/.calc_dist.R` for their details.
    divergence_list <- BiocParallel::bplapply(patient_pair_list,
        function(w){ .calc_dist(mod_list = mod_list, s1 = w[1], s2 = w[2],
            df_list = sample_matrix_list, dist_mat = dist_mat, dens = dens,
            r = r, k = k,
            varapp = FALSE, epapp = FALSE, ep = NA)},BPPARAM=BPPARAM)
    
    # Convert pair-wise distances to a symmetric distance matrix
    divergence_vec <- unlist(divergence_list)
    divergence_matrix <- matrix(0, ncol = length(unique_sample_ids),
                                nrow = length(unique_sample_ids))
    rownames(divergence_matrix) <- unique_sample_ids
    colnames(divergence_matrix) <- unique_sample_ids
    divergence_matrix[lower.tri(divergence_matrix, diag=FALSE)] <- divergence_vec
    divergence_matrix[upper.tri(divergence_matrix)] <- t(divergence_matrix)[upper.tri(divergence_matrix)]

    if(knn_na_values){
        # pad matrix with NA divergences if kNN density estimate cannot be run for
        # some units
        full_samples <- c(unique_sample_ids, knn_withheld_samples)
        num_total <- length(full_samples)
        num_included <- length(unique_sample_ids)
        num_withheld <- length(knn_withheld_samples)
        na_pad_col <- matrix(NA, nrow = num_included, ncol = num_withheld)
        na_pad_row <- matrix(NA, nrow = num_withheld, ncol = num_total)
        divergence_matrix <- cbind(divergence_matrix,na_pad_col)
        divergence_matrix <- rbind(divergence_matrix,na_pad_row)
        rownames(divergence_matrix) <- full_samples
        colnames(divergence_matrix) <- full_samples
    }

    if(dens == "GMM"){
        mod_list <- lapply(mod_list,
            function(x) x[-which(names(x) %in% c("data", "classification", 
                "uncertainty", "density"))])
    }
    if(return_density && dens == "GMM"){
        return(list(dist = divergence_matrix, modlist = mod_list))
    } else{
        return(divergence_matrix)
    }
}

#' @title Calculate divergence between all sample pairs' cell type proportion
#'
#' @description This function calculates a matrix of pairwise divergences
#'   between input samples' cell type proportion.
#'
#' @param cell_sample_ids a vector of the samples IDs each cell comes from. Length
#'   must match the number of element in `cell_type_ids`
#' @param cell_type_ids a vector of use defined cell type
#' @param ep an integer of error term added to 0 proportion. Default ep = 0.
#' @param dist_mat distance metric to calculate the distance. One of
#'   c("KL","JS")
#' @return clusprop_dist a symmetric matrix of divergences 
#' @examples
#' # Bring in small example data of single cell embeddings
#' data(example_SCE_small)
#' sample_id <- SingleCellExperiment::colData(example_SCE_small)$sample_id 
#' cluster_id <- SingleCellExperiment::colData(example_SCE_small)$cluster_id 
#' dist_result <- gloscope_proportion(sample_id, cluster_id, ep = 0.5, 
#'                                    dist_mat = "KL")
#' dist_result
#' @importFrom utils combn
#' @rdname cluster_distance
#' @export

gloscope_proportion <- function(cell_sample_ids, cell_type_ids, ep = 0, 
                             dist_mat = c("KL", "JS")){
    if(length(cell_sample_ids)!=length(cell_type_ids)){
        stop("Lengths of cell id and cell type are not equal!")
    }
    cell_sample_ids <- as.character(cell_sample_ids)
    cell_type_ids <- as.character(cell_type_ids)

    cluster_table <- table(cell_sample_ids, cell_type_ids)
    clusprop <- matrix(cluster_table, ncol = ncol(cluster_table), 
        dimnames = dimnames(cluster_table))
    if(sum(clusprop==0)>0 & ep == 0){
      warning("There are elements haing 0 proportion! You may get invalid results. Please consider setting ep to be e.g 0.5")
    }else if (sum(clusprop==0)>0 & ep != 0){
      warning(paste0("There are elements haing 0 proportion! ep has been set to be ", ep, "."))
    }
    clusprop[which(clusprop==0)] <- ep
    clusprop <- t(apply(clusprop, 1, function(x) x/sum(x)))

    unique_sample_ids <- rownames(clusprop)
    sample_pairs <- utils::combn(unique_sample_ids, 2)
    # Convert patient pairs to a list for BiocParallel::bplapply
    patient_pair_list <- lapply(seq_len(ncol(sample_pairs)), function(i) sample_pairs[,i])
    divergence_list <- lapply(patient_pair_list,
                                              function(w){ .clus_KL(prop1 = clusprop[w[1],], 
                                                                    prop2 = clusprop[w[2],], 
                                                                    dist_mat = dist_mat)})
    
    # Convert pair-wise distances to a symmetric distance matrix
    divergence_vec <- unlist(divergence_list)
    divergence_matrix <- matrix(0, ncol = length(unique_sample_ids),
                                nrow = length(unique_sample_ids))
    rownames(divergence_matrix) <- unique_sample_ids
    colnames(divergence_matrix) <- unique_sample_ids
    divergence_matrix[lower.tri(divergence_matrix, diag=FALSE)] <- divergence_vec
    divergence_matrix[upper.tri(divergence_matrix)] <- t(divergence_matrix)[upper.tri(divergence_matrix)]
    

    return(divergence_matrix)

}
