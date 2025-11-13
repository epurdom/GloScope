#' @title Calculate statistical divergence between all sample pairs
#'
#' @description This function calculates a matrix of pairwise divergences
#'   between input samples of single cell data.
#'
#' @param embedding_matrix a matrix or data.frame of latent embeddings with rows
#'   corresponding to cells and columns to dimensions
#' @param cell_sample_ids a vector of the samples IDs each cell comes from. Length
#'   must match the number of rows in `embedding_matrix`
#' @param dens the density estimation. One of c("GMM","KNN")
#' @param dist_metric distance metric to calculate the distance. One of
#'   c("KL","JS")
#' @param r number of Monte Carlo simulations to generate
#' @param num_components a vector of integers for the number of components to
#'   fit GMMs to, default is c(5,10,15,20)
#' @param k number of nearest neighbours for KNN density estimation, default k =
#'   50.
#' @param GMM_params optional mclust parameters, default is to restrict the fit
#'   model to only VVE
#' @param KNN_params optional arguments for either `FNN:KL.dist` (KL) or
#'   `RANN::nn2` (JS), default is NULL
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
                dens = c("GMM","KNN"), dist_metric = c("KL","JS"),
                r = 10000, num_components = c(5,10,15,20), k = 50,
                GMM_params = list(modelNames = c("VVE"),verbose=FALSE,plot=FALSE),
                KNN_params = NULL,
                BPPARAM = BiocParallel::SerialParam(),
                prefit_density = NULL, return_density = FALSE){
    dens<-match.arg(dens)
    dist_metric<-match.arg(dist_metric)
    # Input safety check
    if(length(cell_sample_ids)!=nrow(embedding_matrix)){
        stop("The number of cells in the embedding matrix does not match the number of sample labels.")
    }

    # Extract the unique sample IDs as characters
    cell_sample_ids <- as.character(cell_sample_ids)
    unique_sample_ids <- unique(cell_sample_ids)
    names(unique_sample_ids) <- unique_sample_ids

    if(is.null(colnames(embedding_matrix))){
      # GMM fitting now requires column names for matrix
      colnames(embedding_matrix)<-paste("Dim",seq_len(ncol(embedding_matrix)),sep="")
    }
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
        mod_list <- .calc_dens(sample_matrix_list, dens = dens, 
        k = k, num_components = num_components,
        GMM_params = GMM_params, BPPARAM = BPPARAM)
    } else {
        # Validate that prefit_density contains densityMclust objects
        for(i in seq_along(prefit_density)){
            obj <- prefit_density[[i]]
            obj_name <- names(prefit_density)[i]
            if(is.null(obj_name)) obj_name <- paste0("[[", i, "]]")
           
            # Check if it's a densityMclust object
            if(!inherits(obj, "densityMclust")){
                # Provide specific guidance if it's a Mclust object
                if(inherits(obj, "Mclust")){
                    stop(
                        sprintf(
                            "prefit_density%s is a 'Mclust' object, but 'densityMclust' is required.\n",
                            obj_name
                        ),
                        "  Use mclust::densityMclust() instead of mclust::Mclust().\n",
                        "  See ?gloscope for more details."
                    )
                } else {
                    stop(
                        sprintf(
                            "prefit_density%s must be a 'densityMclust' object.\n",
                            obj_name
                        ),
                        "  Found object of class: ", class(obj)[1], "\n",
                        "  Use mclust::densityMclust() to create the correct object type.\n",
                    )
                }
            }
        }
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
            df_list = sample_matrix_list, dist_metric = dist_metric, dens = dens,
            r = r, k = k, KNN_params = KNN_params,
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
            function(x){x[-which(names(x) %in% c("data", "classification","uncertainty", "density"))]})
    }
    if(return_density && dens == "GMM"){
        return(list(dist = divergence_matrix, modlist = mod_list))
    } else{
        return(divergence_matrix)
    }
}


