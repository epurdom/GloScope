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
#' @rdname gloscope_proportion
#' @export

gloscope_proportion <- function(cell_sample_ids, cell_type_ids,
                                ep = 0, dist_mat = c("KL", "JS")){
  if(length(cell_sample_ids)!=length(cell_type_ids)){
    stop("Lengths of cell id and cell type are not equal!")
  }
  cell_sample_ids <- as.character(cell_sample_ids)
  cell_type_ids <- as.character(cell_type_ids)
  
  cluster_table <- table(cell_sample_ids, cell_type_ids)
  clusprop <- matrix(cluster_table, ncol = ncol(cluster_table), 
                     dimnames = dimnames(cluster_table))
  if(sum(clusprop==0)>0 & ep == 0){
    warning("There are elements haing 0 proportion! You may get invalid results. Please consider setting ep to be e.g 0.5.")
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
                            function(w){ .calc_prop(prop1 = clusprop[w[1],], 
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

# Helper function to obtain KL divergence using cell type proportion

#' @title Helper function to obtain the KL distance using cell type proportion
#'
#' @description The helper function `.calc_prop` calculates the KL divergence from
#'   two set of proportions.
#'
#' @param prop1 sample1's cluster proportion
#' @param prop2 sample2's cluster proportion
#' @param dist_mat distance metric to calculate the distance. One of
#'   c("KL","JS")
#' @return A single value contains the  KL divergence value calculated for the 
#'    2 samples' cluster proportion.
#'
#' @noRd
.calc_prop <- function(prop1, prop2, dist_mat){
  if(dist_mat == "KL"){
    KLdist <-  sum(prop1*(log(prop1) - log(prop2))) +
      sum(prop2*(log(prop2) - log(prop1)))
  }else if(dist_mat == "JS"){
    KLdist <-  1/2* sum(prop1*(log(prop1) - log(1/2*prop1 + 1/2*prop2))) +
      1/2* sum(prop2*(log(prop2) - log(1/2*prop1 + 1/2*prop2)))
  }
  return(KLdist)
}