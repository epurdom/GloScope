#' @title Calculate divergence between all sample pairs' cell type proportion
#'
#' @description This function calculates a matrix of pairwise divergences
#'   between input samples' cell type proportion.
#'
#' @param cell_sample_ids a vector of the samples IDs each cell comes from. Length
#'   must match the number of element in `cell_type_ids`
#' @param cell_type_ids a vector of use defined cell type
#' @param ep an numeric value added to the summary counts. Default ep = 0 means nothing will be added.
#' @param dist_metric distance metric to calculate the distance.
#' @return clusprop_dist a symmetric matrix of divergences 
#' @details Options for `dist_metric` are as follows: "KL" calculates the
#'   symmetric-KL divergence. "JS" calculates the Jenson-Shannon distance. "TV"
#'   calculates the total variation distance.
#' @examples
#' # Bring in small example data of single cell embeddings
#' data(example_SCE_small)
#' sample_id <- SingleCellExperiment::colData(example_SCE_small)$sample_id 
#' cluster_id <- SingleCellExperiment::colData(example_SCE_small)$cluster_id 
#' dist_result <- gloscopeProp(sample_id, cluster_id, ep = 0.5, 
#'                                    dist_metric = "KL")
#' dist_result
#' @importFrom utils combn
#' @rdname gloscopeProp
#' @export

gloscopeProp <- function(cell_sample_ids, cell_type_ids,
                                ep = 0, dist_metric = c("KL", "JS", "TV")){
  if(length(cell_sample_ids)!=length(cell_type_ids)){
    stop("Lengths of cell id and cell type are not equal!")
  }
  cell_sample_ids <- as.character(cell_sample_ids)
  cell_type_ids <- as.character(cell_type_ids)
  
  cluster_table <- table(cell_sample_ids, cell_type_ids)
  clusprop <- matrix(cluster_table, ncol = ncol(cluster_table), 
                     dimnames = dimnames(cluster_table))
  if(dist_metric %in% c("KL","JS")){
    if(sum(clusprop==0)>0 & ep == 0){
      warning("There are cell-types in some samples having 0 counts! You may get invalid results with this distance type. Please consider setting ep>0, e.g 0.5.")
    }else if (sum(clusprop==0)>0 & ep != 0){
      message(paste0("There are cell-types in some samples having 0 counts. An amount ", ep, " has been added to all counts."))
    }
    clusprop<- clusprop+ep
    
  }
  clusprop <- t(apply(clusprop, 1, function(x) x/sum(x)))
  
  unique_sample_ids <- rownames(clusprop)
  sample_pairs <- utils::combn(unique_sample_ids, 2)
  # Convert patient pairs to a list for BiocParallel::bplapply
  patient_pair_list <- lapply(seq_len(ncol(sample_pairs)), function(i) sample_pairs[,i])
  divergence_list <- lapply(patient_pair_list,
                            function(w){ .calc_prop(prop1 = clusprop[w[1],], 
                                                  prop2 = clusprop[w[2],], 
                                                  dist_metric = dist_metric)})
  
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
#' @param dist_metric distance metric to calculate the distance. One of
#'   c("KL","JS")
#' @return A single value contains the  KL divergence value calculated for the 
#'    2 samples' cluster proportion.
#'
#' @noRd
.calc_prop <- function(prop1, prop2, dist_metric){
  if(dist_metric == "KL"){
    out <-  sum(prop1*(log(prop1) - log(prop2))) +
      sum(prop2*(log(prop2) - log(prop1)))
  }else if(dist_metric == "JS"){
    out <-  1/2* sum(prop1*(log(prop1) - log(1/2*prop1 + 1/2*prop2))) +
      1/2* sum(prop2*(log(prop2) - log(1/2*prop1 + 1/2*prop2)))
  }
  else if(dist_metric == "TV"){
    out <-  1/2* sum(abs(prop1-prop2))
  }
  return(out)
}

