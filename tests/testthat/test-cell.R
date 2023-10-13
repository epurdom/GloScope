library(testthat)
test_that("gloscope warnings for small cell counts", {
	expect_error(gloscope(subsample_data_subset, c(subsample_metadata$sample_id,NA)),
		regexp="The number of cells in the embedding matrix does not match the number of sample labels.",fixed=TRUE)
	sample_size_warnings <- capture_warnings(na_div_matrix <- gloscope(undersized_data_subset,undersized_metadata$sample_id,dens="KNN"))
	expect_match(sample_size_warnings,"^The following samples have fewer than 500 cells.*",all = FALSE)
	expect_match(sample_size_warnings,"^The following samples have fewer than the minimum of 50 cells*",all = FALSE)
	# confirm the NA divergences are added last for samples with less than k cells
	expect_equal(unname(colnames(na_div_matrix))[4],c("CV0234"))
	expect_equal(unname(rownames(na_div_matrix))[4],c("CV0234"))
	expect_equal(all(is.na(na_div_matrix[-4,-4])),FALSE)
	expect_equal(all(is.na(na_div_matrix["CV0234",])),TRUE)
	expect_equal(all(is.na(na_div_matrix[,"CV0234"])),TRUE)
	
	# confirm a warning when the user requests more GMM components than cells
	# Case when there are some valid component numbers specified
	num_component_warnings1 <- capture_warnings(updated_clusters <- get_gmm_num_components_vec(49,c(1, 49)))
	expect_match(num_component_warnings1,"^Unable to fit a GMM with 49 components to a sample with 49 cells. Replacing those parameters with 48 components instead.*",all = FALSE)
	expect_equal(updated_clusters,c(1,48))
	# Case when no valid component numbers specified
	num_component_warnings2 <- capture_warnings(small_gloscope <- gloscope(undersized_data_subset,undersized_metadata$sample_id,dens="GMM",num_components=c(50)))
	expect_match(num_component_warnings2,"^Unable to fit a GMM with 50 components to a sample with 49 cells. Replacing those parameters with 48 components instead.*",all = FALSE)
	expect_equal(dim(small_gloscope),c(4,4)) # Confirm a matrix is returned
})

test_that("gloscope works with KNN",{
  expect_silent(temp_knn<-gloscope(subsample_data_subset,subsample_metadata$sample_id,
      dens = "KNN", dist_mat = "KL",BPPARAM = BiocParallel::SerialParam(RNGseed=2)))
  #test dimensions
  expect_equal(dim(temp_knn),c(4,4))
  #test distances the same as in the past
  knn_expected_answer<-c(2.369535,4.478759,2.215369,3.203009,3.225507,1.520805) #answer got the first time
  expect_equal(round(temp_knn[upper.tri(temp_knn)],6),knn_expected_answer)
  #test diag zero
  expect_equal(unname(diag(temp_knn)),rep(0,4))
  #test row names/colnames
  expect_equal(unname(sort(colnames(temp_knn))),sort(as.character(unique(subsample_metadata$sample_id))))
  expect_equal(unname(sort(rownames(temp_knn))),sort(as.character(unique(subsample_metadata$sample_id))))
  # test symmetry
  expect_equal(isSymmetric(temp_knn),TRUE)
})

test_that("gloscope works with GMM",{
  expect_silent(temp_gmm<-gloscope(subsample_data_subset,subsample_metadata$sample_id,
      dens = "GMM", dist_mat = "KL",BPPARAM = BiocParallel::SerialParam(RNGseed=1)))
  #test dimensions
  expect_equal(dim(temp_gmm),c(4,4))
  #test distances the same as in the past
  gmm_expected_answer <- c(14.220048,13.089617,18.073184,10.097857,11.156613,10.425697) #answer got the first time
  expect_equal(round(temp_gmm[upper.tri(temp_gmm)],6),gmm_expected_answer)
  #test diag zero
  expect_equal(unname(diag(temp_gmm)),rep(0,4))
  #test row names/colnames
  expect_equal(unname(sort(colnames(temp_gmm))),sort(as.character(unique(subsample_metadata$sample_id))))
  expect_equal(unname(sort(rownames(temp_gmm))),sort(as.character(unique(subsample_metadata$sample_id))))
  # test symmetry
  expect_equal(isSymmetric(temp_gmm),TRUE)
})

test_that("different random seeds give different GMM results",{
  expect_false(all(gloscope(subsample_data_subset,subsample_metadata$sample_id,
                                   dens = "GMM", dist_mat = "KL",BPPARAM = BiocParallel::SerialParam(RNGseed=2)) ==
                 gloscope(subsample_data_subset,subsample_metadata$sample_id,
                          dens = "GMM", dist_mat = "KL",BPPARAM = BiocParallel::SerialParam(RNGseed=1))))
})

test_that("plotMDS works with output",{
  expect_silent(dist_mat <- gloscope(embedding_matrix=subsample_data_subset, cell_sample_ids=subsample_metadata$sample_id,dens="KNN"))
  pat_info <- unique(subsample_metadata[,c(1,2)])
  expect_silent(mds_result <- plotMDS(dist_mat = dist_mat,
    metadata_df=pat_info, sample_id="sample_id", group_id="phenotype", k = 2))
})

test_that("GMM density fitting returns `densityMclust` objects",{
  sample_ids <- subsample_metadata$sample_id
  # the following `lapply` creates the necessary input data structure for this fn.
  embeddings_list <- lapply(unique(sample_ids),function(x){subsample_data_subset[(sample_ids==x),]})
  names(embeddings_list) <- unique(sample_ids)
  gmm_density_list <- .calc_dens(embeddings_list, dens = "GMM", BPPARAM = BiocParallel::SerialParam(RNGseed = 2))
  # confirm that each element is in the output is a `densityMclust` object
  expect_true(all(unlist(lapply(gmm_density_list,function(x){class(x)[1]=="densityMclust"}))))
})

test_that("JS divergences are properly implemented",{
  # From Toy Examples 1 and 2 from "On Accuracy of PDF Divergence Estimators
  # and Their Applicability to Representative Data Sampling" by Budka et al. (2011)

  # Example 1
  set.seed(2)
  s1 <- mvnfast::rmvn(10000,mu=c(0,0),sigma=diag(2))
  s2 <- mvnfast::rmvn(10000,mu=c(0.5,-0.5),sigma= matrix(c(0.5,0.1,0.1,0.3),2,2))
  df_list_1 <- list(s1,s2)

  js_knn_1_expected <- 0.17
  mod_list_knn_1 <- .calc_dens(df_list_1, dens="KNN")
  js_knn_1 <- .calc_JS (mod_list_knn_1, df_list_1, 1, 2, dens = "KNN")
  expect_equal(round(js_knn_1,2),js_knn_1_expected)

  js_gmm_1_expected <- 0.18
  mod_list_gmm_1 <- .calc_dens(df_list_1, dens="GMM",
                               BPPARAM = BiocParallel::SerialParam(RNGseed = 2))
  js_gmm_1 <- .calc_JS (mod_list_gmm_1, df_list_1, 1, 2, dens = "GMM")
  expect_equal(round(js_gmm_1,2),js_gmm_1_expected)

  # Example 2
  set.seed(2)
  s1 <- mvnfast::rmvn(10000,mu=c(0,0),sigma= matrix(c(1,0,0,0.1),2,2))
  s2 <- mvnfast::rmvn(10000,mu=c(0,0),sigma= matrix(c(0.1,0,0,1),2,2))
  df_list_2 <- list(s1,s2)

  js_knn_2_expected <- 0.32
  mod_list_knn_2 <- .calc_dens(df_list_1, dens="KNN")
  js_knn_2 <- .calc_JS (mod_list_knn_2, df_list_2, 1, 2, dens = "KNN")
  expect_equal(round(js_knn_2,2),js_knn_2_expected)

  js_gmm_2_expected <- 0.33
  mod_list_gmm_2 <- .calc_dens(df_list_2, dens="GMM",
                               BPPARAM = BiocParallel::SerialParam(RNGseed = 2))
  js_gmm_2 <- .calc_JS (mod_list_gmm_2, df_list_2, 1, 2, dens = "GMM")
  expect_equal(round(js_gmm_2,2),js_gmm_2_expected)
})

test_that("the sKL divergences are properly implemented",{
  # From Toy Examples 1 and 2 from "On Accuracy of PDF Divergence Estimators
  # and Their Applicability to Representative Data Sampling" by Budka et al. (2011)

  # Example 1
  set.seed(2)
  s1 <- mvnfast::rmvn(10000,mu=c(0,0),sigma=diag(2))
  s2 <- mvnfast::rmvn(10000,mu=c(0.5,-0.5),sigma= matrix(c(0.5,0.1,0.1,0.3),2,2))
  df_list_1 <- list(s1,s2)

  kl_knn_1_expected <- 1.45
  mod_list_knn_1 <- .calc_dens(df_list_1, dens="KNN")
  kl_knn_1 <- .calc_kl (mod_list_knn_1, df_list_1, 1, 2, dens = "KNN")
  expect_equal(round(kl_knn_1,2),kl_knn_1_expected)

  kl_gmm_1_expected <- 1.84
  mod_list_gmm_1 <- .calc_dens(df_list_1, dens="GMM",
                               BPPARAM = BiocParallel::SerialParam(RNGseed = 2))
  kl_gmm_1 <- .calc_kl (mod_list_gmm_1, df_list_1, 1, 2, dens = "GMM")
  expect_equal(round(kl_gmm_1,2),kl_gmm_1_expected)

  # Example 2
  set.seed(2)
  s1 <- mvnfast::rmvn(10000,mu=c(0,0),sigma= matrix(c(1,0,0,0.1),2,2))
  s2 <- mvnfast::rmvn(10000,mu=c(0,0),sigma= matrix(c(0.1,0,0,1),2,2))
  df_list_2 <- list(s1,s2)

  kl_knn_2_expected <- 3.28
  mod_list_knn_2 <- .calc_dens(df_list_2, dens="KNN")
  kl_knn_2 <- .calc_kl (mod_list_knn_2, df_list_2, 1, 2, dens = "KNN")
  expect_equal(round(kl_knn_2,2),kl_knn_2_expected)

  kl_gmm_2_expected <- 4.05
  mod_list_gmm_2 <- .calc_dens(df_list_2, dens="GMM",
                               BPPARAM = BiocParallel::SerialParam(RNGseed = 2))
  kl_gmm_2 <- .calc_kl (mod_list_gmm_2, df_list_2, 1, 2, dens = "GMM")
  expect_equal(round(kl_gmm_2,2),kl_gmm_2_expected)
})

test_that("Divergences are properly computed with GloScope inputs and KNN",{
  sample_ids <- subsample_metadata$sample_id
  embeddings_list <- lapply(unique(sample_ids),function(x){subsample_data_subset[(sample_ids==x),]})
  names(embeddings_list) <- unique(sample_ids)
  knn_density_list <- .calc_dens(embeddings_list, dens = "KNN")
  sample_pairs <- utils::combn(unique(sample_ids), 2)
  patient_pair_list <- lapply(seq_len(ncol(sample_pairs)), function(i) sample_pairs[,i])
  w <- as.character(patient_pair_list[[6]]) # pick one unit
  w_div_kl_expected <- 1.52
  w_div_kl <- .calc_kl(knn_density_list, embeddings_list, w[1], w[2], dens = "KNN")
  w_div_js_expected <- 0.04
  w_div_js <- .calc_JS(knn_density_list, embeddings_list, w[1], w[2], dens = "KNN")
  expect_equal(round(w_div_kl,2),w_div_kl_expected)
  expect_equal(round(w_div_js,2),w_div_js_expected)
})


test_that("Divergences are properly computed with GloScope inputs and GMM",{
  sample_ids <- subsample_metadata$sample_id
  embeddings_list <- lapply(unique(sample_ids),function(x){subsample_data_subset[(sample_ids==x),]})
  names(embeddings_list) <- unique(sample_ids)
  gmm_density_list <- .calc_dens(embeddings_list, dens = "GMM", BPPARAM = BiocParallel::SerialParam(RNGseed = 2))
  sample_pairs <- utils::combn(unique(sample_ids), 2)
  patient_pair_list <- lapply(seq_len(ncol(sample_pairs)), function(i) sample_pairs[,i])
  w <- as.character(patient_pair_list[[2]]) # pick one unit
  w_div_kl_expected <- 5.06
  set.seed(2)
  w_div_kl <- .calc_kl(gmm_density_list, embeddings_list, w[1], w[2], dens = "GMM")
  w_div_js_expected <- 0.50
  set.seed(2)
  w_div_js <- .calc_JS(gmm_density_list, embeddings_list, w[1], w[2], dens = "GMM")
  expect_equal(round(w_div_kl,2),w_div_kl_expected)
  expect_equal(round(w_div_js,2),w_div_js_expected)
})

test_that("Divergences using cell type works properly",{
  sample_ids <- subsample_metadata$sample_id
  celltype <- subsample_metadata$cluster_id
  zero_prop_warnings <- capture_warnings(inf_div_matrix <- gloscope_proportion(sample_ids,celltype, dist_mat= "KL"))
  expect_match(zero_prop_warnings,"There are elements haing 0 proportion! You may get invalid results. Please consider setting ep to be e.g 0.5.",all = FALSE)

  set_ep_warnings <- capture_warnings(fix_div_matrix <- gloscope_proportion(sample_ids,celltype, ep = 0.5, dist_mat= "KL"))
  expect_match(set_ep_warnings,"There are elements haing 0 proportion! ep has been set to be 0.5.",all = FALSE)
  
  expect_equal(isSymmetric(fix_div_matrix),TRUE)
  short_sample_ids <- sample_ids[1:100]
  expect_error(gloscope_proportion(short_sample_ids, celltype),
               regexp="Lengths of cell id and cell type are not equal!",fixed=TRUE)
})
