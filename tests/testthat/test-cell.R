library(testthat)
test_that("gloscope warnings for small cells", {
	expect_error(gloscope(full_pca_embeddings_subset, c(sample_ids,NA)),
		regexp="The number of cells in the embedding matrix does not match the number of sample labels.",fixed=TRUE)
	sample_size_warnings <- capture_warnings(gloscope(reliability_pca_embeddings_subset,reliability_sample_ids,dens="KNN"))
	expect_match(sample_size_warnings,"^The following samples have fewer than 500 cells.*",all = FALSE)
	expect_match(sample_size_warnings,"^The following samples have fewer than the minimum of 50 cells*",all = FALSE)
	#expect_warning(gloscope(reliability_pca_embeddings_subset,sample_ids,dens="KNN"),
	#	regexp = "^The following samples have fewer than 50 cells.*")
  })

test_that("gloscope works with KNN",{
  expect_silent(temp_knn<-gloscope(sub_data,sample_id = "patient_id",
      dim_redu = "PC", ndim = 10, dens = "KNN",
      dist_mat = "KL",BPPARAM = BiocParallel::SerialParam()))
  #test dimensions
  expect_equal(dim(temp_knn),c(3,3))
  #test distances the same
  expect_equal(round(temp_knn[upper.tri(temp_knn)],6),c(2.479180, 2.330122, 3.203009))
  #test diag zero
  expect_equal(unname(diag(temp_knn)),rep(0,3))
  #test row names/colnames
  expect_equal(sort(colnames(temp_knn)),sort(unique(sub_data$patient_id)))
  expect_equal(sort(rownames(temp_knn)),sort(unique(sub_data$patient_id)))
})

test_that("gloscope works with GMM",{
  expect_silent(temp_gmm<-gloscope(sub_data,sample_id = "patient_id",
                                   dim_redu = "PC", ndim = 10, dens = "GMM",
                                   dist_mat = "KL",BPPARAM = BiocParallel::SerialParam(RNGseed = 1)))
  #test dimensions
  expect_equal(dim(temp_gmm),c(3,3))
  #test distances the same
  expect_equal(round(temp_gmm[upper.tri(temp_gmm)],6),c(10.196157,  9.458072, 10.453753))
  #test diag zero
  expect_equal(unname(diag(temp_gmm)),rep(0,3))
  #test row names/colnames
  expect_equal(sort(colnames(temp_gmm)),sort(unique(sub_data$patient_id)))
  expect_equal(sort(rownames(temp_gmm)),sort(unique(sub_data$patient_id)))
})
