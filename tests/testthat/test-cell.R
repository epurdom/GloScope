test_that("gloscope", {
	testthat::expect_error(gloscope(full_pca_embeddings_subset, c(sample_ids,NA)), 
		regexp="The number of cells in the embedding matrix does not match the number of sample labels.",fixed=TRUE)
	sample_size_warnings <- testthat::capture_warnings(gloscope(reliability_pca_embeddings_subset,reliability_sample_ids,dens="KNN"))
	testthat::expect_match(sample_size_warnings,"^The following samples have fewer than 500 cells.*",all = FALSE)
	testthat::expect_match(sample_size_warnings,"^The following samples have fewer than the minimum of 50 cells*",all = FALSE)
	#expect_warning(gloscope(reliability_pca_embeddings_subset,sample_ids,dens="KNN"),
	#	regexp = "^The following samples have fewer than 50 cells.*")
})

