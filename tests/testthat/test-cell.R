test_that("gloscope", {
	expect_error(gloscope(pca_embeddings_subset, c(sample_ids,NA)), 
		regexp="The number of cells in the embedding matrix does not match the number of sample labels.",fixed=TRUE),
	expect_warning(gloscope(pca_embeddings_subset,sample_ids,dens="KNN"),
		regexp = "^The following samples have fewer than 50 cells.*")
})

