context("Manifold density divergence")

test_that("calc_den", {
  expect_warning(distMat(sub_data_350, sample_id = "sample_id",dim_redu = "PC", ndim = 10, dens = "KNN", distmat = "KL"), "Some samples have numbers of cells smaller than the minimum cell number 500!")
  expect_error(distMat(sub_data_40, sample_id = "sample_id",dim_redu = "PC", ndim = 10, dens = "KNN", distmat = "KL"), "Some samples have numbers of cells smaller than the valid cell number 50 for KNN downstream analysis!")
  expect_silent(distMat(sub_data, sample_id = "sample_id",dim_redu = "PC", ndim = 10, dens = "KNN", distmat = "KL"))
  
  
})

