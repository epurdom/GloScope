#context("cellNumber")

test_that("distMat", {
  expect_warning(distMat(sub_data_350, sample_id = "patient_id",dim_redu = "PC", ndim = 10, dens = "KNN", dist_mat = "KL", BPPARAM = BiocParallel::SerialParam()), "Some samples have numbers of cells smaller than the minimum cell number 500 to have reliable results!")
  expect_error(distMat(sub_data_40, sample_id = "patient_id",dim_redu = "PC", ndim = 10, dens = "KNN", dist_mat = "KL",BPPARAM = BiocParallel::SerialParam()), "Some samples have numbers of cells smaller than the valid cell number 50 for KNN downstream analysis!")


})

