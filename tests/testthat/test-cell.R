#context("Gloscope")

test_that("gloscope warnings for small cells", {
  expect_warning(gloscope(sub_data_350, sample_id = "patient_id",dim_redu = "PC", ndim = 10, dens = "KNN", dist_mat = "KL", BPPARAM = BiocParallel::SerialParam()), "Some samples have numbers of cells smaller than the minimum cell number 500 to have reliable results!")
  expect_error(gloscope(sub_data_40, sample_id = "patient_id",dim_redu = "PC", ndim = 10, dens = "KNN", dist_mat = "KL",BPPARAM = BiocParallel::SerialParam()), "Some samples have numbers of cells smaller than the valid cell number 50 for KNN downstream analysis!")


  })

test_that("gloscope works with KNN",{
  expect_silent(temp_knn<-gloscope(sub_data,sample_id = "patient_id",
      dim_redu = "PC", ndim = 10, dens = "KNN",
      dist_mat = "KL",BPPARAM = BiocParallel::SerialParam()))
  #test dimensions
  expect_equal(dim(temp_knn),c(3,3))
  #test distances the same as got the first time ran this (check changes in the code haven't changed anything)
  knn_expected_answer<-c(4.478759 ,3.203009, 1.520805) #answer got the first time
  expect_equal(round(temp_knn[upper.tri(temp_knn)],6),knn_expected_answer)
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
  gmm_expected_answer<-round(c(12.996519, 10.759409, 10.517520),6) #answer got the first time
  expect_equal(round(temp_gmm[upper.tri(temp_gmm)],6),gmm_expected_answer)
  #test diag zero
  expect_equal(unname(diag(temp_gmm)),rep(0,3))
  #test row names/colnames
  expect_equal(sort(colnames(temp_gmm)),sort(unique(sub_data$patient_id)))
  expect_equal(sort(rownames(temp_gmm)),sort(unique(sub_data$patient_id)))
})
