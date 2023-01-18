test_that("calc_den", {
  expect_silent(
    mod_list <- calc_dens(df_list, dens = "GMM", k = 50, BPPARAM = BiocParallel::SerialParam()))
  expect_silent(mod_list <- calc_dens(df_list, dens = "KNN", k = 50, BPPARAM = BiocParallel::SerialParam()))

  expect_true(length(mod_list) == length(df_list))

})

