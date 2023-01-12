test_that("setBPParam", {
	# check silent run with default arguments	
	expect_silent(
		setBPParam(BPPARAM=NULL,requested_cores=1)	
	)

	# only run this test is Bioc detects multicore available
	local_bpparam <- BiocParallel::bpparam()
	if (class(local_bpparam) == "MulticoreParam"){
		# check warning if more cores requested than available
		available_cores <- parallel::detectCores()
		requested_cores <- available_cores + 1
		expect_warning(
			setBPParam(BPPARAM=NULL,requested_cores=requested_cores),
			regexp="*.reuqested cores not available.*"
		)
	}
})
