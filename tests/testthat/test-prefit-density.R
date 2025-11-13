library(testthat)
library(mclust)

test_that("gloscope validates prefit_density object types", {
    # Use the small test data from setup_data.R
    sample_ids <- subsample_metadata$sample_id
    
    # Create sample embedding list for testing (needed for creating mclust objects)
    embeddings_list <- lapply(unique(sample_ids), function(x) {
        subsample_data_subset[(sample_ids == x), ]
    })
    names(embeddings_list) <- unique(sample_ids)
    
    # Test 1: Error when prefit_density contains Mclust object (wrong type)
    mclust_obj <- mclust::Mclust(
        embeddings_list[[1]], 
        G = 5, 
        modelNames = "VVE",
        verbose = FALSE
    )
    prefit_wrong_mclust <- list(mclust_obj)
    names(prefit_wrong_mclust) <- names(embeddings_list)[1]
    
    expect_error(
        gloscope(
            subsample_data_subset,
            sample_ids,
            dens = "GMM",
            prefit_density = prefit_wrong_mclust
        ),
        regexp = "is a 'Mclust' object, but 'densityMclust' is required"
    )
    
    # Test 2: Error message suggests correct function (densityMclust)
    expect_error(
        gloscope(
            subsample_data_subset,
            sample_ids,
            dens = "GMM",
            prefit_density = prefit_wrong_mclust
        ),
        regexp = "mclust::densityMclust\\(\\)"
    )
    
    # Test 3: Error when prefit_density contains completely wrong object type
    prefit_wrong_type <- list(matrix(1:10, nrow = 5))
    names(prefit_wrong_type) <- names(embeddings_list)[1]
    
    expect_error(
        gloscope(
            subsample_data_subset,
            sample_ids,
            dens = "GMM",
            prefit_density = prefit_wrong_type
        ),
        regexp = "must be a 'densityMclust' object"
    )
    
    # Test 4: Error message includes object name for named list
    expect_error(
        gloscope(
            subsample_data_subset,
            sample_ids,
            dens = "GMM",
            prefit_density = prefit_wrong_mclust
        ),
        regexp = names(embeddings_list)[1]
    )
    
    # Test 5: Error message includes index [[1]] for unnamed list
    prefit_unnamed <- list(mclust_obj)
    expect_error(
        gloscope(
            subsample_data_subset,
            sample_ids,
            dens = "GMM",
            prefit_density = prefit_unnamed
        ),
        regexp = "\\[\\[1\\]\\]"
    )
    
    # Test 6: Success with correct densityMclust objects
    # Create correct densityMclust objects for all samples
    prefit_correct <- lapply(names(embeddings_list), function(sample_name) {
        mclust::densityMclust(
            embeddings_list[[sample_name]],
            G = 5,
            modelNames = "VVE",
            plot = FALSE,
            verbose = FALSE
        )
    })
    names(prefit_correct) <- names(embeddings_list)
    
    expect_silent(
        result <- gloscope(
            subsample_data_subset,
            sample_ids,
            dens = "GMM",
            prefit_density = prefit_correct,
            BPPARAM = BiocParallel::SerialParam(RNGseed = 2)
        )
    )
    
    # Verify result is a valid divergence matrix
    expect_equal(dim(result), c(length(unique(sample_ids)), length(unique(sample_ids))))
    expect_equal(isSymmetric(result), TRUE)
    expect_equal(unname(diag(result)), rep(0, length(unique(sample_ids))))
})



