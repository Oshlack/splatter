context("splatPopEstimate")
set.seed(42)

# Mock data
bulk_means <- mock_bulk_matrix(n_genes=500, n_samples=100)
bulk_eqtl <- mock_bulk_eqtl(n_genes=500)
counts <- mockSCE()

test_that("splatPopEstimate works", {
    params <- splatPopEstimate(means=bulk_means, 
                               eqtl=bulk_eqtl,
                               counts=counts)
    expect_true(validObject(params))
})

test_that("splatPopEstimate checks on input data", {
    bulk_eqtl_bad <- bulk_eqtl
    bulk_eqtl_bad$gene_id <- NULL
    expect_error(splatPopEstimate(means=bulk_means, 
                                  eqtl=bulk_eqtl_bad, 
                                  counts=counts),
        "Incorrect format for eqtl data.")
    
    bulk_means_bad <- bulk_means
    bulk_means_bad$new_sample <- NA
    expect_error(splatPopEstimate(means=bulk_means_bad, 
                                  eqtl=bulk_eqtl, 
                                  counts=counts),
        "Incorrect format or NAs present in gene.means.")
})

