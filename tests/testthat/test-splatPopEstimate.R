context("splatPopEstimate")

if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
    requireNamespace("preprocessCore", quietly = TRUE)) {

    set.seed(1)

    # Mock data
    bulk.means <- mockBulkMatrix(n.genes=500, n.samples=100)
    bulk.eqtl <- mockBulkeQTL(n.genes=500)

    counts <- scater::mockSCE()
}

test_that("splatPopEstimate works", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    # Check separate functions
    params <- splatPopEstimateMeanCV(newSplatPopParams(), bulk.means)
    expect_true(validObject(params))
    params <- splatPopEstimateEffectSize(newSplatPopParams(), bulk.eqtl)
    expect_true(validObject(params))

    # Check full function
    params <- splatPopEstimate(means=bulk.means,
                               eqtl=bulk.eqtl,
                               counts=counts)
    expect_true(validObject(params))
})

test_that("splatPopEstimate checks on input data", {
    skip_if_not_installed("VariantAnnotation")
    skip_if_not_installed("preprocessCore")
    bulk.eqtl.bad <- bulk.eqtl
    bulk.eqtl.bad$gene_id <- NULL
    expect_error(splatPopEstimate(means=bulk.means,
                                  eqtl=bulk.eqtl.bad,
                                  counts=counts),
        "Incorrect format for eqtl data.")

    bulk.means.bad <- bulk.means
    bulk.means.bad[, 1] <- NA
    expect_error(splatPopEstimate(means=bulk.means.bad,
                                  eqtl=bulk.eqtl,
                                  counts=counts),
        "Incorrect format or NAs present in gene.means.")
})
