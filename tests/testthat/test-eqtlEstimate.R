context("eQTLEstimate")

data(ex_means)
data(ex_pairs)

test_that("eQTLEstimate works", {
    params <- eQTLEstimate()
    expect_true(validObject(params))
})

test_that("eQTLEstimate checks on input data", {
    ex_pairs_bad <- ex_pairs
    ex_pairs_bad$gene_id <- NULL
    expect_error(eQTLEstimate(all.pairs = ex_pairs_bad),
        "Incorrect format for all.pairs. See example data.")
    
    ex_means_bad <- ex_means
    ex_means_bad$new_sample <- NA
    expect_error(eQTLEstimate(gene.means = ex_means_bad),
        "Incorrect format or NAs present in gene.means. See example data.")
})

