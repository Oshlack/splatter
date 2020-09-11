context("popEstimate")
set.seed(1)

# Mock bulk expression matrix
bulk_means <- data.frame(mean = rnbinom(1e3, 1,0.1))
bulk_matrix <- data.frame(abs(t(apply(bulk_means, 1, 
                                      function(x) rnorm(50, x, 1)))))

# Mock eQTL mapping results
bulk_eqtl <- data.frame(list(gene_id=1:1e3, pval_nominal=0.05,
                             slope=rnorm(1e3, 0, 0.5)))

test_that("popEstimate works", {
    params <- popEstimate(gene.means=bulk_matrix, eqtl=bulk_eqtl)
    expect_true(validObject(params))
})

test_that("popEstimate checks on input data", {
    bulk_eqtl_bad <- bulk_eqtl
    bulk_eqtl_bad$gene_id <- NULL
    expect_error(popEstimate(gene.means=bulk_matrix, eqtl=bulk_eqtl_bad),
        "Incorrect format for eqtl data.")
    
    bulk_matrix_bad <- bulk_matrix
    bulk_matrix_bad$new_sample <- NA
    expect_error(popEstimate(gene.means=bulk_matrix_bad, eqtl=bulk_eqtl),
        "Incorrect format or NAs present in gene.means.")
})

