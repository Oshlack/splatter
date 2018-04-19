context("sparseDCEstimate")

library(scater)
data("sc_example_counts")

test_that("sparseDCEstimate works", {
    conditions <- sample(1:2, ncol(sc_example_counts), replace = TRUE)

    params <- sparseDCEstimate(sc_example_counts, conditions,
                               nclusters = 3)
    expect_true(validObject(params))
})
