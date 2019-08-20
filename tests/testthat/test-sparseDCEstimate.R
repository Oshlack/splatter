context("sparseDCEstimate")

library(scater)
set.seed(1)
counts <- counts(mockSCE())

test_that("sparseDCEstimate works", {
    skip_if_not_installed("SparseDC")
    set.seed(1)
    conditions <- sample(1:2, ncol(counts), replace = TRUE)
    params <- sparseDCEstimate(counts, conditions, nclusters = 3)
    expect_true(validObject(params))
})
