context("simpleEstimate")

library(scater)
set.seed(1)
counts <- counts(mockSCE())

test_that("simpleEstimate works", {
    params <- simpleEstimate(counts)
    expect_true(validObject(params))
})
