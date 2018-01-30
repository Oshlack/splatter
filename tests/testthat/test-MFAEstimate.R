context("MFAEstimate")

library(scater)
data("sc_example_counts")

test_that("MFAEstimate works", {
    params <- mfaEstimate(sc_example_counts)
    expect_true(validObject(params))
})
