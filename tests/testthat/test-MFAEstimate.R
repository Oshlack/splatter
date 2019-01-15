context("MFAEstimate")

library(scater)
data("sc_example_counts")

test_that("MFAEstimate works", {
    skip_if_not_installed("mfa")
    params <- mfaEstimate(sc_example_counts)
    expect_true(validObject(params))
})
