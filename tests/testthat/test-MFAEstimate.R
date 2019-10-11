context("MFAEstimate")

if (requireNamespace("mfa", quietly = TRUE)) {
    library(mfa)
    synth <- create_synthetic(C = 20, G = 5, zero_negative = TRUE,
                              model_dropout = TRUE)
}

test_that("MFAEstimate works", {
    skip_if_not_installed("mfa")
    params <- mfaEstimate(synth$X)
    expect_true(validObject(params))
})
