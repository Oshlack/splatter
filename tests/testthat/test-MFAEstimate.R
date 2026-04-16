context("MFAEstimate")

if (requireNamespace("mfa", quietly = TRUE)) {
    library(mfa)
    synth <- create_synthetic(
        C = 20, G = 5, zero_negative = TRUE,
        model_dropout = TRUE
    )
}

test_that("MFAEstimate works", {
    skip_if_not_installed("mfa")

    # Force deprecation warnings to be issued
    withr::local_options(lifecycle_verbosity = "warning")

    params <- expect_warning(mfaEstimate(synth$X), "deprecated")
    expect_true(validObject(params))
})
