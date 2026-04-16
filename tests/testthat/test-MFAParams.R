context("MFAParams")

test_that("printing works", {
    skip_if_not_installed("mfa")

    # Force deprecation warnings to be issued
    withr::local_options(lifecycle_verbosity = "warning")

    params <- expect_warning(newMFAParams(), "deprecated")
    expect_output(show(params), "MFAParams")
})
