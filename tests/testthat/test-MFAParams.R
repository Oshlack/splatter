context("MFAParams")

test_that("printing works", {
    skip_if_not_installed("mfa")
    params <- expect_warning(newMFAParams(), "deprecated")
    expect_output(show(params), "MFAParams")
})
