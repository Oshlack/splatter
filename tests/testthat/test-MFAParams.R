context("MFAParams")

test_that("printing works", {
    skip_if_not_installed("mfa")
    params <- newMFAParams()
    expect_output(show(params), "MFAParams")
})
