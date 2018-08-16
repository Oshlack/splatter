context("MFAParams")

test_that("printing works", {
    skip_if_not_installed("mfa")
    params <- newMFAParams()
    expect_output(show(params), "A Params object of class MFAParams")
})
