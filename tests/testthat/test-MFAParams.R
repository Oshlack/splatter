context("MFAParams")

params <- newMFAParams()

test_that("printing works", {
    expect_output(show(params), "A Params object of class MFAParams")
})
