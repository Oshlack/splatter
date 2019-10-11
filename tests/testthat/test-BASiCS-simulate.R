context("BASiCS simulation")

test_that("BASiCS simulation output is valid", {
    skip_if_not_installed("BASiCS")
    expect_true(validObject(BASiCSSimulate()))
})
