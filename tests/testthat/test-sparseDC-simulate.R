context("SparseDC simulation")

test_that("SparseDC simulation output is valid", {
    skip_if_not_installed("SparseDC")
    expect_true(validObject(sparseDCSimulate()))
})
