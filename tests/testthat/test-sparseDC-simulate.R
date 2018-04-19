context("SparseDC simulation")

test_that("SparseDC simulation output is valid", {
    expect_true(validObject(sparseDCSimulate()))
})
