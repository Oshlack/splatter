context("SparseDCParams")

if (requireNamespace("SparseDC", quietly = TRUE)) {
    params <- newSparseDCParams()
}

test_that("constructor is valid", {
    skip_if_not_installed("SparseDC")
    expect_true(validObject(params))
})

test_that("printing works", {
    skip_if_not_installed("SparseDC")
    expect_output(show(params), "A Params object of class SparseDCParams")
})

test_that("clusts checks work", {
    skip_if_not_installed("SparseDC")
    expect_error(setParam(params, "clusts.c1", 2),
                 "is length 1 it must equal 1")
    expect_error(setParam(params, "clusts.c2", 2),
                 "is length 1 it must equal 1")
    expect_error(setParam(params, "clusts.c1", c(1, 3)),
                 "Cluster labels must be sequential")
    expect_error(setParam(params, "clusts.c2", c(1, 3)),
                 "Cluster labels must be sequential")
})
