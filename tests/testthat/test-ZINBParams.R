context("ZINBParams")

if (requireNamespace("zinbwave", quietly = TRUE)) {
    params <- newZINBParams()
}

test_that("constructor is valid", {
    skip_if_not_installed("zinbwave")
    expect_true(validObject(params))
})

test_that("printing works", {
    skip_if_not_installed("zinbwave")
    expect_output(show(params), "A Params object of class ZINBParams")
})

test_that("nGenes checks work", {
    skip_if_not_installed("zinbwave")
    expect_error(setParam(params, "nGenes", 1),
                 "nGenes cannot be set directly")
    expect_error(setParam(params, "nCells", 1),
                 "nCells cannot be set directly")
})
