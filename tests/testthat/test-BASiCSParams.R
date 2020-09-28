context("BASiCSParams")

if (requireNamespace("BASiCS", quietly = TRUE)) {
    params <- newBASiCSParams()
}

test_that("printing works", {
    skip_if_not_installed("BASiCS")
    expect_output(show(params), "A Params object of class BASiCSParams")
})

test_that("gene.params checks work", {
    skip_if_not_installed("BASiCS")
    expect_error(setParam(params, "gene.params", data.frame(A = 1, B = 1)),
                 "gene.params: Incorrect column names")
    expect_error(setParam(params, "gene.params",
                          data.frame(Mean = 1, Disp = "a")),
                 "gene.params: May only contain the following types: \\{numeric\\}")
})

test_that("cell.params checks work", {
    skip_if_not_installed("BASiCS")
    expect_error(setParam(params, "cell.params", data.frame(A = 1, B = 1)),
                 "cell.params: Incorrect column names")
    expect_error(setParam(params, "cell.params",
                          data.frame(Phi = 1, S = "a")),
                 "cell.params: May only contain the following types: \\{numeric\\}")
})

test_that("nBatches checks work", {
    skip_if_not_installed("BASiCS")
    expect_error(setParam(params, "nCells", 1),
                 "nCells cannot be set directly, set batchCells instead")
    expect_error(setParam(params, "nBatches", 1),
                 "nBatches cannot be set directly, set batchCells instead")
})

test_that("batchCells checks work", {
    skip_if_not_installed("BASiCS")
    pp <- setParam(params, "batchCells", c(10, 10))
    expect_equal(getParam(pp, "nCells"), 20)
    expect_equal(getParam(pp, "nBatches"), 2)
})
