context("Lun2Params")

params <- newLun2Params()

test_that("printing works", {
    expect_output(show(params), "A Params object of class Lun2Params")
})

test_that("nCells checks work", {
    expect_error(setParam(params, "nCells", 1),
                 "nCells cannot be set directly, set cell.plates instead")
    expect_error(setParam(params, "nPlates", 1),
                 "nPlates cannot be set directly, set cell.plates instead")
})

test_that("gene.params checks work", {
    expect_error(setParam(params, "gene.params", data.frame(A = 1, B = 1)),
                 "gene.params: Incorrect column names")
    expect_error(setParam(params, "gene.params",
                          data.frame(Mean = 1, Disp = "a")),
                 "gene.params: May only contain the following types: \\{numeric\\}")
})
