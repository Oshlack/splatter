context("BASiCSParams")

test_that("gene.params checks work", {
    params <- newBASiCSParams()
    expect_error(setParam(params, "gene.params", data.frame(A = 1, B = 1)),
                 "gene.params: Incorrect column names")
    expect_error(setParam(params, "gene.params",
                          data.frame(Mean = 1, Disp = "a")),
                 "gene.params: May only contain the following types: numeric")
})

test_that("cell.params checks work", {
    params <- newBASiCSParams()
    expect_error(setParam(params, "cell.params", data.frame(A = 1, B = 1)),
                 "cell.params: Incorrect column names")
    expect_error(setParam(params, "cell.params",
                          data.frame(Phi = 1, S = "a")),
                 "cell.params: May only contain the following types: numeric")
})
