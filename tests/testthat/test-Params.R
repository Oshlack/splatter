context("Params")

test_that("setParam sets correctly", {
    params <- newSimpleParams()
    params <- setParam(params, "nGenes", 100)
    expect_equal(getParam(params, "nGenes"), 100)
    params <- setParam(params, "mean.rate", 0.5)
    expect_equal(getParam(params, "mean.rate"), 0.5)
})

test_that("getParam gets correctly", {
    params <- newSimpleParams()
    expect_equal(getParam(params, "nGenes"), 10000)
    expect_equal(getParam(params, "mean.rate"), 0.3)
})
