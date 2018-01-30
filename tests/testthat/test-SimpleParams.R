context("SimpleParams")

params <- newSimpleParams()

test_that("printing works", {
    expect_output(show(params), "A Params object of class SimpleParams")
})
