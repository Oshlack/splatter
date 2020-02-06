context("Kersplat simulations")

test_that("kersplatSimulate output is valid", {
    skip_if_not_installed("igraph")
    test.params <- newKersplatParams()
    expect_true(validObject(kersplatSimulate(test.params)))
})
