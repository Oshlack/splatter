context("Kersplat simulations")

test.params <- newKersplatParams()

test_that("kersplatSimulate output is valid", {
    expect_true(validObject(kersplatSimulate(test.params)))
})
