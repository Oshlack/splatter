context("simpleEstimate")

library(scater)
data("sc_example_counts")

test_that("simpleEstimate works", {
    params <- simpleEstimate(sc_example_counts)
    expect_true(validObject(params))
})
