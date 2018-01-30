context("lunEstimate")

library(scater)
data("sc_example_counts")

test_that("lunEstimate works", {
    params <- lunEstimate(sc_example_counts)
    expect_true(validObject(params))
})
