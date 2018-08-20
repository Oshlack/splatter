context("BASiCSEstimate")

library(scater)
data("sc_example_counts")

test_that("BASiCSEstimate works", {
    skip_if_not_installed("BASiCS")
    set.seed(1)
    spike.info <- data.frame(Name = rownames(sc_example_counts)[1:10],
                             Input = rnorm(10, 500, 200),
                             stringsAsFactors = FALSE)
    counts <- sc_example_counts[, 1:20]
    counts <- counts[rowSums(counts) != 0, ]
    params <- BASiCSEstimate(counts[1:100, ],
                             spike.info, verbose = FALSE, progress = FALSE)
    expect_true(validObject(params))
})

test_that("BASiCSEstimate works without spikes", {
    skip_if_not_installed("BASiCS")
    set.seed(1)
    counts <- sc_example_counts[, 1:20]
    counts <- counts[rowSums(counts) != 0, ]
    batch <- sample(1:2, ncol(counts), replace = TRUE)
    params <- BASiCSEstimate(counts[1:100, ], batch = batch,
                             verbose = FALSE, progress = FALSE)
    expect_true(validObject(params))
})
