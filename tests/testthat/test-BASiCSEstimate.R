context("BASiCSEstimate")

library(scater)
set.seed(1)
counts <- counts(mockSCE(ncells = 30, ngenes = 100))

test_that("BASiCSEstimate works", {
    skip_if_not_installed("BASiCS")
    set.seed(1)
    spike.info <- data.frame(Name = rownames(counts)[1:10],
                             Input = rnorm(10, 500, 200),
                             stringsAsFactors = FALSE)
    counts <- counts[rowSums(counts) != 0, ]
    params <- BASiCSEstimate(counts, spike.info, verbose = FALSE,
                             progress = FALSE)
    expect_true(validObject(params))
})

test_that("BASiCSEstimate works without spikes", {
    skip_if_not_installed("BASiCS")
    set.seed(1)
    counts <- counts[rowSums(counts) != 0, ]
    batch <- sample(1:2, ncol(counts), replace = TRUE)
    params <- BASiCSEstimate(counts, batch = batch,
                             verbose = FALSE, progress = FALSE)
    expect_true(validObject(params))
})
