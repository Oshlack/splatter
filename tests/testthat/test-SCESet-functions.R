context("SCESet functions")

test_that("addGeneLengths generate method works", {
    sce <- simpleSimulate()
    expect_silent(addGeneLengths(sce))
    expect_error(addGeneLengths("a"))
    expect_error(addGeneLengths(sce, loc = "a"))
    expect_error(addGeneLengths(sce, scale = "a"))
    expect_error(addGeneLengths(sce, scale = -1))
})

test_that("addGeneLength sample method works", {
    sce <- simpleSimulate()
    lens <- round(runif(100, 100, 10000))
    expect_silent(addGeneLengths(sce, method = "sample", lengths = lens))
    expect_error(addGeneLengths(sce, method = "sample"))
    expect_error(addGeneLengths(sce, method = "sample"), lengths = 0)
    expect_error(addGeneLengths(sce, method = "sample"), lengths = "a")
})
