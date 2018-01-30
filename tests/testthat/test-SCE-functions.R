context("SCE functions")

sce <- simpleSimulate()

test_that("addFeatureStats works with counts", {
    ss <- addFeatureStats(sce)
    expect_true(all(c("MeanCounts", "VarCounts", "CVCounts", "MedCounts",
                      "MADCounts") %in% colnames(rowData(ss))))
    ss <- addFeatureStats(sce, log = TRUE)
    expect_true(all(c("MeanLogCounts", "VarLogCounts", "CVLogCounts",
                      "MedLogCounts", "MADLogCounts") %in%
                    colnames(rowData(ss))))
    ss <- addFeatureStats(sce, no.zeros = TRUE)
    expect_true(all(c("MeanCountsNo0", "VarCountsNo0", "CVCountsNo0",
                      "MedCountsNo0", "MADCountsNo0") %in%
                        colnames(rowData(ss))))
    ss <- addFeatureStats(sce, log = TRUE, no.zeros = TRUE)
    expect_true(all(c("MeanLogCountsNo0", "VarLogCountsNo0", "CVLogCountsNo0",
                      "MedLogCountsNo0", "MADLogCountsNo0") %in%
                        colnames(rowData(ss))))
})

test_that("addGeneLengths generate method works", {
    expect_silent(addGeneLengths(sce))
    expect_error(addGeneLengths("a"))
    expect_error(addGeneLengths(sce, loc = "a"))
    expect_error(addGeneLengths(sce, scale = "a"))
    expect_error(addGeneLengths(sce, scale = -1))
})

test_that("addGeneLength sample method works", {
    lens <- round(runif(100, 100, 10000))
    expect_silent(addGeneLengths(sce, method = "sample", lengths = lens))
    expect_error(addGeneLengths(sce, method = "sample"))
    expect_error(addGeneLengths(sce, method = "sample"), lengths = 0)
    expect_error(addGeneLengths(sce, method = "sample"), lengths = "a")
})
