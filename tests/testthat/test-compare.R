context("Comparison functions")

sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)

comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")

test_that("compareSCEs works", {
    expect_length(comparison, 3)
    expect_true(all(c("RowData", "ColData", "Plots") %in%
                        names(comparison)))
    checkmate::expect_class(comparison$ColData, "data.frame")
    checkmate::expect_class(comparison$RowData, "data.frame")
    expect_length(comparison$Plots, 8)
    expect_true(all(c("Means", "Variances", "MeanVar", "LibrarySizes",
                      "ZerosGene", "ZerosCell", "MeanZeros", "VarGeneCor") %in%
                        names(comparison$Plots)))
    for (plot in names(comparison$Plots)) {
        checkmate::expect_class(comparison$Plots[[plot]], "ggplot")
    }
})

test_that("diffSCEs works", {
    expect_length(difference, 5)
    expect_true(all(c("Reference", "RowData", "ColData", "Plots",
                      "QQPlots") %in% names(difference)))
    checkmate::expect_class(difference$Reference, "SingleCellExperiment")
    checkmate::expect_class(difference$ColData, "data.frame")
    checkmate::expect_class(difference$RowData, "data.frame")
    expect_length(difference$Plots, 7)
    expect_true(all(c("Means", "Variances", "MeanVar", "LibrarySizes",
                      "ZerosGene", "ZerosCell", "MeanZeros") %in%
                        names(difference$Plots)))
    for (plot in names(difference$Plots)) {
        checkmate::expect_class(difference$Plots[[plot]], "ggplot")
    }
    expect_length(difference$QQPlots, 5)
    expect_true(all(c("Means", "Variances", "LibrarySizes", "ZerosGene",
                      "ZerosCell") %in% names(difference$QQPlots)))
    for (plot in names(difference$QQPlots)) {
        checkmate::expect_class(difference$QQPlots[[plot]], "ggplot")
    }
})

# test_that("makeCompPanel works", {
#     panel <- makeCompPanel(comparison)
#     checkmate::expect_class(panel, "ggplot")
# })
#
# test_that("makeDiffPanel works", {
#     panel <- makeDiffPanel(difference)
#     checkmate::expect_class(panel, "ggplot")
# })
#
# test_that("makeOverallPanel works", {
#     panel <- makeOverallPanel(comparison, difference)
#     checkmate::expect_class(panel, "ggplot")
# })
