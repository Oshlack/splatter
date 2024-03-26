#' Compare SingleCellExperiment objects
#'
#' Combine the data from several SingleCellExperiment objects and produce some
#' basic plots comparing them.
#'
#' @param sces named list of SingleCellExperiment objects to combine and
#'        compare.
#' @param point.size size of points in scatter plots.
#' @param point.alpha opacity of points in scatter plots.
#' @param fits whether to include fits in scatter plots.
#' @param colours vector of colours to use for each dataset.
#'
#' @details
#' The returned list has three items:
#'
#' \describe{
#'     \item{\code{RowData}}{Combined row data from the provided
#'     SingleCellExperiments.}
#'     \item{\code{ColData}}{Combined column data from the provided
#'     SingleCellExperiments.}
#'     \item{\code{Plots}}{Comparison plots
#'         \describe{
#'             \item{\code{Means}}{Boxplot of mean distribution.}
#'             \item{\code{Variances}}{Boxplot of variance distribution.}
#'             \item{\code{MeanVar}}{Scatter plot with fitted lines showing the
#'             mean-variance relationship.}
#'             \item{\code{LibrarySizes}}{Boxplot of the library size
#'             distribution.}
#'             \item{\code{ZerosGene}}{Boxplot of the percentage of each gene
#'             that is zero.}
#'             \item{\code{ZerosCell}}{Boxplot of the percentage of each cell
#'             that is zero.}
#'             \item{\code{MeanZeros}}{Scatter plot with fitted lines showing
#'             the mean-zeros relationship.}
#'             \item{\code{VarGeneCor}}{Heatmap of correlation of the 100 most
#'             variable genes.}
#'     }
#'   }
#' }
#'
#' The plots returned by this function are created using
#' \code{\link[ggplot2]{ggplot}} and are only a sample of the kind of plots you
#' might like to consider. The data used to create these plots is also returned
#' and should be in the correct format to allow you to create further plots
#' using \code{\link[ggplot2]{ggplot}}.
#'
#' @return List containing the combined datasets and plots.
#' @examples
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
#' names(comparison)
#' names(comparison$Plots)
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom SingleCellExperiment cpm<- cpm
#' @importFrom rlang .data
#' @importFrom stats cor
#' @export
compareSCEs <- function(sces, point.size = 0.1, point.alpha = 0.1,
                        fits = TRUE, colours = NULL) {
    checkDependencies(deps = "ggplot2")

    checkmate::assertList(
        sces,
        types = "SingleCellExperiment",
        any.missing = FALSE,
        min.len = 1,
        names = "unique"
    )
    checkmate::assertNumber(point.size, finite = TRUE)
    checkmate::assertNumber(point.alpha, lower = 0, upper = 1)
    checkmate::assertLogical(fits, any.missing = FALSE, len = 1)

    if (!is.null(colours)) {
        checkmate::assertCharacter(
            colours,
            any.missing = FALSE,
            len = length(sces)
        )
    } else {
        colours <- scales::hue_pal()(length(sces))
    }

    for (name in names(sces)) {
        sce <- sces[[name]]
        rowData(sce)$Dataset <- name
        colData(sce)$Dataset <- name
        sce <- scuttle::addPerCellQC(sce)
        sce <- scuttle::addPerFeatureQC(sce)
        cpm(sce) <- as.matrix(scuttle::calculateCPM(sce))
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        n.features <- colData(sce)$detected
        colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
        rowData(sce)$PctZero <- 100 - rowData(sce)$detected
        var.genes <- rev(order(rowData(sce)$VarLogCPM))[seq_len(100)]
        var.cpm <- log2(cpm(sce)[var.genes, ] + 1)
        var.cors <- as.data.frame.table(cor(t(var.cpm), method = "spearman"))
        colnames(var.cors) <- c("GeneA", "GeneB", "Correlation")
        var.cors$VarGeneA <- rep(paste0("VarGene", seq_len(100)), 100)
        var.cors$VarGeneB <- rep(paste0("VarGene", seq_len(100)), each = 100)
        var.cors$Dataset <- name
        var.cors <- var.cors[, c(
            "Dataset", "GeneA", "GeneB", "VarGeneA",
            "VarGeneB", "Correlation"
        )]
        metadata(sce)$VarGeneCorrelation <- var.cors

        sces[[name]] <- sce
    }

    features <- rowData(sces[[1]])
    cells <- colData(sces[[1]])
    var.cors <- metadata(sces[[1]])$VarGeneCorrelation

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
            var.cors <- rbindMatched(var.cors, metadata(sce)$VarGeneCorrelation)
        }
    }

    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    var.cors$Dataset <- factor(var.cors$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)

    means <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$MeanLogCPM, colour = .data$Dataset
        )
    ) +
        ggplot2::geom_violin(ggplot2::aes(fill = .data$Dataset),
            draw_quantiles = c(0.25, 0.5, 0.75),
            colour = "white", alpha = 0.3, linewidth = 0.8
        ) +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::ylab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ggplot2::ggtitle("Distribution of mean expression") +
        ggplot2::theme_minimal()

    vars <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$VarLogCPM, colour = .data$Dataset
        )
    ) +
        ggplot2::geom_violin(ggplot2::aes(fill = .data$Dataset),
            draw_quantiles = c(0.25, 0.5, 0.75),
            colour = "white", alpha = 0.3, linewidth = 0.8
        ) +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::ylab(expression(paste("Variance ", log[2], "(CPM + 1)"))) +
        ggplot2::ggtitle("Distribution of variance") +
        ggplot2::theme_minimal()

    mean.var <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$MeanLogCPM, y = .data$VarLogCPM,
            colour = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_point(size = point.size, alpha = point.alpha) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ggplot2::ylab(expression(paste("Variance ", log[2], "(CPM + 1)"))) +
        ggplot2::ggtitle("Mean-variance relationship") +
        ggplot2::theme_minimal()

    libs <- ggplot2::ggplot(
        cells,
        ggplot2::aes(x = .data$Dataset, y = .data$sum, colour = .data$Dataset)
    ) +
        ggplot2::geom_violin(
            ggplot2::aes(fill = .data$Dataset),
            draw_quantiles = c(0.25, 0.5, 0.75),
            colour = "white", alpha = 0.3, linewidth = 0.8
        ) +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_y_continuous(labels = scales::comma) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::ylab("Total counts per cell") +
        ggplot2::ggtitle("Distribution of library sizes") +
        ggplot2::theme_minimal()

    z.gene <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$PctZero,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_violin(
            ggplot2::aes(fill = .data$Dataset),
            draw_quantiles = c(0.25, 0.5, 0.75),
            colour = "white", alpha = 0.3, linewidth = 0.8
        ) +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_y_continuous(limits = c(0, 100)) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::ylab("Percentage zeros per gene") +
        ggplot2::ggtitle("Distribution of zeros per gene") +
        ggplot2::theme_minimal()

    z.cell <- ggplot2::ggplot(
        cells,
        ggplot2::aes(
            x = .data$Dataset, y = .data$PctZero,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_violin(
            ggplot2::aes(fill = .data$Dataset),
            draw_quantiles = c(0.25, 0.5, 0.75),
            colour = "white", alpha = 0.3, linewidth = 0.8
        ) +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_y_continuous(limits = c(0, 100)) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::ylab("Percentage zeros per cell") +
        ggplot2::ggtitle("Distribution of zeros per cell") +
        ggplot2::theme_minimal()

    mean.zeros <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$mean, y = .data$PctZero,
            colour = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_point(size = point.size, alpha = point.alpha) +
        ggplot2::scale_x_log10(labels = scales::comma) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("Mean count") +
        ggplot2::ylab("Percentage zeros") +
        ggplot2::ggtitle("Mean-zeros relationship") +
        ggplot2::theme_minimal()

    var.correlation <- ggplot2::ggplot(
        var.cors,
        ggplot2::aes(
            x = .data$VarGeneA, y = .data$VarGeneB,
            fill = .data$Correlation
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) +
        ggplot2::coord_fixed() +
        ggplot2::facet_wrap(~Dataset) +
        ggplot2::ggtitle("Correlation - 100 variable genes") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank()
        )

    if (fits) {
        mean.var <- mean.var +
            ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
        mean.zeros <- mean.zeros +
            ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
    }

    comparison <- list(
        RowData = features,
        ColData = cells,
        Plots = list(
            Means = means,
            Variances = vars,
            MeanVar = mean.var,
            LibrarySizes = libs,
            ZerosGene = z.gene,
            ZerosCell = z.cell,
            MeanZeros = mean.zeros,
            VarGeneCor = var.correlation
        )
    )

    return(comparison)
}

#' Diff SingleCellExperiment objects
#'
#' Combine the data from several SingleCellExperiment objects and produce some
#' basic plots comparing them to a reference.
#'
#' @param sces named list of SingleCellExperiment objects to combine and
#'        compare.
#' @param ref string giving the name of the SingleCellExperiment to use as the
#'        reference
#' @param point.size size of points in scatter plots.
#' @param point.alpha opacity of points in scatter plots.
#' @param fits whether to include fits in scatter plots.
#' @param colours vector of colours to use for each dataset.
#'
#' @details
#'
#' This function aims to look at the differences between a reference
#' SingleCellExperiment and one or more others. It requires each
#' SingleCellExperiment to have the same dimensions. Properties are compared by
#' ranks, for example when comparing the means the values are ordered and the
#' differences between the reference and another dataset plotted. A series of
#' Q-Q plots are also returned.
#'
#' The returned list has five items:
#'
#' \describe{
#'     \item{\code{Reference}}{The SingleCellExperiment used as the reference.}
#'     \item{\code{RowData}}{Combined feature data from the provided
#'     SingleCellExperiments.}
#'     \item{\code{ColData}}{Combined column data from the provided
#'     SingleCellExperiments.}
#'     \item{\code{Plots}}{Difference plots
#'         \describe{
#'             \item{\code{Means}}{Boxplot of mean differences.}
#'             \item{\code{Variances}}{Boxplot of variance differences.}
#'             \item{\code{MeanVar}}{Scatter plot showing the difference from
#'             the reference variance across expression ranks.}
#'             \item{\code{LibraeySizes}}{Boxplot of the library size
#'             differences.}
#'             \item{\code{ZerosGene}}{Boxplot of the differences in the
#'             percentage of each gene that is zero.}
#'             \item{\code{ZerosCell}}{Boxplot of the differences in the
#'             percentage of each cell that is zero.}
#'             \item{\code{MeanZeros}}{Scatter plot showing the difference from
#'             the reference percentage of zeros across expression ranks.}
#'     }
#'   }
#'   \item{\code{QQPlots}}{Quantile-Quantile plots
#'       \describe{
#'           \item{\code{Means}}{Q-Q plot of the means.}
#'           \item{\code{Variances}}{Q-Q plot of the variances.}
#'           \item{\code{LibrarySizes}}{Q-Q plot of the library sizes.}
#'           \item{\code{ZerosGene}}{Q-Q plot of the percentage of zeros per
#'           gene.}
#'           \item{\code{ZerosCell}}{Q-Q plot of the percentage of zeros per
#'           cell.}
#'       }
#'   }
#' }
#'
#' The plots returned by this function are created using
#' \code{\link[ggplot2]{ggplot}} and are only a sample of the kind of plots you
#' might like to consider. The data used to create these plots is also returned
#' and should be in the correct format to allow you to create further plots
#' using \code{\link[ggplot2]{ggplot}}.
#'
#' @return List containing the combined datasets and plots.
#' @examples
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' names(difference)
#' names(difference$Plots)
#' @importFrom SingleCellExperiment cpm<-
#' @importFrom rlang .data
#' @export
diffSCEs <- function(sces, ref, point.size = 0.1, point.alpha = 0.1,
                     fits = TRUE, colours = NULL) {
    checkDependencies(deps = "ggplot2")

    checkmate::assertList(
        sces,
        types = "SingleCellExperiment",
        any.missing = FALSE,
        min.len = 2,
        names = "unique"
    )
    checkmate::assertString(ref)
    checkmate::assertNumber(point.size, finite = TRUE)
    checkmate::assertNumber(point.alpha, lower = 0, upper = 1)
    checkmate::assertLogical(fits, any.missing = FALSE, len = 1)

    if (!(ref %in% names(sces))) {
        stop("'ref' must be the name of a SingleCellExperiment in 'sces'")
    } else {
        ref.idx <- which(names(sces) == ref)
    }

    if (!is.null(colours)) {
        checkmate::assertCharacter(
            colours,
            any.missing = FALSE,
            len = length(sces) - 1
        )
    } else {
        colours <- scales::hue_pal()(length(sces))
        colours <- colours[-ref.idx]
    }

    ref.dim <- dim(sces[[ref]])

    for (name in names(sces)) {
        sce <- sces[[name]]
        if (!identical(dim(sce), ref.dim)) {
            stop("all datasets in 'sces' must have the same dimensions")
        }
        rowData(sce)$Dataset <- name
        colData(sce)$Dataset <- name
        sce <- scuttle::addPerCellQC(sce)
        sce <- scuttle::addPerFeatureQC(sce)
        cpm(sce) <- as.matrix(scuttle::calculateCPM(sce))
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        n.features <- colData(sce)$detected
        colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
        rowData(sce)$RankCounts <- rank(rowData(sce)$mean)
        rowData(sce)$PctZero <- 100 - rowData(sce)$detected
        sces[[name]] <- sce
    }

    ref.sce <- sces[[ref]]

    ref.means <- sort(rowData(ref.sce)$MeanLogCPM)
    ref.vars <- sort(rowData(ref.sce)$VarLogCPM)
    ref.libs <- sort(colData(ref.sce)$sum)
    ref.z.gene <- sort(rowData(ref.sce)$PctZero)
    ref.z.cell <- sort(colData(ref.sce)$PctZero)

    ref.rank.ord <- order(rowData(ref.sce)$RankCounts)
    ref.vars.rank <- rowData(ref.sce)$VarLogCPM[ref.rank.ord]
    ref.z.gene.rank <- rowData(ref.sce)$PctZero[ref.rank.ord]

    for (name in names(sces)) {
        sce <- sces[[name]]
        rowData(sce)$RefRankMeanLogCPM <- ref.means[
            rank(rowData(sce)$MeanLogCPM)
        ]
        rowData(sce)$RankDiffMeanLogCPM <- rowData(sce)$MeanLogCPM -
            rowData(sce)$RefRankMeanLogCPM
        rowData(sce)$RefRankVarLogCPM <- ref.vars[rank(rowData(sce)$VarLogCPM)]
        rowData(sce)$RankDiffVarLogCPM <- rowData(sce)$VarLogCPM -
            rowData(sce)$RefRankVarLogCPM
        colData(sce)$RefRankLibSize <- ref.libs[rank(colData(sce)$sum)]
        colData(sce)$RankDiffLibSize <- colData(sce)$sum -
            colData(sce)$RefRankLibSize
        rowData(sce)$RefRankZeros <- ref.z.gene[rank(rowData(sce)$PctZero)]
        rowData(sce)$RankDiffZeros <- rowData(sce)$PctZero -
            rowData(sce)$RefRankZeros
        colData(sce)$RefRankZeros <- ref.z.cell[rank(colData(sce)$PctZero)]
        colData(sce)$RankDiffZeros <- colData(sce)$PctZero -
            colData(sce)$RefRankZeros

        rowData(sce)$MeanRankVarDiff <- rowData(sce)$VarLogCPM -
            ref.vars.rank[rowData(sce)$RankCounts]
        rowData(sce)$MeanRankZerosDiff <- rowData(sce)$PctZero -
            ref.z.gene.rank[rowData(sce)$RankCounts]

        sces[[name]] <- sce
    }

    ref.sce <- sces[[ref]]
    sces[[ref]] <- NULL

    features <- rowData(sces[[1]])
    cells <- colData(sces[[1]])

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }

    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)

    means <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$RankDiffMeanLogCPM,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "red") +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::ylab(expression(paste(
            "Rank difference mean ", log[2],
            "(CPM + 1)"
        ))) +
        ggplot2::ggtitle("Difference in mean expression") +
        ggplot2::theme_minimal()

    vars <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$RankDiffVarLogCPM,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "red") +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::ylab(expression(paste(
            "Rank difference variance ", log[2],
            "(CPM + 1)"
        ))) +
        ggplot2::ggtitle("Difference in variance") +
        ggplot2::theme_minimal()

    mean.var <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$RankCounts, y = .data$MeanRankVarDiff,
            colour = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "red") +
        ggplot2::geom_point(size = point.size, alpha = point.alpha) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("Expression rank") +
        ggplot2::ylab(expression(paste(
            "Difference in variance ", log[2],
            "(CPM + 1)"
        ))) +
        ggplot2::ggtitle("Difference in mean-variance relationship") +
        ggplot2::theme_minimal()

    libs <- ggplot2::ggplot(
        cells,
        ggplot2::aes(
            x = .data$Dataset, y = .data$RankDiffLibSize,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "red") +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::ylab(paste("Rank difference library size")) +
        ggplot2::ggtitle("Difference in library sizes") +
        ggplot2::theme_minimal()

    z.gene <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$RankDiffZeros,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "red") +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::ylab(paste("Rank difference percentage zeros")) +
        ggplot2::ggtitle("Difference in zeros per gene") +
        ggplot2::theme_minimal()

    z.cell <- ggplot2::ggplot(
        cells,
        ggplot2::aes(
            x = .data$Dataset, y = .data$RankDiffZeros,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "red") +
        ggplot2::geom_boxplot(notch = TRUE, width = 0.1, size = 0.8) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::ylab(paste("Rank difference percentage zeros")) +
        ggplot2::ggtitle("Difference in zeros per cell") +
        ggplot2::theme_minimal()

    mean.zeros <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$RankCounts, y = .data$MeanRankZerosDiff,
            colour = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "red") +
        ggplot2::geom_point(size = point.size, alpha = point.alpha) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("Expression rank") +
        ggplot2::ylab("Difference in percentage zeros per gene") +
        ggplot2::ggtitle("Difference in mean-zeros relationship") +
        ggplot2::theme_minimal()

    means.qq <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$RefRankMeanLogCPM, y = .data$MeanLogCPM,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red") +
        ggplot2::geom_point(size = point.size) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::xlab(expression(paste(
            "Reference mean ", log[2],
            "(CPM + 1)"
        ))) +
        ggplot2::ylab(expression(paste(
            "Alternative mean ", log[2],
            "(CPM + 1)"
        ))) +
        ggplot2::ggtitle("Ranked means") +
        ggplot2::theme_minimal()

    vars.qq <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$RefRankVarLogCPM, y = .data$VarLogCPM,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red") +
        ggplot2::geom_point(size = point.size) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::xlab(expression(paste(
            "Reference variance ", log[2],
            "(CPM + 1)"
        ))) +
        ggplot2::ylab(expression(paste(
            "Alternative variance ", log[2],
            "(CPM + 1)"
        ))) +
        ggplot2::ggtitle("Ranked variances") +
        ggplot2::theme_minimal()

    libs.qq <- ggplot2::ggplot(
        cells,
        ggplot2::aes(
            x = .data$RefRankLibSize, y = .data$sum,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red") +
        ggplot2::geom_point(size = point.size) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::xlab("Reference library size") +
        ggplot2::ylab("Alternative library size") +
        ggplot2::ggtitle("Ranked library sizes") +
        ggplot2::theme_minimal()

    z.gene.qq <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$RefRankZeros, y = .data$PctZero,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red") +
        ggplot2::geom_point(size = point.size) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::xlab("Reference percentage zeros") +
        ggplot2::ylab("Alternative percentage zeros") +
        ggplot2::ggtitle("Ranked percentage zeros per gene") +
        ggplot2::theme_minimal()

    z.cell.qq <- ggplot2::ggplot(
        cells,
        ggplot2::aes(
            x = .data$RefRankZeros, y = .data$PctZero,
            colour = .data$Dataset
        )
    ) +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red") +
        ggplot2::geom_point(size = point.size) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::xlab("Reference percentage zeros") +
        ggplot2::ylab("Alternative percentage zeros") +
        ggplot2::ggtitle("Ranked percentage zeros per cell") +
        ggplot2::theme_minimal()

    if (fits) {
        mean.var <- mean.var +
            ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
        mean.zeros <- mean.zeros +
            ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
    }

    comparison <- list(
        Reference = ref.sce,
        RowData = features,
        ColData = cells,
        Plots = list(
            Means = means,
            Variances = vars,
            MeanVar = mean.var,
            LibrarySizes = libs,
            ZerosGene = z.gene,
            ZerosCell = z.cell,
            MeanZeros = mean.zeros
        ),
        QQPlots = list(
            Means = means.qq,
            Variances = vars.qq,
            LibrarySizes = libs.qq,
            ZerosGene = z.gene.qq,
            ZerosCell = z.cell.qq
        )
    )

    return(comparison)
}

#' Make comparison panel
#'
#' Combine the plots from \code{compareSCEs} into a single panel.
#'
#' @param comp list returned by \code{\link{compareSCEs}}.
#' @param title title for the panel.
#' @param labels vector of labels for each of the seven plots.
#'
#' @return Combined panel plot
#'
#' @examples
#' \dontrun{
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
#' panel <- makeCompPanel(comparison)
#' }
#'
#' @export
makeCompPanel <- function(comp, title = "Comparison",
                          labels = c(
                              "Means", "Variance",
                              "Mean-variance relationship",
                              "Library size", "Zeros per gene",
                              "Zeros per cell",
                              "Mean-zeros relationship"
                          )) {
    checkDependencies(deps = c("ggplot2", "cowplot"))

    checkmate::assertList(comp, any.missing = FALSE, len = 3)
    checkmate::checkString(title)
    checkmate::checkCharacter(labels, len = 7)

    plots <- list(
        p1 = comp$Plots$Means,
        p2 = comp$Plots$Variances,
        p3 = comp$Plots$MeanVar,
        p4 = comp$Plots$LibrarySizes,
        p5 = comp$Plots$ZerosGene,
        p6 = comp$Plots$ZerosCell,
        p7 = comp$Plots$MeanZeros
    )

    # Remove titles and legends
    for (plot in names(plots)) {
        plots[[plot]] <- plots[[plot]] +
            ggplot2::theme(
                legend.position = "none",
                plot.title = ggplot2::element_blank()
            )
    }

    # Remove x-axis title from some plots
    for (plot in paste0("p", c(1, 2, 4, 5, 6))) {
        plots[[plot]] <- plots[[plot]] +
            ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }

    plots$leg <- cowplot::get_legend(
        plots$p3 + ggplot2::theme(legend.position = "bottom")
    )

    panel <- cowplot::ggdraw() +
        cowplot::draw_label(title, 0.5, 0.98, fontface = "bold", size = 18) +
        cowplot::draw_label(
            labels[1], 0.01, 0.95,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p1, 0.0, 0.74, 0.5, 0.20) +
        cowplot::draw_label(
            labels[2], 0.51, 0.95,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p2, 0.5, 0.74, 0.5, 0.20) +
        cowplot::draw_label(
            labels[3], 0.01, 0.70,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p3, 0.0, 0.49, 0.5, 0.20) +
        cowplot::draw_label(
            labels[4], 0.51, 0.70,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p4, 0.5, 0.49, 0.5, 0.20) +
        cowplot::draw_label(
            labels[5], 0.01, 0.45,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p5, 0.0, 0.24, 0.5, 0.20) +
        cowplot::draw_label(
            labels[6], 0.51, 0.45,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p6, 0.5, 0.24, 0.5, 0.20) +
        cowplot::draw_label(
            labels[7], 0.01, 0.21,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p7, 0.0, 0.00, 0.5, 0.20) +
        cowplot::draw_plot(plots$leg, 0.5, 0.00, 0.5, 0.20)

    return(panel)
}


#' Make difference panel
#'
#' Combine the plots from \code{diffSCEs} into a single panel.
#'
#' @param diff list returned by \code{\link{diffSCEs}}.
#' @param title title for the panel.
#' @param labels vector of labels for each of the seven sections.
#'
#' @return Combined panel plot
#'
#' @examples
#' \dontrun{
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' panel <- makeDiffPanel(difference)
#' }
#'
#' @export
makeDiffPanel <- function(diff, title = "Difference comparison",
                          labels = c(
                              "Means", "Variance", "Library size",
                              "Zeros per cell", "Zeros per gene",
                              "Mean-variance relationship",
                              "Mean-zeros relationship"
                          )) {
    checkDependencies(deps = c("ggplot2", "cowplot"))

    checkmate::assertList(diff, any.missing = FALSE, len = 5)
    checkmate::checkString(title)
    checkmate::checkCharacter(labels, len = 7)

    plots <- list(
        p1 = diff$Plots$Means,
        p2 = diff$QQPlots$Means,
        p3 = diff$Plots$Variances,
        p4 = diff$QQPlots$Variances,
        p5 = diff$Plots$MeanVar,
        p6 = diff$Plots$LibrarySizes,
        p7 = diff$QQPlots$LibrarySizes,
        p8 = diff$Plots$ZerosCell,
        p9 = diff$QQPlots$ZerosCell,
        p10 = diff$Plots$ZerosGene,
        p11 = diff$QQPlots$ZerosGene,
        p12 = diff$Plots$MeanZeros
    )

    # Remove titles and legends
    for (plot in names(plots)) {
        plots[[plot]] <- plots[[plot]] +
            ggplot2::theme(
                legend.position = "none",
                plot.title = ggplot2::element_blank()
            )
    }

    # Remove x-axis title from some plots
    for (plot in paste0("p", c(1, 3, 6, 8, 10))) {
        plots[[plot]] <- plots[[plot]] +
            ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }

    plots$leg <- cowplot::get_legend(
        plots$p5 + ggplot2::theme(legend.position = "bottom")
    )

    panel <- cowplot::ggdraw() +
        cowplot::draw_label(
            title, 0.5, 0.98,
            fontface = "bold", size = 18
        ) +
        cowplot::draw_label(
            labels[1], 0.0, 0.94,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p1, 0.0, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p2, 0.0, 0.32, 0.18, 0.29) +
        cowplot::draw_label(
            labels[2], 0.21, 0.94,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p3, 0.21, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p4, 0.21, 0.32, 0.18, 0.29) +
        cowplot::draw_label(
            labels[6], 0.0, 0.30,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p5, 0.0, 0.0, 0.38, 0.29) +
        cowplot::draw_label(
            labels[3], 0.41, 0.94,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p6, 0.41, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p7, 0.41, 0.32, 0.18, 0.29) +
        cowplot::draw_label(
            labels[4], 0.61, 0.94,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p8, 0.61, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p9, 0.61, 0.32, 0.18, 0.29) +
        cowplot::draw_label(
            labels[7], 0.41, 0.30,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p12, 0.41, 0.0, 0.38, 0.29) +
        cowplot::draw_label(
            labels[5], 0.81, 0.94,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p10, 0.81, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p11, 0.81, 0.32, 0.18, 0.29) +
        cowplot::draw_plot(plots$leg, 0.81, 0.0, 0.2, 0.29)

    return(panel)
}


#' Make overall panel
#'
#' Combine the plots from \code{compSCEs} and \code{diffSCEs} into a
#' single panel.
#'
#' @param comp list returned by \code{\link{compareSCEs}}.
#' @param diff list returned by \code{\link{diffSCEs}}.
#' @param title title for the panel.
#' @param row.labels vector of labels for each of the seven rows.
#'
#' @return Combined panel plot
#'
#' @examples
#' \dontrun{
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
#' difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' panel <- makeOverallPanel(comparison, difference)
#' }
#'
#' @export
makeOverallPanel <- function(comp, diff, title = "Overall comparison",
                             row.labels = c(
                                 "Means", "Variance",
                                 "Mean-variance relationship",
                                 "Library size", "Zeros per cell",
                                 "Zeros per gene",
                                 "Mean-zeros relationship"
                             )) {
    checkDependencies(deps = c("ggplot2", "cowplot"))

    checkmate::assertList(comp, any.missing = FALSE, len = 3)
    checkmate::assertList(diff, any.missing = FALSE, len = 5)
    checkmate::checkString(title)
    checkmate::checkCharacter(row.labels, len = 7)

    plots <- list(
        p1 = comp$Plots$Means,
        p2 = diff$Plots$Means,
        p3 = diff$QQPlots$Means,
        p4 = comp$Plots$Variances,
        p5 = diff$Plots$Variances,
        p6 = diff$QQPlots$Variances,
        p7 = comp$Plots$MeanVar,
        p8 = diff$Plots$MeanVar,
        p9 = comp$Plots$LibrarySizes,
        p10 = diff$Plots$LibrarySizes,
        p11 = diff$QQPlots$LibrarySizes,
        p12 = comp$Plots$ZerosCell,
        p13 = diff$Plots$ZerosCell,
        p14 = diff$QQPlots$ZerosCell,
        p15 = comp$Plots$ZerosGene,
        p16 = diff$Plots$ZerosGene,
        p17 = diff$QQPlots$ZerosGene,
        p18 = comp$Plots$MeanZeros,
        p19 = diff$Plots$MeanZeros
    )

    # Remove titles and legends
    for (plot in names(plots)) {
        plots[[plot]] <- plots[[plot]] +
            ggplot2::theme(
                legend.position = "none",
                plot.title = ggplot2::element_blank()
            )
    }

    # Remove x-axis title from some plots
    for (plot in paste0("p", c(1, 2, 4, 5, 9, 10, 12, 13, 15, 16))) {
        plots[[plot]] <- plots[[plot]] +
            ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }

    plots$leg <- cowplot::get_legend(
        plots$p7 + ggplot2::theme(legend.position = "bottom")
    )

    panel <- cowplot::ggdraw() +
        cowplot::draw_label(
            title, 0.5, 0.995,
            fontface = "bold", size = 18
        ) +
        cowplot::draw_label(
            row.labels[1], 0.01, 0.985,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p1, 0.00, 0.86, 0.32, 0.12) +
        cowplot::draw_plot(plots$p2, 0.34, 0.86, 0.32, 0.12) +
        cowplot::draw_plot(plots$p3, 0.67, 0.86, 0.32, 0.12) +
        cowplot::draw_label(
            row.labels[2], 0.01, 0.845,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p4, 0.00, 0.72, 0.32, 0.12) +
        cowplot::draw_plot(plots$p5, 0.34, 0.72, 0.32, 0.12) +
        cowplot::draw_plot(plots$p6, 0.67, 0.72, 0.32, 0.12) +
        cowplot::draw_label(
            row.labels[3], 0.01, 0.705,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p7, 0.00, 0.58, 0.49, 0.12) +
        cowplot::draw_plot(plots$p8, 0.51, 0.58, 0.49, 0.12) +
        cowplot::draw_label(
            row.labels[4], 0.01, 0.56,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p9, 0.00, 0.44, 0.32, 0.12) +
        cowplot::draw_plot(plots$p10, 0.34, 0.44, 0.32, 0.12) +
        cowplot::draw_plot(plots$p11, 0.67, 0.44, 0.32, 0.12) +
        cowplot::draw_label(
            row.labels[5], 0.01, 0.425,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p12, 0.00, 0.30, 0.32, 0.12) +
        cowplot::draw_plot(plots$p13, 0.34, 0.30, 0.32, 0.12) +
        cowplot::draw_plot(plots$p14, 0.67, 0.30, 0.32, 0.12) +
        cowplot::draw_label(
            row.labels[6], 0.01, 0.285,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p15, 0.00, 0.16, 0.32, 0.12) +
        cowplot::draw_plot(plots$p16, 0.34, 0.16, 0.32, 0.12) +
        cowplot::draw_plot(plots$p17, 0.67, 0.16, 0.32, 0.12) +
        cowplot::draw_label(
            row.labels[7], 0.01, 0.145,
            fontface = "bold", hjust = 0, vjust = 0
        ) +
        cowplot::draw_plot(plots$p18, 0.00, 0.02, 0.49, 0.12) +
        cowplot::draw_plot(plots$p19, 0.51, 0.02, 0.49, 0.12) +
        cowplot::draw_plot(plots$leg, 0.00, 0.00, 1.00, 0.02)

    return(panel)
}

#' Summarise diffSCEs
#'
#' Summarise the results of \code{\link{diffSCEs}}. Calculates the Median
#' Absolute Deviation (MAD), Mean Absolute Error (MAE), Root Mean Squared
#' Error (RMSE) and Kolmogorov-Smirnov (KS) statistics for the various
#' properties and ranks them.
#'
#' @param diff Output from \code{\link{diffSCEs}}
#'
#' @return data.frame with MADs, MAEs, RMSEs, scaled statistics and ranks
#' @examples
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' difference <- diffSCEs(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' summary <- summariseDiff(difference)
#' head(summary)
#' @export
#' @importFrom SummarizedExperiment rowData
summariseDiff <- function(diff) {
    row.stats <- c(
        Mean = "RankDiffMeanLogCPM",
        Variance = "RankDiffVarLogCPM",
        ZerosGene = "RankDiffZeros",
        MeanVar = "MeanRankVarDiff",
        MeanZeros = "MeanRankZerosDiff"
    )

    row.ks.stats <- c(
        Mean = "MeanLogCPM",
        Variance = "VarLogCPM",
        ZerosGene = "PctZero",
        MeanVar = NA,
        MeanZeros = NA
    )

    row.mad <- summariseStats(diff$RowData, "Dataset", row.stats, "MAD")
    row.mae <- summariseStats(diff$RowData, "Dataset", row.stats, "MAE")
    row.rmse <- summariseStats(diff$RowData, "Dataset", row.stats, "RMSE")
    row.ks <- summariseKS(
        diff$RowData,
        SummarizedExperiment::rowData(diff$Reference),
        "Dataset",
        row.ks.stats
    )

    row.list <- list(row.mad, row.mae, row.rmse, row.ks)
    row.list <- lapply(row.list, function(summ) {
        summ[, -c(1, 2)]
    })
    row.summ <- data.frame(
        Dataset = row.mad$Dataset,
        Statistic = row.mad$Statistic
    )
    row.list <- c(row.summ, row.list)
    row.summ <- do.call("cbind", row.list)

    col.stats <- c(
        LibSize = "RankDiffLibSize",
        ZerosCell = "RankDiffZeros"
    )

    col.ks.stats <- c(
        LibSize = "sum",
        ZerosCell = "PctZero"
    )

    col.mad <- summariseStats(diff$ColData, "Dataset", col.stats, "MAD")
    col.mae <- summariseStats(diff$ColData, "Dataset", col.stats, "MAE")
    col.rmse <- summariseStats(diff$ColData, "Dataset", col.stats, "RMSE")
    col.ks <- summariseKS(
        diff$ColData,
        SummarizedExperiment::colData(diff$Reference),
        "Dataset",
        col.ks.stats
    )

    col.list <- list(col.mad, col.mae, col.rmse, col.ks)
    col.list <- lapply(col.list, function(summ) {
        summ[, -c(1, 2)]
    })
    col.summ <- data.frame(
        Dataset = col.mad$Dataset,
        Statistic = col.mad$Statistic
    )
    col.list <- c(col.summ, col.list)
    col.summ <- do.call("cbind", col.list)

    summary <- rbind(row.summ, col.summ)

    return(summary)
}

#' Summarise statistics
#'
#' Summarise columns of a data.frame using a single measure.
#'
#' @param data The data.frame to summarise
#' @param split.col Name of the column used to split the dataset
#' @param stat.cols Names of the columns to summarise. If this vector is named
#' those names will be used in the output.
#' @param measure The measure to use for summarisation.
#'
#' @return data.frame with the summarised measure, scaled and ranked
#'
#' @importFrom stats aggregate
summariseStats <- function(data, split.col, stat.cols,
                           measure = c("MAD", "MAE", "RMSE")) {
    measure <- match.arg(measure)

    if (is.null(names(stat.cols))) {
        names(stat.cols) <- stat.cols
    }

    switch(measure,
        "MAD" = {
            measure_fun <- function(x) {
                median(abs(x))
            }
        },
        "MAE" = {
            measure_fun <- function(x) {
                mean(abs(x))
            }
        },
        "RMSE" = {
            measure_fun <- function(x) {
                sqrt(mean(abs(x^2)))
            }
        }
    )

    summ <- aggregate(
        data[, stat.cols],
        list(Dataset = data[[split.col]]),
        measure_fun
    )
    colnames(summ) <- c(split.col, names(stat.cols))

    tidy.summ <- tidyStatSumm(summ, measure)

    return(tidy.summ)
}

#' Summarise KS
#'
#' Summarise columns of a data.frame compared to a reference using the KS test.
#'
#' @param data The data.frame to summarise
#' @param ref The reference data.frame
#' @param split.col Name of the column used to split the dataset
#' @param stat.cols Names of the columns to summarise. If this vector is named
#' those names will be used in the output.
#'
#' @return data.frame with the summarised measure, scaled and ranked
#' @importFrom stats ks.test
summariseKS <- function(data, ref, split.col, stat.cols) {
    if (is.null(names(stat.cols))) {
        names(stat.cols) <- stat.cols
    }

    splits <- unique(data[[split.col]])

    summ <- expand.grid(
        Dataset = splits,
        Statistic = names(stat.cols),
        stringsAsFactors = FALSE
    )

    ks.res <- mapply(function(split, stat.name) {
        stat <- stat.cols[stat.name]
        if (!is.na(stat)) {
            data.stat <- data[data[[split.col]] == split, stat]
            ref.stat <- ref[[stat]]

            ks <- suppressWarnings(ks.test(ref.stat, data.stat))
            ks.out <- c(KS = unname(ks$statistic), KSPVal = ks$p.value)
        } else {
            ks.out <- c(KS = NA, KSPVal = NA)
        }

        return(ks.out)
    }, summ$Dataset, summ$Statistic)

    summ$KS <- ks.res["KS", ]
    summ$KSPVal <- ks.res["KSPVal", ]

    ks.ranks <- lapply(split(summ, summ$Statistic), function(x) {
        rank(x$KS)
    })
    ks.ranks <- unlist(ks.ranks)

    summ$KSRank <- ks.ranks
    summ$KSRank[is.na(summ$KS)] <- NA

    return(summ)
}

#' Tidy summarised statistics
#'
#' Convert a statistic summary to tidy format and add ranks and scaled values
#'
#' @param stat.summ The summary to convert
#' @param measure The name of the summarisation measure
#'
#' @return tidy data.frame with the summarised measure, scaled and ranked
tidyStatSumm <- function(stat.summ, measure = c("MAD", "MAE", "RMSE")) {
    measure <- match.arg(measure)

    summ.mat <- t(stat.summ[, -1])
    colnames(summ.mat) <- stat.summ[, 1]

    scale.summ <- apply(summ.mat, 1, scale)
    # Check if apply has returned a vector
    if (is.vector(scale.summ)) {
        scale.summ <- t(as.matrix(scale.summ))
    }

    rank.summ <- apply(summ.mat, 1, rank)
    if (is.vector(rank.summ)) {
        rank.summ <- t(as.matrix(rank.summ))
    }

    tidy.summ <- as.data.frame.table(t(summ.mat))
    colnames(tidy.summ) <- c("Dataset", "Statistic", measure)

    tidy.scale <- as.data.frame.table(scale.summ)
    tidy.rank <- as.data.frame.table(rank.summ)

    tidy.summ[[paste0(measure, "Scaled")]] <- tidy.scale[, 3]
    tidy.summ[[paste0(measure, "Rank")]] <- tidy.rank[, 3]

    return(tidy.summ)
}
