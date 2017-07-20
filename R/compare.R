#' Compare SCESet objects
#'
#' Combine the data from several SCESet objects and produce some basic plots
#' comparing them.
#'
#' @param sces named list of SCESet objects to combine and compare.
#' @param point.size size of points in scatter plots.
#' @param point.alpha opacity of points in scatter plots.
#' @param fits whether to include fits in scatter plots.
#' @param colours vector of colours to use for each dataset.
#'
#' @details
#' The returned list has three items:
#'
#' \describe{
#'     \item{\code{FeatureData}}{Combined feature data from the provided
#'     SCESets.}
#'     \item{\code{PhenoData}}{Combined pheno data from the provided SCESets.}
#'     \item{\code{Plots}}{Comparison plots
#'         \describe{
#'             \item{\code{Means}}{Boxplot of mean distribution.}
#'             \item{\code{Variances}}{Boxplot of variance distribution.}
#'             \item{\code{MeanVar}}{Scatter plot with fitted lines showing the
#'             mean-variance relationship.}
#'             \item{\code{LibraySizes}}{Boxplot of the library size
#'             distribution.}
#'             \item{\code{ZerosGene}}{Boxplot of the percentage of each gene
#'             that is zero.}
#'             \item{\code{ZerosCell}}{Boxplot of the percentage of each cell
#'             that is zero.}
#'             \item{\code{MeanZeros}}{Scatter plot with fitted lines showing
#'             the mean-zeros relationship.}
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
#' comparison <- compareSCESets(list(Splat = sim1, Simple = sim2))
#' names(comparison)
#' names(comparison$Plots)
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth geom_boxplot
#' scale_y_continuous scale_y_log10 scale_x_log10 xlab ylab ggtitle
#' theme_minimal scale_colour_manual scale_fill_manual
#' @importFrom scater cpm<-
#' @export
compareSCESets <- function(sces, point.size = 0.1, point.alpha = 0.1,
                           fits = TRUE, colours = NULL) {

    checkmate::assertList(sces, types = "SCESet", any.missing = FALSE,
                          min.len = 1, names = "unique")
    checkmate::assertNumber(point.size, finite = TRUE)
    checkmate::assertNumber(point.alpha, lower = 0, upper = 1)
    checkmate::assertLogical(fits, any.missing = FALSE, len = 1)

    if (!is.null(colours)) {
        checkmate::assertCharacter(colours, any.missing = FALSE,
                                   len = length(sces))
    } else {
        colours <- scales::hue_pal()(length(sces))
    }

    for (name in names(sces)) {
        sce <- sces[[name]]
        fData(sce)$Dataset <- name
        pData(sce)$Dataset <- name
        sce <- scater::calculateQCMetrics(sce)
        cpm(sce) <- edgeR::cpm(counts(sce))
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        sces[[name]] <- sce
    }

    fData.all <- fData(sces[[1]])
    pData.all <- pData(sces[[1]])

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            fData.all <- rbindMatched(fData.all, fData(sce))
            pData.all <- rbindMatched(pData.all, pData(sce))
        }
    }

    fData.all$Dataset <- factor(fData.all$Dataset, levels = names(sces))
    pData.all$Dataset <- factor(pData.all$Dataset, levels = names(sces))

    means <- ggplot(fData.all,
                    aes_string(x = "Dataset", y = "MeanLogCPM",
                               colour = "Dataset")) +
        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        geom_boxplot() +
        scale_colour_manual(values = colours) +
        ylab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ggtitle("Distribution of mean expression") +
        theme_minimal()

    vars <- ggplot(fData.all,
                   aes_string(x = "Dataset", y = "VarLogCPM",
                              colour = "Dataset")) +
        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        geom_boxplot() +
        scale_colour_manual(values = colours) +
        ylab(expression(paste("Variance ", log[2], "(CPM + 1)"))) +
        ggtitle("Distribution of variance") +
        theme_minimal()

    mean.var <- ggplot(fData.all,
                       aes_string(x = "MeanLogCPM", y = "VarLogCPM",
                                  colour = "Dataset", fill = "Dataset")) +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        xlab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ylab(expression(paste("Variance ", log[2], "(CPM + 1)"))) +
        ggtitle("Mean-variance relationship") +
        theme_minimal()

    libs <- ggplot(pData.all,
                   aes_string(x = "Dataset", y = "total_counts",
                              colour = "Dataset")) +
        geom_boxplot() +
        scale_y_continuous(labels = scales::comma) +
        scale_colour_manual(values = colours) +
        ylab("Total counts per cell") +
        ggtitle("Distribution of library sizes") +
        theme_minimal()

    z.gene <- ggplot(fData.all,
                     aes_string(x = "Dataset", y = "pct_dropout",
                                colour = "Dataset")) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        scale_colour_manual(values = colours) +
        ylab("Percentage zeros per gene") +
        ggtitle("Distribution of zeros per gene") +
        theme_minimal()

    z.cell <- ggplot(pData.all,
                     aes_string(x = "Dataset", y = "pct_dropout",
                                colour = "Dataset")) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        scale_colour_manual(values = colours) +
        ylab("Percentage zeros per cell") +
        ggtitle("Distribution of zeros per cell") +
        theme_minimal()

    mean.zeros <- ggplot(fData.all,
                         aes_string(x = "MeanCounts", y = "pct_dropout",
                                    colour = "Dataset", fill = "Dataset")) +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_x_log10(labels = scales::comma) +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        xlab("Mean count") +
        ylab("Percentage zeros") +
        ggtitle("Mean-zeros relationship") +
        theme_minimal()

    if (fits) {
        mean.var <- mean.var + geom_smooth()
        mean.zeros <- mean.zeros + geom_smooth()
    }

    comparison <- list(FeatureData = fData.all,
                       PhenoData = pData.all,
                       Plots = list(Means = means,
                                    Variances = vars,
                                    MeanVar = mean.var,
                                    LibrarySizes = libs,
                                    ZerosGene = z.gene,
                                    ZerosCell = z.cell,
                                    MeanZeros = mean.zeros))

    return(comparison)
}

#' Diff SCESet objects
#'
#' Combine the data from several SCESet objects and produce some basic plots
#' comparing them to a reference.
#'
#' @param sces named list of SCESet objects to combine and compare.
#' @param ref string giving the name of the SCESet to use as the reference
#' @param point.size size of points in scatter plots.
#' @param point.alpha opacity of points in scatter plots.
#' @param fits whether to include fits in scatter plots.
#' @param colours vector of colours to use for each dataset.
#'
#' @details
#'
#' This function aims to look at the differences between a reference SCESet and
#' one or more others. It requires each SCESet to have the same dimensions.
#' Properties are compared by ranks, for example when comparing the means the
#' values are ordered and the differences between the reference and another
#' dataset plotted. A series of Q-Q plots are also returned.
#'
#' The returned list has five items:
#'
#' \describe{
#'     \item{\code{Reference}}{The SCESet used as the reference.}
#'     \item{\code{FeatureData}}{Combined feature data from the provided
#'     SCESets.}
#'     \item{\code{PhenoData}}{Combined pheno data from the provided SCESets.}
#'     \item{\code{Plots}}{Difference plots
#'         \describe{
#'             \item{\code{Means}}{Boxplot of mean differences.}
#'             \item{\code{Variances}}{Boxplot of variance differences.}
#'             \item{\code{MeanVar}}{Scatter plot showing the difference from
#'             the reference variance across expression ranks.}
#'             \item{\code{LibraySizes}}{Boxplot of the library size
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
#' difference <- diffSCESets(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' names(difference)
#' names(difference$Plots)
#' @importFrom ggplot2 ggplot aes_string geom_point geom_boxplot xlab ylab
#' ggtitle theme_minimal geom_hline geom_abline scale_colour_manual
#' scale_fill_manual
#' @importFrom scater cpm<-
#' @export
diffSCESets <- function(sces, ref, point.size = 0.1, point.alpha = 0.1,
                        fits = TRUE, colours = NULL) {

    checkmate::assertList(sces, types = "SCESet", any.missing = FALSE,
                          min.len = 2, names = "unique")
    checkmate::assertString(ref)
    checkmate::assertNumber(point.size, finite = TRUE)
    checkmate::assertNumber(point.alpha, lower = 0, upper = 1)
    checkmate::assertLogical(fits, any.missing = FALSE, len = 1)

    if (!(ref %in% names(sces))) {
        stop("'ref' must be the name of an SCESet in 'sces'")
    } else {
        ref.idx <- which(names(sces) == ref)
    }

    if (!is.null(colours)) {
        checkmate::assertCharacter(colours, any.missing = FALSE,
                                   len = length(sces) - 1)
    } else {
        colours <- scales::hue_pal()(length(sces))
        colours <- colours[-ref.idx]
    }

    ref.dim <- dim(sces[[ref]])

    for (name in names(sces)) {
        sce <- sces[[name]]
        if (!identical(dim(sce), ref.dim)) {
            stop("SCESets must have the same dimensions")
        }
        fData(sce)$Dataset <- name
        pData(sce)$Dataset <- name
        sce <- scater::calculateQCMetrics(sce)
        cpm(sce) <- edgeR::cpm(counts(sce))
        sce <- addFeatureStats(sce, "counts")
        sce <- addFeatureStats(sce, "cpm", log = TRUE)
        sces[[name]] <- sce
    }

    ref.sce <- sces[[ref]]

    ref.means <- sort(fData(ref.sce)$MeanLogCPM)
    ref.vars <- sort(fData(ref.sce)$VarLogCPM)
    ref.libs <- sort(pData(ref.sce)$total_counts)
    ref.z.gene <- sort(fData(ref.sce)$pct_dropout)
    ref.z.cell <- sort(pData(ref.sce)$pct_dropout)

    ref.rank.ord <- order(fData(ref.sce)$exprs_rank)
    ref.vars.rank <- fData(ref.sce)$VarLogCPM[ref.rank.ord]
    ref.z.gene.rank <- fData(ref.sce)$pct_dropout[ref.rank.ord]

    for (name in names(sces)) {
        sce <- sces[[name]]
        fData(sce)$RefRankMeanLogCPM <- ref.means[rank(fData(sce)$MeanLogCPM)]
        fData(sce)$RankDiffMeanLogCPM <- fData(sce)$MeanLogCPM -
            fData(sce)$RefRankMeanLogCPM
        fData(sce)$RefRankVarLogCPM <- ref.vars[rank(fData(sce)$VarLogCPM)]
        fData(sce)$RankDiffVarLogCPM <- fData(sce)$VarLogCPM -
            fData(sce)$RefRankVarLogCPM
        pData(sce)$RefRankLibSize <- ref.libs[rank(pData(sce)$total_counts)]
        pData(sce)$RankDiffLibSize <- pData(sce)$total_counts -
            pData(sce)$RefRankLibSize
        fData(sce)$RefRankZeros <- ref.z.gene[rank(fData(sce)$pct_dropout)]
        fData(sce)$RankDiffZeros <- fData(sce)$pct_dropout -
            fData(sce)$RefRankZeros
        pData(sce)$RefRankZeros <- ref.z.cell[rank(pData(sce)$pct_dropout)]
        pData(sce)$RankDiffZeros <- pData(sce)$pct_dropout -
            pData(sce)$RefRankZeros

        fData(sce)$MeanRankVarDiff <- fData(sce)$VarLogCPM -
            ref.vars.rank[fData(sce)$exprs_rank]
        fData(sce)$MeanRankZerosDiff <- fData(sce)$pct_dropout -
            ref.z.gene.rank[fData(sce)$exprs_rank]

        sces[[name]] <- sce
    }

    ref.sce <- sces[[ref]]
    sces[[ref]] <- NULL

    fData.all <- fData(sces[[1]])
    pData.all <- pData(sces[[1]])

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            fData.all <- rbindMatched(fData.all, fData(sce))
            pData.all <- rbindMatched(pData.all, pData(sce))
        }
    }

    fData.all$Dataset <- factor(fData.all$Dataset, levels = names(sces))
    pData.all$Dataset <- factor(pData.all$Dataset, levels = names(sces))

    means <- ggplot(fData.all,
                    aes_string(x = "Dataset", y = "RankDiffMeanLogCPM",
                               colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "red") +
        geom_boxplot() +
        scale_colour_manual(values = colours) +
        ylab(expression(paste("Rank difference mean ", log[2], "(CPM + 1)"))) +
        ggtitle("Difference in mean expression") +
        theme_minimal()

    vars <- ggplot(fData.all,
                    aes_string(x = "Dataset", y = "RankDiffVarLogCPM",
                               colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "red") +
        geom_boxplot() +
        scale_colour_manual(values = colours) +
        ylab(expression(paste("Rank difference variance ", log[2],
                              "(CPM + 1)"))) +
        ggtitle("Difference in variance") +
        theme_minimal()

    mean.var <- ggplot(fData.all,
                       aes_string(x = "exprs_rank", y = "MeanRankVarDiff",
                                  colour = "Dataset", fill = "Dataset")) +
        geom_hline(yintercept = 0, colour = "red") +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        xlab("Expression rank") +
        ylab(expression(paste("Difference in variance ", log[2],
                              "(CPM + 1)"))) +
        ggtitle("Difference in mean-variance relationship") +
        theme_minimal()

    libs <- ggplot(pData.all,
                   aes_string(x = "Dataset", y = "RankDiffLibSize",
                              colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "red") +
        geom_boxplot() +
        scale_colour_manual(values = colours) +
        ylab(paste("Rank difference libray size")) +
        ggtitle("Difference in library sizes") +
        theme_minimal()

    z.gene <- ggplot(fData.all,
                     aes_string(x = "Dataset", y = "RankDiffZeros",
                                colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "red") +
        geom_boxplot() +
        scale_colour_manual(values = colours) +
        ylab(paste("Rank difference percentage zeros")) +
        ggtitle("Difference in zeros per gene") +
        theme_minimal()

    z.cell <- ggplot(pData.all,
                     aes_string(x = "Dataset", y = "RankDiffZeros",
                                colour = "Dataset")) +
        geom_hline(yintercept = 0, colour = "red") +
        geom_boxplot() +
        scale_colour_manual(values = colours) +
        ylab(paste("Rank difference percentage zeros")) +
        ggtitle("Difference in zeros per cell") +
        theme_minimal()

    mean.zeros <- ggplot(fData.all,
                       aes_string(x = "exprs_rank", y = "MeanRankZerosDiff",
                                  colour = "Dataset", fill = "Dataset")) +
        geom_hline(yintercept = 0, colour = "red") +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        scale_fill_manual(values = colours) +
        xlab("Expression rank") +
        ylab("Difference in percentage zeros per gene") +
        ggtitle("Difference in mean-zeros relationship") +
        theme_minimal()

    means.qq <- ggplot(fData.all,
                       aes_string(x = "RefRankMeanLogCPM", y = "MeanLogCPM",
                                  colour = "Dataset")) +
        geom_abline(intercept = 0, slope = 1, colour = "red") +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        xlab(expression(paste("Reference mean ", log[2], "(CPM + 1)"))) +
        ylab(expression(paste("Alternative mean ", log[2], "(CPM + 1)"))) +
        ggtitle("Ranked means") +
        theme_minimal()

    vars.qq <- ggplot(fData.all,
                      aes_string(x = "RefRankVarLogCPM", y = "VarLogCPM",
                                 colour = "Dataset")) +
        geom_abline(intercept = 0, slope = 1, colour = "red") +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        xlab(expression(paste("Reference variance ", log[2], "(CPM + 1)"))) +
        ylab(expression(paste("Alternative variance ", log[2], "(CPM + 1)"))) +
        ggtitle("Ranked variances") +
        theme_minimal()

    libs.qq <- ggplot(pData.all,
                      aes_string(x = "RefRankLibSize", y = "total_counts",
                                 colour = "Dataset")) +
        geom_abline(intercept = 0, slope = 1, colour = "red") +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        xlab("Reference library size") +
        ylab("Alternative library size") +
        ggtitle("Ranked library sizes") +
        theme_minimal()

    z.gene.qq <- ggplot(fData.all,
                        aes_string(x = "RefRankZeros", y = "pct_dropout",
                                   colour = "Dataset")) +
        geom_abline(intercept = 0, slope = 1, colour = "red") +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        xlab("Reference percentage zeros") +
        ylab("Alternative percentage zeros") +
        ggtitle("Ranked percentage zeros per gene") +
        theme_minimal()

    z.cell.qq <- ggplot(pData.all,
                        aes_string(x = "RefRankZeros", y = "pct_dropout",
                                   colour = "Dataset")) +
        geom_abline(intercept = 0, slope = 1, colour = "red") +
        geom_point(size = point.size, alpha = point.alpha) +
        scale_colour_manual(values = colours) +
        xlab("Reference percentage zeros") +
        ylab("Alternative percentage zeros") +
        ggtitle("Ranked percentage zeros per cell") +
        theme_minimal()

    if (fits) {
        mean.var <- mean.var + geom_smooth()
        mean.zeros <- mean.zeros + geom_smooth()
    }

    comparison <- list(Reference = ref.sce,
                       FeatureData = fData.all,
                       PhenoData = pData.all,
                       Plots = list(Means = means,
                                    Variances = vars,
                                    MeanVar = mean.var,
                                    LibrarySizes = libs,
                                    ZerosGene = z.gene,
                                    ZerosCell = z.cell,
                                    MeanZeros = mean.zeros),
                       QQPlots = list(Means = means.qq,
                                      Variances = vars.qq,
                                      LibrarySizes = libs.qq,
                                      ZerosGene = z.gene.qq,
                                      ZerosCell = z.cell.qq))

    return(comparison)
}

#' Make comparison panel
#'
#' Combine the plots from \code{compareSCESets} into a single panel.
#'
#' @param comp list returned by \code{\link{compareSCESets}}.
#' @param title title for the panel.
#' @param labels vector of labels for each of the seven plots.
#'
#' @return Combined panel plot
#'
#' @examples
#' \dontrun{
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' comparison <- compareSCESets(list(Splat = sim1, Simple = sim2))
#' panel <- makeCompPanel(comparison)
#' }
#'
#' @importFrom ggplot2 theme element_blank
#' @export
makeCompPanel <- function(comp, title = "Comparison",
                          labels = c("Means", "Variance",
                                     "Mean-variance relationship",
                                     "Library size", "Zeros per gene",
                                     "Zeros per cell",
                                     "Mean-zeros relationship")) {

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The `cowplot` package is required to make panels.")
    }

    checkmate::assertList(comp, any.missing = FALSE, len = 3)
    checkmate::checkString(title)
    checkmate::checkCharacter(labels, len = 7)

    plots <- list(p1 = comp$Plots$Means, p2 = comp$Plots$Variances,
                  p3 = comp$Plots$MeanVar, p4 = comp$Plots$LibrarySizes,
                  p5 = comp$Plots$ZerosGene, p6 = comp$Plots$ZerosCell,
                  p7 = comp$Plots$MeanZeros)

    # Remove titles and legends
    for (plot in names(plots)) {
        plots[[plot]] <- plots[[plot]] + theme(legend.position = "none",
                                               plot.title = element_blank())
    }

    # Remove x-axis title from some plots
    for (plot in paste0("p", c(1, 2, 4, 5, 6))) {
        plots[[plot]] <- plots[[plot]] + theme(axis.title.x = element_blank())
    }

    plots$leg <- cowplot::get_legend(plots$p3 +
                                         theme(legend.position = "bottom"))

    panel <- cowplot::ggdraw() +
        cowplot::draw_label(title, 0.5, 0.98,
                            fontface = "bold", size = 18) +
        cowplot::draw_label(labels[1], 0.01, 0.95,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p1,  0.0, 0.74, 0.5, 0.20) +
        cowplot::draw_label(labels[2], 0.51, 0.95,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p2,  0.5, 0.74, 0.5, 0.20) +
        cowplot::draw_label(labels[3], 0.01, 0.70,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p3,  0.0, 0.49, 0.5, 0.20) +
        cowplot::draw_label(labels[4], 0.51, 0.70,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p4,  0.5, 0.49, 0.5, 0.20) +
        cowplot::draw_label(labels[5], 0.01, 0.45,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p5,  0.0, 0.24, 0.5, 0.20) +
        cowplot::draw_label(labels[6], 0.51, 0.45,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p6,  0.5, 0.24, 0.5, 0.20) +
        cowplot::draw_label(labels[7], 0.01, 0.21,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p7,  0.0, 0.00, 0.5, 0.20) +
        cowplot::draw_plot(plots$leg, 0.5, 0.00, 0.5, 0.20)

    return(panel)
}


#' Make difference panel
#'
#' Combine the plots from \code{diffSCESets} into a single panel.
#'
#' @param diff list returned by \code{\link{diffSCESets}}.
#' @param title title for the panel.
#' @param labels vector of labels for each of the seven sections.
#'
#' @return Combined panel plot
#'
#' @examples
#' \dontrun{
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' difference <- diffSCESets(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' panel <- makeDiffPanel(difference)
#' }
#'
#' @importFrom ggplot2 theme element_blank
#' @export
makeDiffPanel <- function(diff, title = "Difference comparison",
                          labels = c("Means", "Variance", "Library size",
                                     "Zeros per cell", "Zeros per gene",
                                     "Mean-variance relationship",
                                     "Mean-zeros relationship")) {

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The `cowplot` package is required to make panels.")
    }

    checkmate::assertList(diff, any.missing = FALSE, len = 5)
    checkmate::checkString(title)
    checkmate::checkCharacter(labels, len = 7)

    plots <- list(p1 = diff$Plots$Means, p2 = diff$QQPlots$Means,
                  p3 = diff$Plots$Variances, p4 = diff$QQPlots$Variances,
                  p5 = diff$Plots$MeanVar, p6 = diff$Plots$LibrarySizes,
                  p7 = diff$QQPlots$LibrarySizes, p8 = diff$Plots$ZerosCell,
                  p9 = diff$QQPlots$ZerosCell, p10 = diff$Plots$ZerosGene,
                  p11 = diff$QQPlots$ZerosGene, p12 = diff$Plots$MeanZeros)

    # Remove titles and legends
    for (plot in names(plots)) {
        plots[[plot]] <- plots[[plot]] + theme(legend.position = "none",
                                               plot.title = element_blank())
    }

    # Remove x-axis title from some plots
    for (plot in paste0("p", c(1, 3, 6, 8, 10))) {
        plots[[plot]] <- plots[[plot]] + theme(axis.title.x = element_blank())
    }

    plots$leg <- cowplot::get_legend(plots$p5 +
                                         theme(legend.position = "bottom"))

    panel <- cowplot::ggdraw() +
        cowplot::draw_label(title, 0.5, 0.98,
                            fontface = "bold", size = 18) +
        cowplot::draw_label(labels[1], 0.0, 0.94,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p1,  0.0, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p2,  0.0, 0.32, 0.18, 0.29) +
        cowplot::draw_label(labels[2], 0.21, 0.94,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p3,  0.21, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p4,  0.21, 0.32, 0.18, 0.29) +
        cowplot::draw_label(labels[6], 0.0, 0.30,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p5,  0.0, 0.0, 0.38, 0.29) +
        cowplot::draw_label(labels[3], 0.41, 0.94,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p6,  0.41, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p7,  0.41, 0.32, 0.18, 0.29) +
        cowplot::draw_label(labels[4], 0.61, 0.94,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p8,  0.61, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p9,  0.61, 0.32, 0.18, 0.29) +
        cowplot::draw_label(labels[7], 0.41, 0.30,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p12, 0.41, 0.0, 0.38, 0.29) +
        cowplot::draw_label(labels[5], 0.81, 0.94,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p10, 0.81, 0.64, 0.18, 0.29) +
        cowplot::draw_plot(plots$p11, 0.81, 0.32, 0.18, 0.29) +
        cowplot::draw_plot(plots$leg, 0.81, 0.0, 0.2, 0.29)

    return(panel)
}


#' Make overall panel
#'
#' Combine the plots from \code{compSCESets} and \code{diffSCESets} into a
#' single panel.
#'
#' @param comp list returned by \code{\link{compareSCESets}}.
#' @param diff list returned by \code{\link{diffSCESets}}.
#' @param title title for the panel.
#' @param row.labels vector of labels for each of the seven rows.
#'
#' @return Combined panel plot
#'
#' @examples
#' \dontrun{
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' comparison <- compSCESets(list(Splat = sim1, Simple = sim2))
#' difference <- diffSCESets(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' panel <- makeOverallPanel(comparison, difference)
#' }
#'
#' @importFrom ggplot2 theme element_blank
#' @export
makeOverallPanel <- function(comp, diff, title = "Overall comparison",
                             row.labels = c("Means", "Variance",
                                            "Mean-variance relationship",
                                            "Library size", "Zeros per cell",
                                            "Zeros per gene",
                                            "Mean-zeros relationship")) {

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The `cowplot` package is required to make panels.")
    }

    checkmate::assertList(comp, any.missing = FALSE, len = 3)
    checkmate::assertList(diff, any.missing = FALSE, len = 5)
    checkmate::checkString(title)
    checkmate::checkCharacter(row.labels, len = 7)

    plots <- list(p1 = comp$Plots$Means, p2 = diff$Plots$Means,
                  p3 = diff$QQPlots$Means, p4 = comp$Plots$Variances,
                  p5 = diff$Plots$Variances, p6 = diff$QQPlots$Variances,
                  p7 = comp$Plots$MeanVar, p8 = diff$Plots$MeanVar,
                  p9 = comp$Plots$LibrarySizes, p10 = diff$Plots$LibrarySizes,
                  p11 = diff$QQPlots$LibrarySizes, p12 = comp$Plots$ZerosCell,
                  p13 = diff$Plots$ZerosCell, p14 = diff$QQPlots$ZerosCell,
                  p15 = comp$Plots$ZerosGene, p16 = diff$Plots$ZerosGene,
                  p17 = diff$QQPlots$ZerosGene, p18 = comp$Plots$MeanZeros,
                  p19 = diff$Plots$MeanZeros)

    # Remove titles and legends
    for (plot in names(plots)) {
        plots[[plot]] <- plots[[plot]] + theme(legend.position = "none",
                                               plot.title = element_blank())
    }

    # Remove x-axis title from some plots
    for (plot in paste0("p", c(1, 2, 4, 5, 9, 10, 12, 13, 15, 16))) {
        plots[[plot]] <- plots[[plot]] + theme(axis.title.x = element_blank())
    }

    plots$leg <- cowplot::get_legend(plots$p7 +
                                         theme(legend.position = "bottom"))

    panel <- cowplot::ggdraw() +
        cowplot::draw_label(title, 0.5, 0.995,
                            fontface = "bold", size = 18) +
        cowplot::draw_label(row.labels[1], 0.01, 0.985,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p1,  0.00, 0.86, 0.32, 0.12) +
        cowplot::draw_plot(plots$p2,  0.34, 0.86, 0.32, 0.12) +
        cowplot::draw_plot(plots$p3,  0.67, 0.86, 0.32, 0.12) +
        cowplot::draw_label(row.labels[2], 0.01, 0.845,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p4,  0.00, 0.72, 0.32, 0.12) +
        cowplot::draw_plot(plots$p5,  0.34, 0.72, 0.32, 0.12) +
        cowplot::draw_plot(plots$p6,  0.67, 0.72, 0.32, 0.12) +
        cowplot::draw_label(row.labels[3], 0.01, 0.705,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p7,  0.00, 0.58, 0.49, 0.12) +
        cowplot::draw_plot(plots$p8,  0.51, 0.58, 0.49, 0.12) +
        cowplot::draw_label(row.labels[4], 0.01, 0.56,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p9,  0.00, 0.44, 0.32, 0.12) +
        cowplot::draw_plot(plots$p10, 0.34, 0.44, 0.32, 0.12) +
        cowplot::draw_plot(plots$p11, 0.67, 0.44, 0.32, 0.12) +
        cowplot::draw_label(row.labels[5], 0.01, 0.425,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p12, 0.00, 0.30, 0.32, 0.12) +
        cowplot::draw_plot(plots$p13, 0.34, 0.30, 0.32, 0.12) +
        cowplot::draw_plot(plots$p14, 0.67, 0.30, 0.32, 0.12) +
        cowplot::draw_label(row.labels[6], 0.01, 0.285,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p15, 0.00, 0.16, 0.32, 0.12) +
        cowplot::draw_plot(plots$p16, 0.34, 0.16, 0.32, 0.12) +
        cowplot::draw_plot(plots$p17, 0.67, 0.16, 0.32, 0.12) +
        cowplot::draw_label(row.labels[7], 0.01, 0.145,
                            fontface = "bold", hjust = 0, vjust = 0) +
        cowplot::draw_plot(plots$p18, 0.00, 0.02, 0.49, 0.12) +
        cowplot::draw_plot(plots$p19, 0.51, 0.02, 0.49, 0.12) +
        cowplot::draw_plot(plots$leg, 0.00, 0.00, 1.00, 0.02)

    return(panel)
}

#' Summarise diffSCESets
#'
#' Summarise the results of \code{\link{diffSCESets}}. Calculates the Median
#' Absolute Deviation (MAD), Mean Absolute Error (MAE) and Root Mean Squared
#' Error (RMSE) for the various properties and ranks them.
#'
#' @param diff Output from \code{\link{diffSCESets}}
#'
#' @return data.frame with MADs, MAEs, RMSEs, scaled statistics and ranks
#' @examples
#' sim1 <- splatSimulate(nGenes = 1000, batchCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' difference <- diffSCESets(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' summary <- summariseDiff(difference)
#' head(summary)
#' @export
summariseDiff <- function(diff) {

    datasets <- unique(diff$PhenoData$Dataset)

    fData.mads <- sapply(datasets, function(dataset) {
        df <- diff$FeatureData[diff$FeatureData$Dataset == dataset, ]
        mean <- median(abs(df$RankDiffMeanLogCPM))
        var <- median(abs(df$RankDiffVarLogCPM))
        zeros <- median(abs(df$RankDiffZeros))
        mean.var <- median(abs(df$MeanRankVarDiff))
        mean.zeros <- median(abs(df$MeanRankZerosDiff))
        return(c(Mean = mean, Variance = var, ZerosGene = zeros,
                 MeanVar = mean.var, MeanZeros = mean.zeros))
    })
    fData.mads.z <- t(scale(t(fData.mads)))

    pData.mads <- sapply(datasets, function(dataset) {
        df <- diff$PhenoData[diff$PhenoData$Dataset == dataset, ]
        lib.size <- median(abs(df$RankDiffLibSize))
        zeros <- median(abs(df$RankDiffZeros))
        return(c(LibSize = lib.size, ZerosCell = zeros))
    })
    pData.mads.z <- t(scale(t(pData.mads)))

    mads <- data.frame(Dataset = datasets, t(fData.mads), t(pData.mads))
    mads.z <- data.frame(Dataset = datasets, t(fData.mads.z), t(pData.mads.z))

    fData.ranks <- matrixStats::rowRanks(fData.mads)
    pData.ranks <- matrixStats::rowRanks(pData.mads)

    ranks.mads <- data.frame(Dataset = datasets, t(fData.ranks), t(pData.ranks))
    colnames(ranks.mads) <- paste0(colnames(mads), "Rank")

    fData.maes <- sapply(datasets, function(dataset) {
        df <- diff$FeatureData[diff$FeatureData$Dataset == dataset, ]
        mean <- mean(abs(df$RankDiffMeanLogCPM))
        var <- mean(abs(df$RankDiffVarLogCPM))
        zeros <- mean(abs(df$RankDiffZeros))
        mean.var <- mean(abs(df$MeanRankVarDiff))
        mean.zeros <- mean(abs(df$MeanRankZerosDiff))
        return(c(Mean = mean, Variance = var, ZerosGene = zeros,
                 MeanVar = mean.var, MeanZeros = mean.zeros))
    })
    fData.maes.z <- t(scale(t(fData.maes)))

    pData.maes <- sapply(datasets, function(dataset) {
        df <- diff$PhenoData[diff$PhenoData$Dataset == dataset, ]
        lib.size <- mean(abs(df$RankDiffLibSize))
        zeros <- mean(abs(df$RankDiffZeros))
        return(c(LibSize = lib.size, ZerosCell = zeros))
    })
    pData.maes.z <- t(scale(t(pData.maes)))

    maes <- data.frame(Dataset = datasets, t(fData.maes), t(pData.maes))
    maes.z <- data.frame(Dataset = datasets, t(fData.maes.z), t(pData.maes.z))

    fData.ranks <- matrixStats::rowRanks(fData.maes)
    pData.ranks <- matrixStats::rowRanks(pData.maes)

    ranks.maes <- data.frame(Dataset = datasets, t(fData.ranks), t(pData.ranks))
    colnames(ranks.maes) <- paste0(colnames(mads), "Rank")

    fData.rmse <- sapply(datasets, function(dataset) {
        df <- diff$FeatureData[diff$FeatureData$Dataset == dataset, ]
        mean <- sqrt(mean(df$RankDiffMeanLogCPM ^ 2))
        var <- sqrt(mean(df$RankDiffVarLogCPM ^ 2))
        zeros <- sqrt(mean(df$RankDiffZeros ^ 2))
        mean.var <- sqrt(mean(df$MeanRankVarDiff ^ 2))
        mean.zeros <- sqrt(mean(df$MeanRankZerosDiff ^ 2))
        return(c(Mean = mean, Variance = var, ZerosGene = zeros,
                 MeanVar = mean.var, MeanZeros = mean.zeros))
    })
    fData.rmse.z <- t(scale(t(fData.rmse)))

    pData.rmse <- sapply(datasets, function(dataset) {
        df <- diff$PhenoData[diff$PhenoData$Dataset == dataset, ]
        lib.size <- sqrt(mean(df$RankDiffLibSize ^ 2))
        zeros <- sqrt(mean(df$RankDiffZeros ^ 2))
        return(c(LibSize = lib.size, ZerosCell = zeros))
    })
    pData.rmse.z <- t(scale(t(pData.rmse)))

    rmse <- data.frame(Dataset = datasets, t(fData.rmse), t(pData.rmse))
    rmse.z <- data.frame(Dataset = datasets, t(fData.rmse.z), t(pData.rmse.z))

    fData.ranks <- matrixStats::rowRanks(fData.rmse)
    pData.ranks <- matrixStats::rowRanks(pData.rmse)

    ranks.rmse <- data.frame(Dataset = datasets, t(fData.ranks), t(pData.ranks))
    colnames(ranks.rmse) <- paste0(colnames(rmse), "Rank")

    mads <- stats::reshape(mads, varying = 2:8, direction = "long",
                           idvar = "Dataset", timevar = "Statistic",
                           times = colnames(mads)[2:8], v.names = "MAD")

    mads.z <- stats::reshape(mads.z, varying = 2:8, direction = "long",
                             idvar = "Dataset", timevar = "Statistic",
                             times = colnames(mads)[2:8],
                             v.names = "MADScaled")

    ranks.mads <- stats::reshape(ranks.mads, varying = 2:8, direction = "long",
                                 idvar = "Dataset", timevar = "Statistic",
                                 times = colnames(ranks.mads)[2:8],
                                 v.names = "Rank")

    maes <- stats::reshape(maes, varying = 2:8, direction = "long",
                           idvar = "Dataset", timevar = "Statistic",
                           times = colnames(maes)[2:8], v.names = "MAE")

    maes.z <- stats::reshape(maes.z, varying = 2:8, direction = "long",
                             idvar = "Dataset", timevar = "Statistic",
                             times = colnames(mads)[2:8],
                             v.names = "MAEScaled")

    ranks.maes <- stats::reshape(ranks.maes, varying = 2:8, direction = "long",
                                 idvar = "Dataset", timevar = "Statistic",
                                 times = colnames(ranks.maes)[2:8],
                                 v.names = "Rank")

    rmse <- stats::reshape(rmse, varying = 2:8, direction = "long",
                           idvar = "Dataset", timevar = "Statistic",
                           times = colnames(mads)[2:8], v.names = "RMSE")

    rmse.z <- stats::reshape(rmse.z, varying = 2:8, direction = "long",
                             idvar = "Dataset", timevar = "Statistic",
                             times = colnames(mads)[2:8],
                             v.names = "RMSEScaled")

    ranks.rmse <- stats::reshape(ranks.rmse, varying = 2:8, direction = "long",
                                 idvar = "Dataset", timevar = "Statistic",
                                 times = colnames(ranks.rmse)[2:8],
                                 v.names = "Rank")

    summary <- data.frame(mads, MADScaled = mads.z$MADScaled,
                          MADRank = ranks.mads$Rank,
                          MAE = maes$MAE, MAEScaled = maes.z$MAEScaled,
                          MAERank = ranks.maes$Rank,
                          RMSE = rmse$RMSE, RMSEScaled = rmse.z$RMSEScaled,
                          RMSERank = ranks.rmse$Rank)
    row.names(summary) <- NULL

    return(summary)
}
