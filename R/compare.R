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
#'             the mean-dropout relationship.}
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
#' sim1 <- splatSimulate(nGenes = 1000, groupCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' comparison <- compareSCESets(list(Splat = sim1, Simple = sim2))
#' names(comparison)
#' names(comparison$Plots)
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth geom_boxplot
#' scale_y_continuous scale_y_log10 scale_x_log10 xlab ylab ggtitle
#' theme_minimal
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
                   aes_string(x = "Dataset", y = "VarCPM",
                              colour = "Dataset")) +
        #geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        geom_boxplot() +
        scale_y_log10(labels = scales::comma) +
        scale_colour_manual(values = colours) +
        ylab("CPM Variance") +
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
        ggtitle("Mean-dropout relationship") +
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
#' sim1 <- splatSimulate(nGenes = 1000, groupCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' difference <- diffSCESets(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' names(difference)
#' names(difference$Plots)
#' @importFrom ggplot2 ggplot aes_string geom_point geom_boxplot xlab ylab
#' ggtitle theme_minimal geom_hline geom_abline
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
    }

    if (!is.null(colours)) {
        checkmate::assertCharacter(colours, any.missing = FALSE,
                                   len = length(sces) - 1)
    } else {
        colours <- scales::hue_pal()(length(sces))
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
#' Combine the plots from \code{combineSCESets} into a single panel.
#'
#' @param comp list returned by \code{\link{combineSCESets}}.
#'
#' @return Combined panel plot
#'
#' @examples
#' sim1 <- splatSimulate(nGenes = 1000, groupCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' comparison <- compareSCESets(list(Splat = sim1, Simple = sim2))
#' panel <- makeCompPanel(comparison)
#'
#' @importFrom ggplot2 theme
#' @export
makeCompPanel <- function(comp,
                          labels = c("A", "B", "C", "D", "E", "F", "G", "")) {

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The `cowplot` package is required to make panels.")
    }

    p1 <- comp$Plots$Means + theme(legend.position = "none",
                                   axis.title.x = element_blank())
    p2 <- comp$Plots$Variances + theme(legend.position = "none",
                                       axis.title.x = element_blank())
    p3 <- comp$Plots$MeanVar + theme(legend.position = "none")
    p4 <- comp$Plots$LibrarySizes + theme(legend.position = "none",
                                          axis.title.x = element_blank())
    p5 <- comp$Plots$ZerosGene + theme(legend.position = "none",
                                       axis.title.x = element_blank())
    p6 <-  comp$Plots$ZerosCell + theme(legend.position = "none",
                                        axis.title.x = element_blank())
    p7 <- comp$Plots$MeanZeros + theme(legend.position = "none")
    leg <- cowplot::get_legend(p1 + theme(legend.position = "bottom"))

    #panel <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, leg,
    #                            ncol = 2,
    #                            labels = c("A", "B", "C", "D", "E", "F", "G", ""))

    panel <- cowplot::ggdraw() +
        cowplot::draw_label(labels[1], 0.01, 0.97, fontface = "bold") +
        cowplot::draw_plot(p1,  0.0, 0.75, 0.5, 0.20) +
        cowplot::draw_label(labels[2], 0.51, 0.97, fontface = "bold") +
        cowplot::draw_plot(p2,  0.5, 0.75, 0.5, 0.20) +
        cowplot::draw_label(labels[3], 0.01, 0.72, fontface = "bold") +
        cowplot::draw_plot(p3,  0.0, 0.50, 0.5, 0.20) +
        cowplot::draw_label(labels[4], 0.51, 0.72, fontface = "bold") +
        cowplot::draw_plot(p4,  0.5, 0.50, 0.5, 0.20) +
        cowplot::draw_label(labels[5], 0.01, 0.47, fontface = "bold") +
        cowplot::draw_plot(p5,  0.0, 0.25, 0.5, 0.20) +
        cowplot::draw_label(labels[6], 0.51, 0.47, fontface = "bold") +
        cowplot::draw_plot(p6,  0.5, 0.25, 0.5, 0.20) +
        cowplot::draw_label(labels[7], 0.01, 0.22, fontface = "bold") +
        cowplot::draw_plot(p7,  0.0, 0.00, 0.5, 0.20) +
        cowplot::draw_label(labels[8], 0.51, 0.22, fontface = "bold") +
        cowplot::draw_plot(leg, 0.5, 0.00, 0.5, 0.20)

    return(panel)
}


#' Make difference panel
#'
#' Combine the plots from \code{diffSCESets} into a single panel.
#'
#' @param comp list returned by \code{\link{diffSCESets}}.
#'
#' @return Combined panel plot
#'
#' @examples
#' sim1 <- splatSimulate(nGenes = 1000, groupCells = 20)
#' sim2 <- simpleSimulate(nGenes = 1000, nCells = 20)
#' difference <- diffSCESets(list(Splat = sim1, Simple = sim2), ref = "Simple")
#' panel <- makeCompPanel(difference)
#'
#' @importFrom ggplot2 theme
#' @export
makeDiffPanel <- function(diff,
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H",
                                     "I", "J", "K", "L", "")) {

    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("The `cowplot` package is required to make panels.")
    }

    p1 <- diff$Plots$Means + theme(legend.position = "none",
                                   axis.title.x = element_blank())
    p2 <- diff$QQPlots$Means + theme(legend.position = "none")
    p3 <- diff$Plots$Variances + theme(legend.position = "none",
                                       axis.title.x = element_blank())
    p4 <- diff$QQPlots$Variances + theme(legend.position = "none")
    p5 <- diff$Plots$MeanVar + theme(legend.position = "none")
    p6 <- diff$Plots$LibrarySizes + theme(legend.position = "none",
                                          axis.title.x = element_blank())
    p7 <- diff$QQPlots$LibrarySizes + theme(legend.position = "none")
    p8 <- diff$Plots$ZerosCell + theme(legend.position = "none",
                                       axis.title.x = element_blank())
    p9 <- diff$QQPlots$ZerosCell + theme(legend.position = "none")
    p10 <- diff$Plots$ZerosGene + theme(legend.position = "none",
                                       axis.title.x = element_blank())
    p11 <- diff$QQPlots$ZerosGene + theme(legend.position = "none")
    p12 <- diff$Plots$MeanZeros + theme(legend.position = "none")

    leg <- cowplot::get_legend(p1 + theme(legend.position = "bottom"))

#    panel <- cowplot::ggdraw() +
#        #cowplot::draw_label(labels[1], 0.01, 0.97, fontface = "bold") +
#        cowplot::draw_plot(p1,  0.0, 0.86, 0.5, 0.14) +
#        #cowplot::draw_label(labels[2], 0.51, 0.97, fontface = "bold") +
#        cowplot::draw_plot(p2,  0.5, 0.86, 0.5, 0.14) +
#        #cowplot::draw_label(labels[3], 0.01, 0.72, fontface = "bold") +
#        cowplot::draw_plot(p3,  0.0, 0.72, 0.5, 0.14) +
#        #cowplot::draw_label(labels[4], 0.51, 0.72, fontface = "bold") +
#        cowplot::draw_plot(p4,  0.5, 0.72, 0.5, 0.14) +
#        #cowplot::draw_label(labels[5], 0.01, 0.47, fontface = "bold") +
#        cowplot::draw_plot(p5,  0.0, 0.58, 1.0, 0.14) +
#        #cowplot::draw_label(labels[6], 0.51, 0.47, fontface = "bold") +
#        cowplot::draw_plot(p6,  0.0, 0.44, 0.5, 0.14) +
#        #cowplot::draw_label(labels[7], 0.01, 0.22, fontface = "bold") +
#        cowplot::draw_plot(p7,  0.5, 0.44, 0.5, 0.14) +
#        #cowplot::draw_label(labels[8], 0.51, 0.22, fontface = "bold") +
#        cowplot::draw_plot(p8,  0.0, 0.30, 0.5, 0.14) +
#        #cowplot::draw_label(labels[9], 0.51, 0.22, fontface = "bold") +
#        cowplot::draw_plot(p9,  0.5, 0.30, 0.5, 0.14) +
#        #cowplot::draw_label(labels[10], 0.51, 0.22, fontface = "bold") +
#        cowplot::draw_plot(p10, 0.0, 0.16, 0.5, 0.14) +
#        #cowplot::draw_label(labels[11], 0.51, 0.22, fontface = "bold") +
#        cowplot::draw_plot(p11, 0.5, 0.16, 0.5, 0.14) +
#        #cowplot::draw_label(labels[12], 0.51, 0.22, fontface = "bold") +
#        cowplot::draw_plot(p12, 0.0, 0.02, 1.0, 0.14) +
#        #cowplot::draw_label(labels[13], 0.51, 0.22, fontface = "bold") +
#        cowplot::draw_plot(leg, 0.0, 0.00, 1.0, 0.02)

#    panel <- cowplot::ggdraw() +
#            #cowplot::draw_label(labels[1], 0.01, 0.97, fontface = "bold") +
#            cowplot::draw_plot(p1,  0.0, 0.8, 0.3, 0.20) +
#            #cowplot::draw_label(labels[2], 0.51, 0.97, fontface = "bold") +
#            cowplot::draw_plot(p2,  0.3, 0.8, 0.3, 0.20) +
#            #cowplot::draw_label(labels[3], 0.01, 0.72, fontface = "bold") +
#            cowplot::draw_plot(p3,  0.0, 0.6, 0.3, 0.20) +
#            #cowplot::draw_label(labels[4], 0.51, 0.72, fontface = "bold") +
#            cowplot::draw_plot(p4,  0.3, 0.6, 0.3, 0.20) +
#            #cowplot::draw_label(labels[5], 0.01, 0.47, fontface = "bold") +
#            cowplot::draw_plot(p5,  0.6, 0.6, 0.4, 0.40) +
#            #cowplot::draw_label(labels[6], 0.51, 0.47, fontface = "bold") +
#            cowplot::draw_plot(p6,  0.0, 0.4, 0.3, 0.20) +
#           #cowplot::draw_label(labels[7], 0.01, 0.22, fontface = "bold") +
#            cowplot::draw_plot(p7,  0.3, 0.4, 0.3, 0.20) +
#            #cowplot::draw_label(labels[8], 0.51, 0.22, fontface = "bold") +
#            cowplot::draw_plot(p8,  0.0, 0.2, 0.3, 0.20) +
#            #cowplot::draw_label(labels[9], 0.51, 0.22, fontface = "bold") +
#            cowplot::draw_plot(p9,  0.3, 0.2, 0.3, 0.20) +
#            #cowplot::draw_label(labels[10], 0.51, 0.22, fontface = "bold") +
#            cowplot::draw_plot(p12, 0.6, 0.2, 0.4, 0.40) +
#            #cowplot::draw_label(labels[11], 0.51, 0.22, fontface = "bold") +
#            cowplot::draw_plot(p10, 0.0, 0.0, 0.3, 0.20) +
#            #cowplot::draw_label(labels[12], 0.51, 0.22, fontface = "bold") +
#            cowplot::draw_plot(p11, 0.3, 0.0, 0.3, 0.20) +
#            #cowplot::draw_label(labels[13], 0.51, 0.22, fontface = "bold") +
#            cowplot::draw_plot(leg, 0.6, 0.0, 0.4, 0.20)

    panel <- cowplot::ggdraw() +
        #cowplot::draw_label(labels[1], 0.01, 0.97, fontface = "bold") +
        cowplot::draw_plot(p1,  0.0, 0.6, 0.2, 0.4) +
        #cowplot::draw_label(labels[2], 0.51, 0.97, fontface = "bold") +
        cowplot::draw_plot(p2,  0.0, 0.2, 0.2, 0.4) +
        #cowplot::draw_label(labels[3], 0.01, 0.72, fontface = "bold") +
        cowplot::draw_plot(p3,  0.2, 0.6, 0.2, 0.4) +
        #cowplot::draw_label(labels[4], 0.51, 0.72, fontface = "bold") +
        cowplot::draw_plot(p4,  0.2, 0.2, 0.2, 0.4) +
        #cowplot::draw_label(labels[5], 0.01, 0.47, fontface = "bold") +
        cowplot::draw_plot(p5,  0.0, 0.0, 0.4, 0.2) +
        #cowplot::draw_label(labels[6], 0.51, 0.47, fontface = "bold") +
        cowplot::draw_plot(p6,  0.4, 0.6, 0.2, 0.4) +
        #cowplot::draw_label(labels[7], 0.01, 0.22, fontface = "bold") +
        cowplot::draw_plot(p7,  0.4, 0.2, 0.2, 0.4) +
        #cowplot::draw_label(labels[8], 0.51, 0.22, fontface = "bold") +
        cowplot::draw_plot(p8,  0.6, 0.6, 0.2, 0.4) +
        #cowplot::draw_label(labels[9], 0.51, 0.22, fontface = "bold") +
        cowplot::draw_plot(p9,  0.6, 0.2, 0.2, 0.4) +
        #cowplot::draw_label(labels[10], 0.51, 0.22, fontface = "bold") +
        cowplot::draw_plot(p12, 0.4, 0.0, 0.4, 0.2) +
        #cowplot::draw_label(labels[11], 0.51, 0.22, fontface = "bold") +
        cowplot::draw_plot(p10, 0.8, 0.6, 0.2, 0.4) +
        #cowplot::draw_label(labels[12], 0.51, 0.22, fontface = "bold") +
        cowplot::draw_plot(p11, 0.8, 0.2, 0.2, 0.4) +
        #cowplot::draw_label(labels[13], 0.51, 0.22, fontface = "bold") +
        cowplot::draw_plot(leg, 0.8, 0.0, 0.2, 0.2)

    return(panel)
}
