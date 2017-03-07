#' Compare SCESet objects
#'
#' Combine the data from several SCESet objects and produce some basic plots
#' comparing.
#'
#' @param sces named list of SCESet objects to combine and compare.
#'
#' @details
#' The return list has three items:
#'
#' \describe{
#'     \item{\code{FeatureData}}{Combined feature data from the provided
#'     SCESets.}
#'     \item{\code{PhenoData}}{Combined pheno data from the provided SCESets.}
#'     \item{\code{Plots}}{Comparison plots
#'         \describe{
#'             \item{\code{Means}}{Violin plot of mean distribution.}
#'             \item{\code{Variances}}{Violin plot of variance distribution.}
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
#' geom_violin scale_y_continuous scale_y_log10 xlab ylab ggtitle theme_minimal
#' @importFrom scater cpm<-
#' @export
compareSCESets <- function(sces) {

    checkmate::assertList(sces, types = "SCESet", any.missing = FALSE,
                          min.len = 1, names = "unique")

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
                    aes_string(x = "Dataset", y = "mean_log_cpm",
                               colour = "Dataset")) +
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        ylab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ggtitle("Distribution of mean expression") +
        theme_minimal()

    vars <- ggplot(fData.all,
                   aes_string(x = "Dataset", y = "var_cpm",
                              colour = "Dataset")) +
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        scale_y_log10(labels = scales::comma) +
        ylab("CPM Variance") +
        ggtitle("Distribution of variance") +
        theme_minimal()

    mean.var <- ggplot(fData.all,
                       aes_string(x = "mean_log_cpm", y = "var_log_cpm",
                                  colour = "Dataset", fill = "Dataset")) +
        geom_point(alpha = 0.2) +
        geom_smooth() +
        xlab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ylab(expression(paste("Variance ", log[2], "(CPM + 1)"))) +
        ggtitle("Mean-variance relationship") +
        theme_minimal()

    libs <- ggplot(pData.all,
                   aes_string(x = "Dataset", y = "total_counts",
                              colour = "Dataset")) +
        geom_boxplot() +
        scale_y_continuous(labels = scales::comma) +
        ylab("Total counts per cell") +
        ggtitle("Distribution of library sizes") +
        theme_minimal()

    z.gene <- ggplot(fData.all,
                     aes_string(x = "Dataset", y = "pct_dropout",
                                colour = "Dataset")) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        ylab("Percentage zeros per gene") +
        ggtitle("Distribution of zeros per gene") +
        theme_minimal()

    z.cell <- ggplot(pData.all,
                     aes_string(x = "Dataset", y = "pct_dropout",
                                colour = "Dataset")) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        ylab("Percentage zeros per cell") +
        ggtitle("Distribution of zeros per cell") +
        theme_minimal()

    mean.zeros <- ggplot(fData.all,
                         aes_string(x = "mean_log_cpm", y = "pct_dropout",
                                    colour = "Dataset", fill = "Dataset")) +
        geom_point(alpha = 0.2) +
        geom_smooth() +
        xlab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ylab(expression(paste("Percentage zeros"))) +
        ggtitle("Mean-dropout relationship") +
        theme_minimal()

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
