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
    fData.all$Dataset <- names(sces)[[1]]
    pData.all <- pData(sces[[1]])
    pData.all$Dataset <- names(sces)[[1]]

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            fData.all <- rbindMatched(fData.all, fData(sce))
            pData.all <- rbindMatched(pData.all, pData(sce))
        }
    }

    means <- ggplot(fData.all,
                    aes(x = Dataset, y = mean_log_cpm, colour = Dataset)) +
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        ylab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ggtitle("Distribution of mean expression") +
        theme_minimal() +
        theme(legend.position = "none")

    vars <- ggplot(fData.all,
                   aes(x = Dataset, y = var_cpm, colour = Dataset)) +
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        scale_y_log10(labels = scales::comma) +
        ylab("CPM Variance") +
        ggtitle("Distribution of variance") +
        theme_minimal()

    mean.var <- ggplot(fData.all,
                       aes(x = mean_log_cpm, y = var_log_cpm, colour = Dataset,
                           fill = Dataset)) +
        geom_point() +
        geom_smooth() +
        xlab(expression(paste("Mean ", log[2], "(CPM + 1)"))) +
        ylab(expression(paste("Variance ", log[2], "(CPM + 1)"))) +
        ggtitle("Mean-variance relationship") +
        theme_minimal()

    libs <- ggplot(pData.all,
                   aes(x = Dataset, y = total_counts, colour = Dataset)) +
        geom_boxplot() +
        scale_y_continuous(labels = scales::comma) +
        ylab("Total counts per cell") +
        ggtitle("Distribution of library sizes") +
        theme_minimal()

    z.gene <- ggplot(fData.all,
                     aes(x = Dataset, y = pct_dropout, colour = Dataset)) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        ylab("Percentage zeros per gene") +
        ggtitle("Distribution of zeros per gene") +
        theme_minimal()

    z.cell <- ggplot(pData.all,
                     aes(x = Dataset, y = pct_dropout, colour = Dataset)) +
        geom_boxplot() +
        scale_y_continuous(limits = c(0, 100)) +
        ylab("Percentage zeros per cell") +
        ggtitle("Distribution of zeros per cell") +
        theme_minimal()

    comparison <- list(FeatureData = fData.all,
                       PhenoData = pData.all,
                       Plots = list(Means = means,
                                    Variances = vars,
                                    MeanVar = mean.var,
                                    LibrarySizes = libs,
                                    ZerosGene = z.gene,
                                    ZerosCell = z.cell))

    return(comparison)
}

rbindMatched <- function(df1, df2) {
    common.names <- intersect(colnames(df1), colnames(df2))
    combined <- rbind(df1[, common.names], df2[, common.names])

    return(combined)
}

addFeatureStats <- function(sce, value = c("counts", "cpm", "tpm", "fpkm"),
                         log = FALSE, offset = 1, no.zeros = FALSE) {

    value <- match.arg(value)

    switch(value,
           counts = {
               values = scater::counts(sce)
           },
           cpm = {
               values = scater::cpm(sce)
           },
           tpm = {
               values = scater::tpm(sce)
           },
           fpkm = {
               values = scater::fpkm(sce)
           }
    )

    suffix = value

    if (no.zeros) {
        values[values == 0] <- NA
        suffix = paste0(suffix, "_no0")
    }

    if (log) {
        values = log2(values + offset)
        suffix = paste0("log_", suffix)
    }

    mean.str <- paste0("mean_", suffix)
    var.str  <- paste0("var_",  suffix)
    cv.str   <- paste0("cv_",   suffix)
    med.str  <- paste0("med_",  suffix)
    mad.str  <- paste0("mad_",  suffix)

    Biobase::fData(sce)[, mean.str] <- rowMeans(values, na.rm = TRUE)
    Biobase::fData(sce)[, var.str]  <- matrixStats::rowVars(values,
                                                            na.rm = TRUE)
    Biobase::fData(sce)[, cv.str]   <- sqrt(Biobase::fData(sce)[, var.str]) /
        Biobase::fData(sce)[, mean.str]
    Biobase::fData(sce)[, med.str]  <- matrixStats::rowMedians(values,
                                                               na.rm = TRUE)
    Biobase::fData(sce)[, mad.str]  <- matrixStats::rowMads(values,
                                                            na.rm = TRUE)
    return(sce)
}
