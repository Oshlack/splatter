#' Estimate population/eQTL simulation parameters
#'
#' Estimate simulation parameters for the eQTL population simulation from
#' real data. See the individual estimation functions for more details on
#' how this is done.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param means Matrix of real gene means across a population, where
#'        each row is a gene and each column is an individual in the population.
#' @param eqtl data.frame with all or top eQTL pairs from a real eQTL analysis.
#'        Must include columns: 'gene_id', 'pval_nominal', and 'slope'.
#' @param params SplatPopParams object containing parameters for the
#'        simulation of the mean expression levels for the population.
#'        See \code{\link{SplatPopParams}} for details.
#'
#' @seealso
#' \code{\link{splatPopEstimateEffectSize}},
#' \code{\link{splatPopEstimateMeanCV}}
#'
#' @return SplatPopParams object containing the estimated parameters.
#'
#' @examples
#'
#' if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
#'     requireNamespace("preprocessCore", quietly = TRUE)) {
#'     # Load example data
#'     library(scater)
#'
#'     sce <- mockSCE()
#'     params <- splatPopEstimate(sce)
#' }
#'
#' @export
splatPopEstimate <- function(counts = NULL, means = NULL, eqtl = NULL,
                             params = newSplatPopParams()) {

    checkmate::assertClass(params, "SplatPopParams")

    # Estimate single-cell parameters using base splatEstimate function
    if (!is.null(counts)) {
        params <- splatEstimate(counts, params)
    }

    # Get parameters for eQTL Effect Size distribution
    if (!is.null(eqtl)) {
        params <- splatPopEstimateEffectSize(params, eqtl)
    }

    # Get parameters for population wide gene mean and variance distributions
    if (!is.null(means)) {
        params <- splatPopEstimateMeanCV(params, means)
    }

    return(params)

}

#' Estimate eQTL Effect Size parameters
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate eQTL (eSNP-eGene) effect sizes.
#'
#' @param eqtl data.frame with all or top eQTL pairs from a real eQTL analysis.
#'        Must include columns: gene_id, pval_nominal, and slope.
#' @param params SplatPopParams object containing parameters for the
#'        simulation of the mean expression levels for the population.
#'        See \code{\link{SplatPopParams}} for details.
#'
#' @details
#' Parameters for the gamma distribution are estimated by fitting the top eSNP-
#' eGene pair effect sizes using \code{\link[fitdistrplus]{fitdist}}. The
#' maximum goodness-of-fit estimation method is used to minimise the
#' Cramer-von Mises distance. This can fail in some situations, in which case
#' the method of moments estimation method is used instead.
#'
#' @return params object with estimated values.
#'
splatPopEstimateEffectSize <- function(params, eqtl) {

    # Test input eSNP-eGene pairs
    if (!("gene_id" %in% names(eqtl) &
          "pval_nominal" %in% names(eqtl) &
          "slope" %in% names(eqtl))) {
        stop("Incorrect format for eqtl data.")}

    # Select top eSNP for each gene (i.e. lowest p.value)
    eqtl.top <- eqtl[order(eqtl$gene_id, eqtl$pval_nominal), ]
    eqtl.top <- eqtl.top[!duplicated(eqtl.top$gene_id), ]

    # Fit absolute value of effect sizes to gamma distribution
    e.sizes <- abs(eqtl.top$slope)
    fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mge", gof = "CvM")
    if (fit$convergence > 0) {
        warning("Fitting effect sizes using the Goodness of Fit method failed,",
                " using the Method of Moments instead")
        fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mme")
    }

    params <- setParams(params, eqtl.ES.shape = unname(fit$estimate["shape"]),
                        eqtl.ES.rate = unname(fit$estimate["rate"]))

    return(params)
}

#' Estimate gene mean and gene mean variance parameters
#'
#' @param params SplatPopParams object containing parameters for the
#'        simulation of the mean expression levels for the population.
#'        See \code{\link{SplatPopParams}} for details.
#' @param means data.frame of real gene means across a population, where
#'        each row is a gene and each column is an individual in the population.
#'
#' @details
#' Parameters for the mean gamma distribution are estimated by fitting the mean
#' (across the population) expression of genes that meet the criteria (<50% of
#' samples have exp <0.1) and parameters for the cv gamma distribution are
#' estimated for each bin of mean expression using the cv of expression across
#' the population for genes in that bin. Both are fit using
#' \code{\link[fitdistrplus]{fitdist}}. The "Nelder-Mead" method is used to fit
#' the mean gamma distribution and the maximum goodness-of-fit estimation
#' method is used to minimise the Cramer-von Mises distance for the CV
#' distribution.
#'
#' @return params object with estimated values.
#' @importFrom stats quantile
#' @importFrom grDevices boxplot.stats
#' @importFrom matrixStats rowMedians
#'
splatPopEstimateMeanCV <- function(params, means) {

    # Test input gene means
    if ((anyNA(means) | !(validObject(rowSums(means))))) {
        stop("Incorrect format or NAs present in gene.means. See example data.")
    }

    # Remove genes with low variance/low means
    abv.thr <- data.frame(perc = (rowSums(means >= 0.1)/ncol(means)))

    means.use <- means[abv.thr$perc > 0.5, ]

    # Calculate mean expression parameters
    row.means <- rowMedians(means.use)
    names(row.means) <- row.names(means.use)
    mfit <- fitdistrplus::fitdist(row.means, "gamma",
                                  optim.method = "Nelder-Mead")

    # Calculate CV parameters for genes based on 10 mean expression bins
    nbins <- getParam(params, "pop.cv.bins")
    bins <- split(row.means, cut(row.means, quantile(row.means,(0:nbins)/nbins),
                             include.lowest = TRUE))
    cvparams <- data.frame(start = character(), end = character(),
                           shape = character(), rate = character(),
                           stringsAsFactors = FALSE)

    for(b in names(bins)){
        re.brack.paren <- "\\[|\\]|\\)|\\("
        min.max <- strsplit(gsub(re.brack.paren, "", unlist(b)), split = ",")

        b.gene.means <- means.use[row.means > as.numeric(min.max[[1]][1]) &
                                      row.means < as.numeric(min.max[[1]][2]), ]

        cv <- apply(b.gene.means, 1, co.var)
        cv[is.na(cv)] <- 0
        cv <- cv[!cv %in% boxplot.stats(cv)$out]
        cvfit <- fitdistrplus::fitdist(cv, "gamma", method = "mge", gof = "CvM")
        cvparams <- rbind(cvparams,
                          list(start = as.numeric(as.numeric(min.max[[1]][1])),
                               end = as.numeric(as.numeric(min.max[[1]][2])),
                               shape = cvfit$estimate["shape"],
                               rate = cvfit$estimate["rate"]),
                          stringsAsFactors = FALSE)
    }

    cvparams[1, "start"] <- 0
    cvparams[nrow(cvparams), "end"] <- 1e100
    params <- setParams(params,
                        pop.mean.shape = unname(mfit$estimate["shape"]),
                        pop.mean.rate = unname(mfit$estimate["rate"]),
                        pop.cv.param = cvparams)

    return(params)
}
