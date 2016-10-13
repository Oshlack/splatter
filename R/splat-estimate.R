#' Estimate Splatter simulation parameters
#'
#' Estimate simulation parameters for the Splatter simulation from a real
#' dataset. See the individual estimation functions for more details on how this
#' is done.
#'
#' @param counts either a counts matrix or an SCESet object containing count
#'        data to estimate parameters from.
#' @param params SplatParams object to store estimated values in.
#'
#' @seealso
#' \code{\link{estSplatMean}},  \code{\link{estSplatLib}},
#' \code{\link{estSplatOutlier}}, \code{\link{estSplatBCV}},
#' \code{\link{estSplatDropout}}
#'
#' @return SplatParams object containing the estimated parameters.
#'
#' @examples
#' data("sc_example_counts")
#' params <- estimateSplatParams(sc_example_counts)
#' params
#' @export
estimateSplatParams <- function(counts, params = newSplatParams()) {
    UseMethod("estimateSplatParams")
}

#' @rdname estimateSplatParams
#' @export
estimateSplatParams.SCESet <- function(counts, params = newSplatParams()) {
    counts <- scater::counts(counts)
    estimateSplatParams(counts, params)
}

#' @rdname estimateSplatParams
#' @importFrom stats median
#' @export
estimateSplatParams.matrix <- function(counts, params = newSplatParams()) {

    checkmate::assertClass(params, "SplatParams")

    # Normalise for library size and remove all zero genes
    lib.sizes <- colSums(counts)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

    params <- estSplatMean(norm.counts, params)
    params <- estSplatLib(counts, params)
    params <- estSplatOutlier(norm.counts, params)
    params <- estSplatBCV(counts, params)
    params <- estSplatDropout(norm.counts, params)

    params <- setParams(params, nGenes = nrow(counts),
                        groupCells = ncol(counts))

    return(params)
}

#' Estimate Splatter mean parameters
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate gene expression means using the 'moment matching estimation' method
#' of \code{\link[fitdistrplus]{fitdist}}.
#'
#' @param norm.counts library size normalised counts matrix.
#' @param params SplatParams object to store estimated values in.
#'
#' @return SplatParams object with estimated values.
estSplatMean <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)
    means <- means[means != 0]

    fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")

    params <- setParams(params, mean.shape = unname(fit$estimate["shape"]),
                        mean.rate = unname(fit$estimate["rate"]))

    return(params)
}

#' Estimate Splatter library size parameters
#'
#' A log-normal distribution is fitted to the library sizes and the estimated
#' parameters are added to the params object. See
#' \code{\link[fitdistrplus]{fitdist}} for details on the fitting.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param params splatParams object to store estimated values in.
#'
#' @return splatParams object with estimated values.
estSplatLib <- function(counts, params) {

    lib.sizes <- colSums(counts)
    fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")

    params <- setParams(params, lib.loc = unname(fit$estimate["meanlog"]),
                        lib.scale = unname(fit$estimate["sdlog"]))

    return(params)
}

#' Estimate Splatter expression outlier parameters
#'
#' Parameters are estimated by comparing means of individual genes to the
#' median mean expression level.
#'
#' @param norm.counts library size normalised counts matrix.
#' @param params SplatParams object to store estimated values in.
#'
#' @details
#' Expression outlier genes are detected using the Median Absolute Deviation
#' (MAD) from median method. If the log2 mean expression of a gene is greater
#' than two MADs from the median log2 mean expression it is designated as a
#' outlier. The proportion of outlier genes is used to estimate the outlier
#' probability. The low outlier probability is estimated as the proportion of
#' outlier genes that have a log2 mean less than the median log2 mean. Factors
#' for each outlier gene are calculated by dividing mean expression by the
#' median mean expression. A log-normal distribution is then fitted to these
#' factors in order to estimate the outlier factor location and scale
#' parameters. See \code{\link[fitdistrplus]{fitdist}} for details on the
#' fitting.
#'
#' @return SplatParams object with estimated values.
estSplatOutlier <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)
    lmeans <- log(means)

    med <- median(lmeans)
    mad <- mad(lmeans)

    lo.bound <- med - 2 * mad
    hi.bound <- med + 2 * mad

    lo.outs <- which(lmeans < lo.bound)
    hi.outs <- which(lmeans > hi.bound)

    prob <- (length(lo.outs) + length(hi.outs)) / nrow(norm.counts)
    lo.prob <- length(lo.outs) / (length(lo.outs) + length(hi.outs))

    facs <- means[c(lo.outs, hi.outs)] / median(means)
    fit <- fitdistrplus::fitdist(facs, "lnorm")

    params <- setParams(params, out.prob = prob, out.loProb = lo.prob,
                        out.facLoc = unname(fit$estimate["meanlog"]),
                        out.facScale = unname(fit$estimate["sdlog"]))

    return(params)
}

#' Estimate Splatter Biological Coefficient of Variation parameters
#'
#' Parameters are estimated using the \code{estimateDisp} function in the
#' \code{edgeR} package. Specifically the common dispersion and prior degrees
#' of freedom. See \code{\link{estimateDisp}} for details.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param params SplatParams object to store estimated values in.
#'
#' @return SplatParams object with estimated values.
estSplatBCV <- function(counts, params) {

    # Add dummy design matrix to avoid print statement
    design <- matrix(1, ncol(counts), 1)
    disps <- edgeR::estimateDisp(counts, design = design)

    params <- setParams(params, bcv.common = disps$common.dispersion,
                        bcv.df = disps$prior.df)

    return(params)
}

#' Estimate Splatter dropout parameters
#'
#' Estimate the midpoint and shape parameters for the logistic function used
#' when simulating dropout. Also estimates whether dropout is likely to be
#' present in the dataset.
#'
#' @param norm.counts library size normalised counts matrix.
#' @param params SplatParams object to store estimated values in.
#'
#' @details
#' Logistic function parameters are estimated by fitting a logistic function
#' to the relationship between log2 mean gene expression and the proportion of
#' zeros in each gene. See \code{\link[stats]{nls}} for details of fitting. The
#' presence of dropout is determined by comparing the observed number of zeros
#' in each gene to the expected number of zeros from a negative binomial
#' distribution with the gene mean and a dispersion of 0.1. If the maximum
#' difference between the observed number of zeros and the expected number is
#' greater than 10 percent of the number of cells
#' (\code{max(obs.zeros - exp.zeros) > 0.1 * ncol(norm.counts)}) then dropout is
#' considered to be present in the dataset. This is a somewhat crude measure
#' but should give a reasonable indication. A more accurate approach is to look
#' at a plot of log2 mean expression vs the difference between observed and
#' expected number of zeros across all genes.
#'
#' @return SplatParams object with estimated values.
estSplatDropout <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)

    x <- log(means)

    obs.zeros <- rowSums(norm.counts == 0)

    y <- obs.zeros / ncol(norm.counts)

    df <- data.frame(x, y)

    fit <- nls(y ~ logistic(x, x0 = x0, k = k), data = df,
               start = list(x0 = 0, k = -1))

    exp.zeros <- dnbinom(0, mu = means, size = 1 / 0.1) * ncol(norm.counts)

    present <- max(obs.zeros - exp.zeros) > 0.1 * ncol(norm.counts)

    mid <- summary(fit)$coefficients["x0", "Estimate"]
    shape <- summary(fit)$coefficients["k", "Estimate"]

    params <- setParams(params, dropout.present = present, dropout.mid = mid,
                        dropout.shape = shape)

    return(params)
}