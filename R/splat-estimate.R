#' Estimate Splat simulation parameters
#'
#' Estimate simulation parameters for the Splat simulation from a real
#' dataset. See the individual estimation functions for more details on how this
#' is done.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param params SplatParams object to store estimated values in.
#'
#' @seealso
#' \code{\link{splatEstMean}},  \code{\link{splatEstLib}},
#' \code{\link{splatEstOutlier}}, \code{\link{splatEstBCV}},
#' \code{\link{splatEstDropout}}
#'
#' @return SplatParams object with estimated values.
#'
#' @examples
#' # Load example data
#' library(scater)
#' set.seed(1)
#' sce <- mockSCE()
#'
#' params <- splatEstimate(sce)
#' params
#' @export
splatEstimate <- function(counts, params = newSplatParams()) {
    UseMethod("splatEstimate")
}

#' @rdname splatEstimate
#' @export
splatEstimate.SingleCellExperiment <- function(counts,
                                               params = newSplatParams()) {
    counts <- getCounts(counts)
    splatEstimate(counts, params)
}

#' @rdname splatEstimate
#' @importFrom stats median
#' @export
splatEstimate.matrix <- function(counts, params = newSplatParams()) {

    checkmate::assertClass(params, "SplatParams")

    # Normalise for library size and remove all zero genes
    lib.sizes <- colSums(counts)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

    params <- splatEstMean(norm.counts, params)
    params <- splatEstLib(counts, params)
    params <- splatEstOutlier(norm.counts, params)
    params <- splatEstBCV(counts, params)
    params <- splatEstDropout(norm.counts, params)

    params <- setParams(params, nGenes = nrow(counts),
                        batchCells = ncol(counts))

    return(params)
}

#' Estimate Splat mean parameters
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate gene expression means.
#'
#' @param norm.counts library size normalised counts matrix.
#' @param params SplatParams object to store estimated values in.
#'
#' @details
#' Parameters for the gamma distribution are estimated by fitting the mean
#' normalised counts using \code{\link[fitdistrplus]{fitdist}}. The 'maximum
#' goodness-of-fit estimation' method is used to minimise the Cramer-von Mises
#' distance. This can fail in some situations, in which case the 'method of
#' moments estimation' method is used instead. Prior to fitting the means are
#' winsorized by setting the top and bottom 10 percent of values to the 10th
#' and 90th percentiles.
#'
#' @return SplatParams object containing the estimated parameters.
splatEstMean <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)
    means <- means[means != 0]

    means <- winsorize(means, q = 0.1)

    fit <- fitdistrplus::fitdist(means, "gamma", method = "mge",
                                 gof = "CvM")
    if (fit$convergence > 0) {
        warning("Fitting means using the Goodness of Fit method failed, ",
                "using the Method of Moments instead")
        fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")
    }

    params <- setParams(params, mean.shape = unname(fit$estimate["shape"]),
                        mean.rate = unname(fit$estimate["rate"]))

    return(params)
}

#' Estimate Splat library size parameters
#'
#' The Shapiro-Wilks test is used to determine if the library sizes are
#' normally distributed. If so a normal distribution is fitted to the library
#' sizes, if not (most cases) a log-normal distribution is fitted and the
#' estimated parameters are added to the params object. See
#' \code{\link[fitdistrplus]{fitdist}} for details on the fitting.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param params splatParams object to store estimated values in.
#'
#' @return SplatParams object with estimated values.
#'
#' @importFrom stats shapiro.test
splatEstLib <- function(counts, params) {

    lib.sizes <- colSums(counts)

    if (length(lib.sizes) > 5000) {
        message("NOTE: More than 5000 cells provided. ",
                "5000 sampled library sizes will be used to test normality.")
        lib.sizes.sampled <- sample(lib.sizes, 5000, replace = FALSE)
    } else {
        lib.sizes.sampled <- lib.sizes
    }

    norm.test <- shapiro.test(lib.sizes.sampled)
    lib.norm <- norm.test$p.value > 0.2

    if (lib.norm) {
        fit <- fitdistrplus::fitdist(lib.sizes, "norm")
        lib.loc <- unname(fit$estimate["mean"])
        lib.scale <- unname(fit$estimate["sd"])
        message("NOTE: Library sizes have been found to be normally ",
                "distributed instead of log-normal. You may want to check ",
                "this is correct.")
    } else {
        fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")
        lib.loc <- unname(fit$estimate["meanlog"])
        lib.scale <- unname(fit$estimate["sdlog"])
    }

    params <- setParams(params, lib.loc = lib.loc, lib.scale = lib.scale,
                        lib.norm = lib.norm)

    return(params)
}

#' Estimate Splat expression outlier parameters
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
#' than two MADs above the median log2 mean expression it is designated as an
#' outlier. The proportion of outlier genes is used to estimate the outlier
#' probability. Factors for each outlier gene are calculated by dividing mean
#' expression by the median mean expression. A log-normal distribution is then
#' fitted to these factors in order to estimate the outlier factor location and
#' scale parameters using \code{\link[fitdistrplus]{fitdist}}.
#'
#' @return SplatParams object with estimated values.
splatEstOutlier <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)
    lmeans <- log(means)

    med <- median(lmeans)
    mad <- mad(lmeans)

    bound <- med + 2 * mad

    outs <- which(lmeans > bound)

    prob <- length(outs) / nrow(norm.counts)

    params <- setParams(params, out.prob = prob)

    if (length(outs) > 1) {
        facs <- means[outs] / median(means)
        fit <- fitdistrplus::fitdist(facs, "lnorm")

        params <- setParams(params,
                            out.facLoc = unname(fit$estimate["meanlog"]),
                            out.facScale = unname(fit$estimate["sdlog"]))
    }

    return(params)
}

#' Estimate Splat Biological Coefficient of Variation parameters
#'
#' Parameters are estimated using the \code{\link[edgeR]{estimateDisp}} function
#' in the \code{edgeR} package.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param params SplatParams object to store estimated values in.
#'
#' @details
#' The \code{\link[edgeR]{estimateDisp}} function is used to estimate the common
#' dispersion and prior degrees of freedom. See
#' \code{\link[edgeR]{estimateDisp}} for details. When estimating parameters on
#' simulated data we found a broadly linear relationship between the true
#' underlying common dispersion and the \code{edgR} estimate, therefore we
#' apply a small correction, \code{disp = 0.1 + 0.25 * edgeR.disp}.
#'
#' @return SplatParams object with estimated values.
splatEstBCV <- function(counts, params) {

    # Add dummy design matrix to avoid print statement
    design <- matrix(1, ncol(counts), 1)
    disps <- edgeR::estimateDisp(counts, design = design)

    params <- setParams(params,
                        bcv.common = 0.1 + 0.25 * disps$common.dispersion,
                        bcv.df = disps$prior.df)

    return(params)
}

#' Estimate Splat dropout parameters
#'
#' Estimate the midpoint and shape parameters for the logistic function used
#' when simulating dropout.
#'
# #' Also estimates whether dropout is likely to be
# #' present in the dataset.
#'
#' @param norm.counts library size normalised counts matrix.
#' @param params SplatParams object to store estimated values in.
#'
#' @details
#' Logistic function parameters are estimated by fitting a logistic function
#' to the relationship between log2 mean gene expression and the proportion of
#' zeros in each gene. See \code{\link[stats]{nls}} for details of fitting.
#' Note this is done on the experiment level, more granular (eg. group or cell)
#' level dropout is not estimated.
#'
# #' The
# #' presence of dropout is determined by comparing the observed number of zeros
# #' in each gene to the expected number of zeros from a negative binomial
# #' distribution with the gene mean and a dispersion of 0.1. If the maximum
# #' difference between the observed number of zeros and the expected number is
# #' greater than 10 percent of the number of cells
# #' (\code{max(obs.zeros - exp.zeros) > 0.1 * ncol(norm.counts)}) then dropout
# #' is considered to be present in the dataset. This is a somewhat crude
# #' measure but should give a reasonable indication. A more accurate approach
# #' is to look at a plot of log2 mean expression vs the difference between
# #' observed and expected number of zeros across all genes.
#'
#' @return SplatParams object with estimated values.
#'
#' @importFrom stats dnbinom nls
splatEstDropout <- function(norm.counts, params) {

    means <- rowMeans(norm.counts)

    x <- log(means)

    obs.zeros <- rowSums(norm.counts == 0)

    y <- obs.zeros / ncol(norm.counts)

    df <- data.frame(x, y)

    fit <- tryCatch({
        nls(y ~ logistic(x, x0 = x0, k = k), data = df,
            start = list(x0 = 0, k = -1))
    },
    error = function(err) {
        warning("Fitting dropout using the Gauss-Newton method failed, ",
                "using the Golub-Pereyra algorithm instead")
        nls(y ~ logistic(x, x0 = x0, k = k), data = df,
            start = list(x0 = 0, k = -1), algorithm = "plinear")
    })

    #exp.zeros <- dnbinom(0, mu = means, size = 1 / 0.1) * ncol(norm.counts)
    #present <- max(obs.zeros - exp.zeros) > 0.1 * ncol(norm.counts)

    mid <- summary(fit)$coefficients["x0", "Estimate"]
    shape <- summary(fit)$coefficients["k", "Estimate"]

    params <- setParams(params, dropout.mid = mid, dropout.shape = shape)

    return(params)
}
