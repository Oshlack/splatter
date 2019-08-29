#' Estimate Splotch simulation parameters
#'
#' Estimate simulation parameters for the Splotch simulation from a real
#' dataset. See the individual estimation functions for more details on how this
#' is done.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param params SplotchParams object to store estimated values in.
#' @param verbose logical. Whether to print progress messages.
#'
#' @seealso
#' \code{\link{splotchEstMean}}, \code{\link{splotchEstLib}}
#'
#' @return SplotchParams object containing the estimated parameters.
#'
#' @examples
#' # Load example data
#' library(scater)
#' set.seed(1)
#' sce <- mockSCE()
#'
#' params <- splotchEstimate(sce)
#' params
#' @export
splotchEstimate <- function(counts, params = newSplotchParams(),
                            verbose = TRUE) {
    UseMethod("splotchEstimate")
}

#' @rdname splotchEstimate
#' @export
splotchEstimate.SingleCellExperiment <- function(counts,
                                                 params = newSplotchParams(),
                                                 verbose = TRUE) {
    counts <- BiocGenerics::counts(counts)
    splotchEstimate(counts, params, verbose)
}

#' @rdname splotchEstimate
#' @importFrom stats median
#' @export
splotchEstimate.matrix <- function(counts, params = newSplotchParams(),
                                   verbose = TRUE) {

    checkmate::assertClass(params, "SplotchParams")
    checkmate::assertFlag(verbose)

    # Normalise for library size and remove all zero genes
    lib.sizes <- colSums(counts)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

    params <- splotchEstMean(norm.counts, params, verbose)
    params <- splotchEstBCV(counts, params, verbose)
    params <- splotchEstLib(counts, params, verbose)

    params <- setParams(params, nGenes = nrow(counts), nCells = ncol(counts))

    return(params)
}

splotchEstMean <- function(norm.counts, params, verbose) {

    if (verbose) {message("Estimating mean parameters...")}

    means <- rowMeans(norm.counts)
    means <- means[means != 0]
    non.zero <- rowSums(norm.counts > 0)

    fit <- selectFit(means, "gamma", non.zero, verbose)

    params <- setParams(params, mean.shape = unname(fit$estimate["shape"]),
                        mean.rate = unname(fit$estimate["rate"]))

    if (verbose) {message("Estimating expression outlier parameters...")}
    lmeans <- log(means)
    med <- median(lmeans)
    mad <- mad(lmeans)

    bound <- med + 2 * mad

    outs <- which(lmeans > bound)

    prob <- length(outs) / nrow(norm.counts)

    params <- setParam(params, "mean.outProb", prob)

    if (length(outs) > 1) {
        facs <- means[outs] / median(means)
        #fit <- selectFit(facs, "lnorm", verbose = verbose)
        fit <- fitdistrplus::fitdist(facs, "lnorm")

        params <- setParams(params,
                            mean.outLoc = unname(fit$estimate["meanlog"]),
                            mean.outScale = unname(fit$estimate["sdlog"]))
    }

    params <- setParams(params, mean.dens = density(lmeans))

    return(params)
}

splotchEstBCV <- function(counts, params, verbose) {

    if (verbose) {message("Estimating BCV parameters...")}

    # Add dummy design matrix to avoid print statement
    design <- matrix(1, ncol(counts), 1)
    disps <- edgeR::estimateDisp(counts, design = design)
    raw <- disps$common.dispersion

    mean.rate <- getParam(params, "mean.rate")
    mean.shape <- getParam(params, "mean.shape")

    # Caculate factors for correction based on mean parameters
    # Coefficents come from fitting
    #   RawEst ~
    #       (A1 * mean.rate + A2 * mean.shape + A3 * mean.rate *
    #           E_mean.shape + A4) *
    #       ((B1 * mean.rate + B2 * mean.shape + B3 * mean.rate *
    #           mean.shape + B4) ^ SimBCVCommon) +
    #       (C1 * mean.rate + C2 * mean.shape + C3 * mean.rate *
    #           mean.shape + C4)
    # Using minpack.lm::nlsLM
    A <- -0.6 * mean.rate - 2.9 * mean.shape +
        0.4 * mean.rate * mean.shape + 9.5
    B <- 0.15 * mean.rate + 0.25 * mean.shape -
        0.1 * mean.rate * mean.shape + 1.2
    C <- 0.9 * mean.rate + 4.5 * mean.shape -
        0.6 * mean.rate * mean.shape - 10.6
    Y <- (raw - C) / A

    message("Raw: ", raw, " A: ", A, " B: ", B, " C: ", C, " Y: ", Y)

    # Check if Y <= 0 to avoid problems when taking log
    if (Y <= 0) {
        Y <- 0.0001
    }

    corrected <- log(Y, base = B)

    # Dispersion cannot be negative so apply a simpler correction to those.
    # Coefficients come from fitting
    #   SimBCVCommon ~ EstRaw
    # Using lm (negative values only)
    if (corrected < 0) {
        warning("Exponential corrected BCV is negative.",
                "Using linear correction.")
        corrected <- -0.1 + 0.1 * raw
    }

    if (corrected < 0) {
        warning("Linear corrected BCV is negative.",
                "Using existing bcv.common.")
        corrected <- getParam(params, "bcv.common")
    }

    params <- setParams(params, bcv.common = corrected)

    return(params)
}

#' @importFrom stats density
splotchEstLib <- function(counts, params, verbose) {

    if (verbose) {message("Estimating library size parameters...")}

    lib.sizes <- colSums(counts)

    fit <- selectFit(lib.sizes, "lnorm", verbose = verbose)

    lib.loc <- unname(fit$estimate["meanlog"])
    lib.scale <- unname(fit$estimate["sdlog"])

    dens <- density(lib.sizes)

    params <- setParams(params, lib.loc = lib.loc, lib.scale = lib.scale,
                        lib.dens = dens)

    return(params)
}

selectFit <- function(data, distr, weights = NULL, verbose = TRUE) {

    checkmate::assertNumeric(data, finite = TRUE, any.missing = FALSE)
    checkmate::assertString(distr)
    checkmate::assertNumeric(weights, finite = TRUE, any.missing = TRUE,
                             len = length(data), null.ok = TRUE)
    checkmate::assertFlag(verbose)

    # Sink output that sometimes happens when fitting
    sink(tempfile())
    on.exit(sink())

    fits <- list()

    try(
        fits$`MLE` <- fitdistrplus::fitdist(data, distr, method = "mle"),
        silent = TRUE
    )

    try(
        fits$`MME` <- fitdistrplus::fitdist(data, distr, method = "mme"),
        silent = TRUE
    )

    try(
        fits$`QME` <- fitdistrplus::fitdist(data, distr, method = "qme",
                                            probs = c(1/3, 2/3)),
        silent = TRUE
    )

    try(
        fits$`MGE (CvM)` <- fitdistrplus::fitdist(data, distr, method = "mge",
                                                  gof = "CvM"),
        silent = TRUE
    )

    try(
        fits$`MGE (KS)` <- fitdistrplus::fitdist(data, distr, method = "mge",
                                                 gof = "KS"),
        silent = TRUE
    )

    try(
        fits$`MGE (AD)` <- fitdistrplus::fitdist(data, distr, method = "mge",
                                                 gof = "AD"),
        silent = TRUE
    )

    if (!is.null(weights)) {
        try(suppressWarnings(
            fits$`Weighted MLE` <- fitdistrplus::fitdist(data, distr,
                                                         method = "mle",
                                                         weights = weights)
            ), silent = TRUE
        )

        try(suppressWarnings(
            fits$`Weighted MME` <- fitdistrplus::fitdist(data, distr,
                                                         method = "mme",
                                                         weights = weights)
            ), silent = TRUE
        )

        try(suppressWarnings(
            fits$`Weighted QME` <- fitdistrplus::fitdist(data, distr,
                                                         method = "qme",
                                                         probs = c(1/3, 2/3),
                                                         weights = weights)
            ), silent = TRUE
        )
    }

    scores <- fitdistrplus::gofstat(fits)$cvm

    # Flatten in case scores is a list
    scores.flat <- unlist(scores)
    selected <- which(scores.flat == min(scores.flat, na.rm = TRUE))

    if (verbose) {
        # Work around to get name in case scores is a list
        name <- names(fits)[names(scores) == names(scores.flat)[selected]]
        message("Selected ", name, " fit")
    }

    return(fits[[selected]])
}
