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
#' data("sc_example_counts")
#'
#' params <- splotchEstimate(sc_example_counts)
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

    bound <- med + 1 * mad

    outs <- which(lmeans > bound)

    prob <- length(outs) / nrow(norm.counts)

    params <- setParam(params, "mean.outProb", prob)

    if (length(outs) > 1) {
        facs <- means[outs] / median(means)
        fit <- selectFit(facs, "lnorm", verbose = verbose)

        params <- setParams(params,
                            mean.outLoc = unname(fit$estimate["meanlog"]),
                            mean.outScale = unname(fit$estimate["sdlog"]))
    }

    return(params)
}

splotchEstLib <- function(counts, params, verbose) {

    if (verbose) {message("Estimating library size parameters...")}

    lib.sizes <- colSums(counts)

    fit <- selectFit(lib.sizes, "lnorm", verbose = verbose)

    lib.loc <- unname(fit$estimate["meanlog"])
    lib.scale <- unname(fit$estimate["sdlog"])

    params <- setParams(params, lib.loc = lib.loc, lib.scale = lib.scale)

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
