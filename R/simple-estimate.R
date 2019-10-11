#' Estimate simple simulation parameters
#'
#' Estimate simulation parameters for the simple simulation from a real dataset.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param params SimpleParams object to store estimated values in.
#'
#' @details
#' The \code{nGenes} and \code{nCells} parameters are taken from the size of the
#' input data. The mean parameters are estimated by fitting a gamma distribution
#' to the library size normalised mean expression level using
#' \code{\link[fitdistrplus]{fitdist}}. See \code{\link{SimpleParams}} for more
#' details on the parameters.
#'
#' @return SimpleParams object containing the estimated parameters.
#'
#' @examples
#' # Load example data
#' library(scater)
#' set.seed(1)
#' sce <- mockSCE()
#'
#' params <- simpleEstimate(sce)
#' params
#' @export
simpleEstimate <- function(counts, params = newSimpleParams()) {
    UseMethod("simpleEstimate")
}

#' @rdname simpleEstimate
#' @export
simpleEstimate.SingleCellExperiment <- function(counts,
                                                params = newSimpleParams()) {
    counts <- getCounts(counts)
    simpleEstimate(counts, params)
}

#' @rdname simpleEstimate
#' @importFrom stats median
#' @export
simpleEstimate.matrix <- function(counts, params = newSimpleParams()) {

    checkmate::assertClass(params, "SimpleParams")

    # Normalise for library size and remove all zero genes
    lib.sizes <- colSums(counts)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(counts) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

    means <- rowMeans(norm.counts)

    means.fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")

    params <- setParams(params, nGenes = nrow(counts), nCells = ncol(counts),
                        mean.shape = unname(means.fit$estimate["shape"]),
                        mean.rate = unname(means.fit$estimate["rate"]))

    return(params)
}
