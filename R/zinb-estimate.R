#' Estimate ZINB-WaVE simulation parameters
#'
#' Estimate simulation parameters for the ZINB-WaVE simulation from a real
#' dataset.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param design.samples design matrix of sample-level covariates.
#' @param design.genes design matrix of gene-level covariates.
#' @param common.disp logical. Whether or not a single dispersion for all
#'        features is estimated.
#' @param iter.init number of iterations to use for initialization.
#' @param iter.opt number of iterations to use for optimization.
#' @param stop.opt stopping criterion for optimization.
#' @param params ZINBParams object to store estimated values in.
#' @param verbose logical. Whether to print progress messages.
#' @param BPPARAM A \code{\link{BiocParallelParam}} instance giving the parallel
#'        back-end to be used. Default is \code{\link{SerialParam}} which uses a
#'        single core.
#' @param ... additional arguments passes to \code{\link[zinbwave]{zinbFit}}.
#'
#' @details
#' The function is a wrapper around \code{\link[zinbwave]{zinbFit}} that takes
#' the fitted model and inserts it into a \code{\link{ZINBParams}} object. See
#' \code{\link{ZINBParams}} for more details on the parameters and
#' \code{\link[zinbwave]{zinbFit}} for details of the estimation procedure.
#'
#' @return ZINBParams object containing the estimated parameters.
#'
#' @examples
#' if (requireNamespace("zinbwave", quietly = TRUE)) {
#'     library(scater)
#'     set.seed(1)
#'     sce <- mockSCE(ncells = 20, ngenes = 100)
#'
#'     params <- zinbEstimate(sce)
#'     params
#' }
#'
#' @importFrom BiocParallel SerialParam
#' @export
zinbEstimate <- function(counts, design.samples = NULL, design.genes = NULL,
                         common.disp = TRUE, iter.init = 2, iter.opt = 25,
                         stop.opt = 1e-04, params = newZINBParams(),
                         verbose = TRUE, BPPARAM = SerialParam(), ...) {
    UseMethod("zinbEstimate")
}

#' @rdname zinbEstimate
#' @export
zinbEstimate.SingleCellExperiment <- function(counts, design.samples = NULL,
                                              design.genes = NULL,
                                              common.disp = TRUE,
                                              iter.init = 2, iter.opt = 25,
                                              stop.opt = 1e-04,
                                              params = newZINBParams(),
                                              verbose = TRUE,
                                              BPPARAM = SerialParam(), ...) {
    counts <- getCounts(counts)
    zinbEstimate(counts, design.samples, design.genes, common.disp,
                 iter.init, iter.opt, stop.opt, params, verbose, BPPARAM, ...)
}

#' @rdname zinbEstimate
#' @export
zinbEstimate.matrix <- function(counts, design.samples = NULL,
                                design.genes = NULL, common.disp = TRUE,
                                iter.init = 2, iter.opt = 25, stop.opt = 1e-04,
                                params = newZINBParams(), verbose = TRUE,
                                BPPARAM = SerialParam(), ...) {

    checkmate::assertClass(params, "ZINBParams")

    if (verbose) {message("Removing all zero genes...")}
    counts <- counts[rowSums(counts) > 0, ]

    args.list <- list(Y = counts,
                      commondispersion = common.disp,
                      verbose = verbose,
                      nb.repeat.initialize = iter.init,
                      maxiter.optimize = iter.opt,
                      stop.epsilon.optimize = stop.opt,
                      BPPARAM = BPPARAM)

    if (!is.null(design.samples)) {
        args.list$X <- design.samples
    }

    if (!is.null(design.genes)) {
        args.list$V <- design.genes
    }

    args.list <- c(args.list, list(...))

    if (verbose) {message("Fitting model...")}
    model <- do.call(zinbwave::zinbFit, args.list)

    params <- setParams(params, model = model)

    return(params)
}
