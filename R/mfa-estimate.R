#' Estimate mfa simulation parameters
#'
#' Estimate simulation parameters for the mfa simulation from a real dataset.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param params MFAParams object to store estimated values in.
#'
#' @details
#' The \code{nGenes} and \code{nCells} parameters are taken from the size of the
#' input data. The dropout lambda parameter is estimate using
#' \code{\link[mfa]{empirical_lambda}}. See \code{\link{MFAParams}} for more
#' details on the parameters.
#'
#' @return MFAParams object containing the estimated parameters.
#'
#' @examples
#' # Load example data
#' if (requireNamespace("mfa", quietly = TRUE)) {
#'     library(mfa)
#'     synth <- create_synthetic(C = 20, G = 5, zero_negative = TRUE,
#'                               model_dropout = TRUE)
#'
#'     params <- mfaEstimate(synth$X)
#'     params
#' }
#' @export
mfaEstimate <- function(counts, params = newMFAParams()) {
    UseMethod("mfaEstimate")
}

#' @rdname mfaEstimate
#' @export
mfaEstimate.SingleCellExperiment <- function(counts,
                                             params = newMFAParams()) {
    counts <- getCounts(counts)
    mfaEstimate(counts, params)
}

#' @rdname mfaEstimate
#' @export
mfaEstimate.matrix <- function(counts, params = newMFAParams()) {

    checkmate::assertClass(params, "MFAParams")

    dropout.lambda <- mfa::empirical_lambda(t(counts))

    params <- setParams(params, nGenes = nrow(counts), nCells = ncol(counts),
                        dropout.lambda = dropout.lambda)

    return(params)
}
