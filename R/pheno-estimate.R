#' Estimate PhenoPath simulation parameters
#'
#' Estimate simulation parameters for the PhenoPath simulation from a real
#' dataset.
#'
#' @param counts either a counts matrix or an SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param params PhenoParams object to store estimated values in.
#'
#' @details
#' The \code{nGenes} and \code{nCells} parameters are taken from the size of the
#' input data. The total number of genes is evenly divided into the four types.
#' See \code{\link{PhenoParams}} for more details on the parameters.
#'
#' @return PhenoParams object containing the estimated parameters.
#'
#' @examples
#' if (requireNamespace("phenopath", quietly = TRUE)) {
#'     # Load example data
#'     library(scater)
#'     set.seed(1)
#'     sce <- mockSCE()
#'
#'     params <- phenoEstimate(sce)
#'     params
#' }
#' @export
phenoEstimate <- function(counts, params = newPhenoParams()) {
    UseMethod("phenoEstimate")
}

#' @rdname phenoEstimate
#' @export
phenoEstimate.SingleCellExperiment <- function(counts,
                                               params = newPhenoParams()) {
    counts <- getCounts(counts)
    phenoEstimate(counts, params)
}

#' @rdname phenoEstimate
#' @export
phenoEstimate.matrix <- function(counts, params = newPhenoParams()) {

    checkmate::assertClass(params, "PhenoParams")

    nGenes <- nrow(counts)
    quarter <- floor(nGenes / 4)

    params <- setParams(params, nCells = ncol(counts),
                        n.de = nGenes - 3 * quarter,
                        n.pst = quarter, n.pst.beta = quarter,
                        n.de.pst.beta = quarter)

    return(params)
}
