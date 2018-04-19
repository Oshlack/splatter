#' Estimate SparseDC simulation parameters
#'
#' Estimate simulation parameters for the SparseDC simulation from a real
#' dataset.
#'
#' @param counts either a counts matrix or an SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param conditions numeric vector giving the condition each cell belongs to.
#' @param nclusters number of cluster present in the dataset.
#' @param norm logical, whether to libray size normalise counts before
#'        estimation. Set this to FALSE if counts is already normalised.
#' @param params PhenoParams object to store estimated values in.
#'
#' @details
#' The \code{nGenes} and \code{nCells} parameters are taken from the size of the
#' input data. The counts are preprocessed using
#' \code{\link[SparseDC]{pre_proc_data}} and then parameters are estimated using
#' \code{\link[SparseDC]{sparsedc_cluster}} using lambda values calculated using
#' \code{\link[SparseDC]{lambda1_calculator}} and
#' \code{\link[SparseDC]{lambda2_calculator}}.
#'
#' See \code{\link{SparseDCParams}} for more details on the parameters.
#'
#' @return SparseParams object containing the estimated parameters.
#'
#' @examples
#' # Load example data
#' library(scater)
#' data("sc_example_counts")
#'
#' conditions <- sample(1:2, ncol(sc_example_counts), replace = TRUE)
#'
#' params <- sparseDCEstimate(sc_example_counts[1:500, ], conditions,
#'                            nclusters = 3)
#' params
#' @export
sparseDCEstimate <- function(counts, conditions, nclusters, norm = TRUE,
                             params = newSparseDCParams()) {
    UseMethod("sparseDCEstimate")
}

#' @rdname sparseDCEstimate
#' @export
sparseDCEstimate.SingleCellExperiment <- function(counts, conditions, nclusters,
                                                  norm = TRUE,
                                                  params = newSparseDCParams()) {
    counts <- BiocGenerics::counts(counts)
    sparseDCEstimate(counts, params)
}

#' @rdname sparseDCEstimate
#' @export
sparseDCEstimate.matrix <- function(counts, conditions, nclusters, norm = TRUE,
                                    params = newSparseDCParams()) {

    checkmate::assertClass(params, "SparseDCParams")
    checkmate::assertIntegerish(conditions, lower = 1, upper = 2,
                                any.missing = FALSE, len = ncol(counts))

    counts1 <- counts[, conditions == 1]
    counts2 <- counts[, conditions == 2]

    pre.data <- SparseDC::pre_proc_data(counts1, counts2, norm = norm,
                                        log = TRUE, center = TRUE)

    lambda1 <- SparseDC::lambda1_calculator(pre.data[[1]], pre.data[[2]],
                                            nclusters)
    lambda2 <- SparseDC::lambda2_calculator(pre.data[[1]], pre.data[[2]],
                                            nclusters)

    dummy <- utils::capture.output(
    sdc.res <- SparseDC::sparsedc_cluster(pre.data[[1]], pre.data[[2]],
                                          nclusters, lambda1 = lambda1,
                                          lambda2 = lambda2)
    )

    markers.n <- round(mean(c(
        colSums(sdc.res$centers1 != 0),
        colSums(sdc.res$centers2 != 0)
    )))

    markers.diff <- round(mean(colSums(
        (sdc.res$centers1 - sdc.res$centers2) != 0
    )))

    params <- setParams(params,
                        nGenes = nrow(counts),
                        nCells = round(ncol(counts) / 2),
                        markers.n = markers.n,
                        markers.shared = markers.n - markers.diff,
                        clusts.c1 = sort(unique(sdc.res$clusters1)),
                        clusts.c2 = sort(unique(sdc.res$clusters2)))

    return(params)
}
