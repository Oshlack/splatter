#' Estimate scDD simulation parameters
#'
#' Estimate simulation parameters for the scDD simulation from a real dataset.
#'
#' @param counts either a counts matrix or an SCESet object containing count
#'        data to estimate parameters from.
#' @param conditions Vector giving the condition that each cell belongs to.
#'        Conditions can be 1 or 2.
#' @param params SCDDParams object to store estimated values in.
#'
#' @details
#' This function is just a wrapper around \code{\link[scDD]{preprocess}} that
#' takes the output and converts it to a SCDDParams object. See
#' \code{\link[scDD]{preprocess}} for details.
#'
#' @return SCDDParams object containing the estimated parameters.
#'
#' @examples
#' data("sc_example_counts")
#' conditions <- sample(1:2, ncol(sc_example_counts), replace = TRUE)
#' params <- scDDEstimate(sc_example_counts, conditions)
#' params
#' @export
scDDEstimate <- function(counts, conditions, params = newSCDDParams()) {
    UseMethod("scDDEstimate")
}

#' @rdname scDDEstimate
#' @export
scDDEstimate.SCESet <- function(counts, conditions, params = newSCDDParams()) {
    counts <- scater::counts(counts)
    scDDEstimate(counts, conditions, params)
}

#' @rdname scDDEstimate
#' @importFrom methods as
#' @export
scDDEstimate.matrix <- function(counts, conditions, params = newSCDDParams()) {

    if (!requireNamespace("scDD", quietly = TRUE)) {
        stop("The scDD simulation requires the 'scDD' package. ",
             "See https://github.com/kdkorthauer/scDD for installation.")
    }

    checkmate::assertClass(params, "SCDDParams")
    checkmate::assertIntegerish(conditions, len = ncol(counts), lower = 1,
                                upper = 2)

    counts.list <- list(Cond1 = counts[, conditions == 1],
                        Cond2 = counts[, conditions == 2])

    processed <- scDD::preprocess(counts.list, c("Cond1", "Cond2"),
                                  median_norm = TRUE)

    assays <- S4Vectors::SimpleList(NormCounts = processed)

    colData <- S4Vectors::DataFrame(condition = conditions,
                                    row.names = colnames(processed))

    SCdat <- SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                        colData = colData)

    params <- setParams(params, nCells = round(dim(SCdat)[2] / 2),
                        SCdat = SCdat)

    return(params)
}