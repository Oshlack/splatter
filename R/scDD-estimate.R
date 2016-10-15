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
#' params <- scDDEstimate(sc_example_counts)
#' params
#' @export
scDDEstimate <- function(counts, conditions, params = newSCDDParams()) {
    UseMethod("scDDEstimate")
}

#' @rdname scDDEstimate
#' @export
scDDEstimate.SCESet <- function(counts, conditions, params = newSCDDParams()) {
    counts <- scater::counts(counts)
    scDDEstimate(counts, params)
}

#' @rdname scDDEstimate
#' @export
scDDEstimate.matrix <- function(counts, conditions, params = newSCDDParams()) {

    checkmate::assertClass(params, "SCDDParams")
    checkmate::assertIntegerish(conditions, len = ncol(counts), lower = 1,
                                upper = 2)

    counts.list <- list(Cond1 = counts[, conditions == 1],
                        Cond2 = counts[, conditions == 2])

    processed <- scDD::preprocess(counts.list, c("Cond1", "Cond2"),
                                  median_norm = TRUE)

    names(conditions) <- colnames(processed)
    pheno <- as(data.frame(condition = conditions), "AnnotatedDataFrame")

    SCDat <- Biobase::ExpressionSet(assayData = processed, phenoData = pheno)

    params <- setParams(params, nCells = dim(SCDat)[2], SCDat = SCDat)

    return(params)
}