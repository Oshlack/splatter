#' Estimate scDD simulation parameters
#'
#' Estimate simulation parameters for the scDD simulation from a real dataset.
#'
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param conditions Vector giving the condition that each cell belongs to.
#'        Conditions can be 1 or 2.
#' @param condition String giving the column that represents biological group of
#'        interest.
#' @param params SCDDParams object to store estimated values in.
#' @param verbose logical. Whether to show progress messages.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} instance
#'        giving the parallel back-end to be used. Default is
#'        \code{\link[BiocParallel]{SerialParam}} which uses a single core.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' This function applies \code{\link[scDD]{preprocess}} to the counts then uses
#' \code{\link[scDD]{scDD}} to estimate the numbers of each gene type to
#' simulate. The output is then converted to a SCDDParams object. See
#' \code{\link[scDD]{preprocess}} and \code{\link[scDD]{scDD}} for details.
#'
#' @return SCDDParams object containing the estimated parameters.
#'
#' @examples
#' \dontrun{
#' data("sc_example_counts")
#' conditions <- sample(1:2, ncol(sc_example_counts), replace = TRUE)
#' params <- scDDEstimate(sc_example_counts, conditions = conditions)
#' params
#' }
#' @importFrom BiocParallel SerialParam
#' @export
scDDEstimate <- function(counts, #conditions, condition,
                         params = newSCDDParams(), verbose = TRUE,
                         BPPARAM = SerialParam(), ...) {

    if (!requireNamespace("scDD", quietly = TRUE)) {
        stop("The scDD simulation requires the 'scDD' package.")
    }

    UseMethod("scDDEstimate")
}

#' @rdname scDDEstimate
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
scDDEstimate.matrix <- function(counts, params = newSCDDParams(),
                                verbose = TRUE, BPPARAM = SerialParam(),
                                conditions, ...) {

    checkmate::assertMatrix(counts, mode = "numeric", any.missing = FALSE,
                            min.rows = 1, min.cols = 1, row.names = "unique",
                            col.names = "unique")
    checkmate::assertIntegerish(conditions, len = ncol(counts), lower = 1,
                                upper = 2)

    counts <- SingleCellExperiment(assays = list(counts = counts),
                                   colData = data.frame(condition = conditions))
    scDDEstimate.default(counts, params, verbose, BPPARAM,
                         condition = "condition")
}

#' @rdname scDDEstimate
#' @export
scDDEstimate.SingleCellExperiment <- function(counts,
                                              params = newSCDDParams(),
                                              verbose = TRUE,
                                              BPPARAM = SerialParam(),
                                              condition = "condition", ...) {
    scDDEstimate.default(counts, params, verbose, BPPARAM,
                         condition = condition)
}

#' @rdname scDDEstimate
#' @importFrom methods as
#' @export
scDDEstimate.default <- function(counts,
                                 params = newSCDDParams(), verbose = TRUE,
                                 BPPARAM = SerialParam(), condition, ...) {

    checkmate::assertClass(params, "SCDDParams")
    checkmate::assertClass(counts, "SingleCellExperiment")
    checkmate::assertCharacter(condition, min.chars = 1, any.missing = FALSE,
                               len = 1)
    if (!(condition %in% colnames(SummarizedExperiment::colData(counts)))) {
        stop("'condition' must be the name of a column in `colData(counts)`")
    }

    if (verbose) {
        processed <- scDD::preprocess(counts, condition, median_norm = TRUE)
    } else {
        suppressMessages(
        processed <- scDD::preprocess(counts, condition, median_norm = TRUE)
        )
    }

    if (verbose) {
        SCdat <- scDD::scDD(processed, testZeroes = FALSE, param = BPPARAM,
                            condition = condition)
    } else {
        dummy <- utils::capture.output(suppressMessages(
        SCdat <- scDD::scDD(processed, testZeroes = FALSE, param = BPPARAM,
                            condition = condition)
        ))
    }

    res <- scDD::results(SCdat)
    res <- res[!is.na(res$DDcategory), ]
    dd.cats <- table(res$DDcategory)

    not.dd <- res$DDcategory == "NS"
    nDE <- ifelse("DE" %in% names(dd.cats), dd.cats["DE"], 0)
    nDP <- ifelse("DP" %in% names(dd.cats), dd.cats["DP"], 0)
    nDM <- ifelse("DM" %in% names(dd.cats), dd.cats["DM"], 0)
    nDB <- ifelse("DB" %in% names(dd.cats), dd.cats["DB"], 0)
    nEP <- sum(res$Clusters.c1[not.dd] > 1 & res$Clusters.c2[not.dd] > 1)
    nEE <- nrow(counts) - nDE - nDP - nDM - nDB - nEP

    params <- setParams(params,
                        nCells = round(dim(SCdat)[2] / 2),
                        SCdat = SCdat,
                        nDE = nDE,
                        nDP = nDP,
                        nDM = nDM,
                        nDB = nDB,
                        nEE = nEE,
                        nEP = nEP)

    return(params)
}
