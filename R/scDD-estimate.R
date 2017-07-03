#' Estimate scDD simulation parameters
#'
#' Estimate simulation parameters for the scDD simulation from a real dataset.
#'
#' @param counts either a counts matrix or an SCESet object containing count
#'        data to estimate parameters from.
#' @param conditions Vector giving the condition that each cell belongs to.
#'        Conditions can be 1 or 2.
#' @param params SCDDParams object to store estimated values in.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} instance
#'        giving the parallel back-end to be used. Default is
#'        \code{\link[BiocParallel]{SerialParam}} which uses a single core.
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
#' data("sc_example_counts")
#' conditions <- sample(1:2, ncol(sc_example_counts), replace = TRUE)
#' params <- scDDEstimate(sc_example_counts, conditions)
#' params
#' @importFrom BiocParallel SerialParam
#' @export
scDDEstimate <- function(counts, conditions, params = newSCDDParams(),
                         BPPARAM = SerialParam()) {
    UseMethod("scDDEstimate")
}

#' @rdname scDDEstimate
#' @export
scDDEstimate.SCESet <- function(counts, conditions, params = newSCDDParams(),
                                BPPARAM = SerialParam()) {
    counts <- scater::counts(counts)
    scDDEstimate(counts, conditions, params)
}

#' @rdname scDDEstimate
#' @importFrom methods as
#' @export
scDDEstimate.matrix <- function(counts, conditions, params = newSCDDParams(),
                                BPPARAM = SerialParam()) {

    if (!requireNamespace("scDD", quietly = TRUE)) {
        stop("The scDD simulation requires the 'scDD' package.")
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

    SCdat <- scDD::scDD(SCdat, testZeroes = FALSE, param = BPPARAM)

    res <- scDD::results(SCdat)
    res <- res[!is.na(res$DDcategory), ]
    dd.cats <- table(res$DDcategory)

    not.dd <- res$DDcategory == "NC" | res$DDcategory == "NS"
    nDE <- ifelse("DE" %in% names(dd.cats), dd.cats["DE"], 0)
    nDP <- ifelse("DP" %in% names(dd.cats), dd.cats["DP"], 0)
    nDM <- ifelse("DM" %in% names(dd.cats), dd.cats["DM"], 0)
    nDB <- ifelse("DB" %in% names(dd.cats), dd.cats["DB"], 0)
    nEP <- sum(res$Clusters.c1[not.dd] > 1 | res$Clusters.c2[not.dd] > 1)
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
