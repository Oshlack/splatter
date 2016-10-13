#' Simple simulation
#'
#' Simulate counts from a simple negative binomial distribution without
#' simulated library sizes, differential expression etc.
#'
#' @param params splatParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' Uses the following parameters: \code{nCells}, \code{nGenes},
#' \code{mean.shape}, \code{mean.rate}, \code{bcv.common}.
#'
#' Gene means are simulated from a gamma distribution with
#' \code{shape = mean.shape} and \code{rate = mean.rate}. Counts are then
#' simulated from a negative binomial distribution with \code{mu = means} and
#' \code{size = 1 / bcv.common}.
#'
#' Parameters are set in the tiered manner described in \code{\link{splat}}.
#'
#' @return SCESet containing simulated counts
#' @examples
#' sim <- simSimple()
#' @export
#' @importFrom stats rgamma rnbinom
simpleSimulate <- function(params = newSimpleParams(), verbose = TRUE, ...) {

    checkmate::assertClass(params, "SimpleParams")
    params <- setParams(params, ...)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    count.disp <- getParam(params, "count.disp")

    if (verbose) {message("Simulating means...")}
    means <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

    if (verbose) {message("Simulating counts...")}
    counts <- matrix(rnbinom(nGenes * nCells, mu = means,
                             size = 1 / count.disp),
                     nrow = nGenes, ncol = nCells)

    if (verbose) {message("Creating final SCESet...")}
    cell.names <- paste0("Cell", 1:nCells)
    gene.names <- paste0("Gene", 1:nGenes)

    rownames(counts) <- gene.names
    colnames(counts) <- cell.names
    phenos <- new("AnnotatedDataFrame", data = data.frame(Cell = cell.names))
    rownames(phenos) <- cell.names
    features <- new("AnnotatedDataFrame",
                    data = data.frame(Gene = gene.names, GeneMean = means))
    rownames(features) <- gene.names
    sim <- newSCESet(countData = counts, phenoData = phenos,
                     featureData = features)

    return(sim)
}