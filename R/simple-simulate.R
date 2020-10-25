#' Simple simulation
#'
#' Simulate counts from a simple negative binomial distribution without
#' simulated library sizes, differential expression etc.
#'
#' @param params SimpleParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' Gene means are simulated from a gamma distribution with
#' \code{shape = mean.shape} and \code{rate = mean.rate}. Counts are then
#' simulated from a negative binomial distribution with \code{mu = means} and
#' \code{size = 1 / counts.disp}. See \code{\link{SimpleParams}} for more
#' details of the parameters.
#'
#' @return SingleCellExperiment containing simulated counts
#' @examples
#' sim <- simpleSimulate()
#' # Override default parameters
#' sim <- simpleSimulate(nGenes = 1000, nCells = 50)
#' @export
#' @importFrom stats rgamma rnbinom
#' @importFrom SingleCellExperiment SingleCellExperiment
simpleSimulate <- function(params = newSimpleParams(), sparsify = TRUE,
                           verbose = TRUE, ...) {

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
    counts <- matrix(rnbinom(
            as.numeric(nGenes) * as.numeric(nCells),
            mu = means, size = 1 / count.disp),
    nrow = nGenes, ncol = nCells)

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    rownames(counts) <- gene.names
    colnames(counts) <- cell.names
    cells <- data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names, GeneMean = means)
    rownames(features) <- gene.names

    sim <- SingleCellExperiment(assays = list(counts = counts),
                                rowData = features,
                                colData = cells,
                                metadata = list(Params = params))

    if (sparsify) {
        if (verbose) {message("Sparsifying assays...")}
        assays(sim) <- sparsifyMatrices(assays(sim), auto = TRUE,
                                        verbose = verbose)
    }

    return(sim)
}
