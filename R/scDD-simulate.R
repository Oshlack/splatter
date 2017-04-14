#' scDD simulation
#'
#' Simulate counts using the scDD method.
#'
#' @param params SCDDParams object containing simulation parameters.
#' @param plots logical. whether to generate scDD fold change and validation
#' plots.
#' @param plot.file File path to save plots as PDF.
#' @param verbose logical. Whether to print progress messages
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} instance
#'        giving the parallel back-end to be used. Default is
#'        \code{\link[BiocParallel]{SerialParam}} which uses a single core.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around \code{\link[scDD]{simulateSet}} that
#' takes a \code{\link{SCDDParams}}, runs the simulation then converts the
#' output to an \code{\link[scater]{SCESet}} object. See
#' \code{\link[scDD]{simulateSet}} for more details of how the simulation works.
#'
#' @return SCESet containing simulated counts
#'
#' @references
#' Korthauer KD, Chu L-F, Newton MA, Li Y, Thomson J, Stewart R, et al. A
#' statistical approach for identifying differential distributions in
#' single-cell RNA-seq experiments. Genome Biology (2016).
#'
#' Paper: \url{10.1186/s13059-016-1077-y}
#'
#' Code: \url{https://github.com/kdkorthauer/scDD}
#'
#' @examples
#' \dontrun{
#' sim <- scDDSimulate()
#' }
#' @export
#' @importFrom scater newSCESet
#' @importFrom BiocParallel SerialParam
scDDSimulate <- function(params = newSCDDParams(), plots = FALSE,
                         plot.file = NULL, verbose = TRUE,
                         BPPARAM = SerialParam(), ...) {

    checkmate::assertClass(params, "SCDDParams")
    params <- setParams(params, ...)

    nCells <- getParam(params, "nCells")

    if (verbose) {message("Simulating counts with scDD...")}
    # Restore the default varInflation to NULL
    varInflation <- getParam(params, "varInflation")
    if (all(varInflation == 1)) {
        varInflation <- NULL
    }
    if (verbose) {
        scDD.sim <- scDD::simulateSet(SCdat = getParam(params, "SCdat"),
                                      numSamples = nCells,
                                      nDE = getParam(params, "nDE"),
                                      nDP = getParam(params, "nDP"),
                                      nDM = getParam(params, "nDM"),
                                      nDB = getParam(params, "nDB"),
                                      nEE = getParam(params, "nEE"),
                                      nEP = getParam(params, "nEP"),
                                      sd.range = getParam(params, "sd.range"),
                                      modeFC = getParam(params, "modeFC"),
                                      plots = plots,
                                      plot.file = plot.file,
                                      random.seed = getParam(params, "seed"),
                                      varInflation = varInflation,
                                      condition = getParam(params, "condition"),
                                      param = BPPARAM)
    } else {
        suppressMessages(
        scDD.sim <- scDD::simulateSet(SCdat = getParam(params, "SCdat"),
                                      numSamples = nCells,
                                      nDE = getParam(params, "nDE"),
                                      nDP = getParam(params, "nDP"),
                                      nDM = getParam(params, "nDM"),
                                      nDB = getParam(params, "nDB"),
                                      nEE = getParam(params, "nEE"),
                                      nEP = getParam(params, "nEP"),
                                      sd.range = getParam(params, "sd.range"),
                                      modeFC = getParam(params, "modeFC"),
                                      plots = plots,
                                      plot.file = plot.file,
                                      random.seed = getParam(params, "seed"),
                                      varInflation = varInflation,
                                      condition = getParam(params, "condition"),
                                      param = BPPARAM)
        )
    }

    counts <- scDD.sim[[1]]
    foldchanges <- scDD.sim[[2]]
    de.status <- rownames(counts)

    if (verbose) {message("Creating SCESet...")}
    cell.names <- paste0("Cell", seq_len(nCells * 2))
    gene.names <- paste0("Gene", seq_len(getParam(params, "nGenes")))

    rownames(counts) <- gene.names
    colnames(counts) <- cell.names
    phenos <- new("AnnotatedDataFrame",
                  data = data.frame(Cell = cell.names,
                                    Condition = rep(1:2, each = nCells)))
    rownames(phenos) <- cell.names
    features <- new("AnnotatedDataFrame",
                    data = data.frame(Gene = gene.names, DEStatus = de.status,
                                      FoldChange = foldchanges))
    rownames(features) <- gene.names
    sim <- newSCESet(countData = counts, phenoData = phenos,
                     featureData = features)

    if (verbose) {message("Done!")}

    return(sim)
}
