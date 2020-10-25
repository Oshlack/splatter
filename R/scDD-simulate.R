#' scDD simulation
#'
#' Simulate counts using the scDD method.
#'
#' @param params SCDDParams object containing simulation parameters.
#' @param plots logical. whether to generate scDD fold change and validation
#' plots.
#' @param plot.file File path to save plots as PDF.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#' @param BPPARAM A \code{\link{BiocParallelParam}} instance giving the parallel
#'        back-end to be used. Default is \code{\link{SerialParam}} which uses a
#'        single core.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around \code{\link[scDD]{simulateSet}} that
#' takes a \code{\link{SCDDParams}}, runs the simulation then converts the
#' output to a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' See \code{\link[scDD]{simulateSet}} for more details about how the simulation
#' works.
#'
#' @return SingleCellExperiment containing simulated counts
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
#' @importFrom BiocParallel SerialParam
#' @importFrom SingleCellExperiment SingleCellExperiment
scDDSimulate <- function(params = newSCDDParams(), plots = FALSE,
                         plot.file = NULL, sparsify = TRUE, verbose = TRUE,
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
        suppressMessages({
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
        })
    }

    counts <- SummarizedExperiment::assays(scDD.sim)$normcounts
    foldchanges <- SummarizedExperiment::rowData(scDD.sim)$FC
    de.status <- SummarizedExperiment::rowData(scDD.sim)$Category

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(nCells * 2))
    gene.names <- paste0("Gene", seq_len(getParam(params, "nGenes")))

    rownames(counts) <- gene.names
    colnames(counts) <- cell.names

    cells <- data.frame(Cell = cell.names,
                        Condition = rep(seq_len(2), each = nCells))
    rownames(cells) <- cell.names

    features <- data.frame(Gene = gene.names, DEStatus = de.status,
                           FoldChange = foldchanges)
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

    if (verbose) {message("Done!")}

    return(sim)
}
