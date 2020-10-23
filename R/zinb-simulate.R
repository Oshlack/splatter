#' ZINB-WaVE simulation
#'
#' Simulate counts using the ZINB-WaVE method.
#'
#' @param params ZINBParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around \code{\link[zinbwave]{zinbSim}} that
#' takes a \code{\link{ZINBParams}}, runs the simulation then converts the
#' output to a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' See \code{\link[zinbwave]{zinbSim}} and the ZINB-WaVE paper for
#' more details about how the simulation works.
#'
#' @return SingleCellExperiment containing simulated counts
#'
#' @references
#' Campbell K, Yau C. Uncovering genomic trajectories with heterogeneous genetic
#' and environmental backgrounds across single-cells and populations. bioRxiv
#' (2017).
#'
#' Risso D, Perraudeau F, Gribkova S, Dudoit S, Vert J-P. ZINB-WaVE: A general
#' and flexible method for signal extraction from single-cell RNA-seq data
#' bioRxiv (2017).
#'
#' Paper: \url{10.1101/125112}
#'
#' Code: \url{https://github.com/drisso/zinbwave}
#'
#' @examples
#' if (requireNamespace("zinbwave", quietly = TRUE)) {
#'     sim <- zinbSimulate()
#' }
#'
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
zinbSimulate <- function(params = newZINBParams(), sparsify = TRUE,
                         verbose = TRUE, ...) {

    checkmate::assertClass(params, "ZINBParams")
    params <- setParams(params, ...)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    model <- getParam(params, "model")
    seed <- getParam(params, "seed")

    if (verbose) {message("Simulating counts...")}
    zinb.sim <- zinbwave::zinbSim(model, seed)

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    for (item in c("counts", "dataNB", "dataDropouts")) {
        rownames(zinb.sim[[item]]) <- gene.names
        colnames(zinb.sim[[item]]) <- cell.names
    }

    cells <- data.frame(Cell = cell.names)
    rownames(cells) <- cell.names

    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names

    sim <- SingleCellExperiment(assays = list(counts = zinb.sim$counts,
                                              TrueCounts = zinb.sim$dataNB,
                                              Dropouts = zinb.sim$dataDropouts),
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
