#' BASiCS simulation
#'
#' Simulate counts using the BASiCS method.
#'
#' @param params BASiCSParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around \code{\link[BASiCS]{BASiCS_Sim}} that
#' takes a \code{\link{BASiCSParams}}, runs the simulation then converts the
#' output to an \code{\link[scater]{SCESet}} object. See
#' \code{\link[BASiCS]{BASiCS_Sim}} for more details of how the simulation
#' works.
#'
#' @return SingleCellExperiment containing simulated counts
#'
#' @references
#' Vallejos CA, Marioni JC, Richardson S. BASiCS: Bayesian Analysis of
#' Single-Cell Sequencing data. PLoS Comput. Biol. (2015).
#'
#' Paper: \url{10.1371/journal.pcbi.1004333}
#'
#' Code: \url{https://github.com/catavallejos/BASiCS}
#'
#' @examples
#' sim <- BASiCSSimulate()
#' @export
#' @importFrom scater newSCESet
BASiCSSimulate <- function(params = newBASiCSParams(), verbose = TRUE, ...) {

    checkmate::assertClass(params, "BASiCSParams")
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    if (verbose) {message("Getting parameters...")}
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    nSpikes <- getParam(params, "nSpikes")
    nBatches <- getParam(params, "nBatches")
    batch.cells <- getParam(params, "batchCells")
    gene.params <- getParam(params, "gene.params")

    # Sample gene.params if necessary
    if (nrow(gene.params) != nGenes) {
        warning("Number of gene.params not equal to nGenes, ",
                "gene.params will be sampled.")
        selected <- sample(nrow(gene.params), nGenes, replace = TRUE)
        gene.params <- gene.params[selected, ]
    }

    mu <- gene.params$Mean
    delta <- gene.params$Delta

    if (nSpikes > 0) {
        spike.mu <- getParam(params, "spike.means")
        if (length(spike.mu) != nSpikes) {
            warning("Number of spike-in means not equal to nSpikes, ",
                    "spike.means will be sampled.")
            selected <- sample(length(spike.mu), nSpikes, replace = TRUE)
            spike.mu <- spike.mu[selected]
        }
    } else {
        # Create dummy spike-ins to get around BASiCS_Sim...
        spike.mu <- c(10, 10)
    }

    cell.params <- getParam(params, "cell.params")
    if (nrow(cell.params) != nCells) {
        warning("Number of cell.params not equal to nCells, ",
                "cell.params will be sampled.")
        selected <- sample(nrow(cell.params), nCells, replace = TRUE)
        cell.params <- cell.params[selected, ]
    }

    thetas <- getParam(params, "theta")

    batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                      b = batch.cells)
    batches <- unlist(batches)

    if (verbose) {message("Simulating counts with BASiCS...")}
    counts.list <- list()
    for (batch in seq_len(nBatches)) {
        batch.cells <- batches == batch
        phi <- cell.params[batch.cells, "Phi"]
        phi <- (phi / sum(phi)) * sum(batch.cells)
        s <- cell.params[batch.cells, "S"]
        theta <- thetas[batch]
        BASiCS.sim <- suppressMessages(
                          BASiCS::BASiCS_Sim(mu, spike.mu, delta, phi, s, theta)
                      )
        batch.counts <- SummarizedExperiment::assay(BASiCS.sim)
        counts.list[[batch]] <- batch.counts
    }

    counts <- do.call(cbind, counts.list)

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))
    if (nSpikes > 0) {
        gene.names <- c(gene.names, paste0("Spike", seq_len(nSpikes)))
    } else {
        # Remove dummy spikes
        counts <- counts[1:(nrow(counts) - 2), ]
        spike.mu <- numeric()
    }

    rownames(counts) <- gene.names
    colnames(counts) <- cell.names

    cells <- data.frame(Cell = cell.names,
                        Phi = cell.params[, "Phi"],
                        S = cell.params[, "S"],
                        Batch = batches,
                        BatchTheta = thetas[batches])
    rownames(cells) <- cell.names

    features <- data.frame(Gene = gene.names,
                           Mean = c(mu, spike.mu),
                           Delta = c(delta, rep(NA, nSpikes)),
                           IsSpike = c(rep(FALSE, nGenes), rep(TRUE, nSpikes)))
    rownames(features) <- gene.names

    sim <- SingleCellExperiment(assays = list(counts = counts),
                                rowData = features,
                                colData = cells,
                                metadata = list(params = params))

    if (verbose) {message("Done!")}
    return(sim)
}
