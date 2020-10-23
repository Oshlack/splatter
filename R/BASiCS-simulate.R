#' BASiCS simulation
#'
#' Simulate counts using the BASiCS method.
#'
#' @param params BASiCSParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around \code{\link[BASiCS]{BASiCS_Sim}} that
#' takes a \code{\link{BASiCSParams}}, runs the simulation then converts the
#' output to a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object.
#' See \code{\link[BASiCS]{BASiCS_Sim}} for more details of how the simulation
#' works.
#'
#' @return SingleCellExperiment containing simulated counts
#'
#' @references
#' Vallejos CA, Marioni JC, Richardson S. BASiCS: Bayesian Analysis of
#' Single-Cell Sequencing data. PLoS Computational Biology (2015).
#'
#' Paper: \url{10.1371/journal.pcbi.1004333}
#'
#' Code: \url{https://github.com/catavallejos/BASiCS}
#'
#' @examples
#' if (requireNamespace("BASiCS", quietly = TRUE)) {
#'     sim <- BASiCSSimulate()
#' }
#' @export
BASiCSSimulate <- function(params = newBASiCSParams(), sparsify = TRUE,
                           verbose = TRUE, ...) {

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

    if (nSpikes == 0 && nBatches == 1) {
        stop("If there are no spikes there must be multiple batches")
    }

    # Sample gene.params if necessary
    if (nrow(gene.params) != nGenes) {
        warning("Number of gene.params not equal to nGenes, ",
                "gene.params will be sampled.")
        selected <- sample(nrow(gene.params), nGenes, replace = TRUE)
        gene.params <- gene.params[selected, ]
    }

    mu <- gene.params$Mean
    delta <- gene.params$Delta

    has.spikes <- nSpikes > 0
    spike.mu <- NULL
    if (has.spikes) {
        spike.mu <- getParam(params, "spike.means")
        if (length(spike.mu) != nSpikes) {
            warning("Number of spike-in means not equal to nSpikes, ",
                    "spike.means will be sampled.")
            selected <- sample(length(spike.mu), nSpikes, replace = TRUE)
            spike.mu <- spike.mu[selected]

            params <- setParam(params, "spike.mu", spike.mu)
        }
    }

    cell.params <- getParam(params, "cell.params")
    if (nrow(cell.params) != nCells) {
        warning("Number of cell.params not equal to nCells, ",
                "cell.params will be sampled.")
        selected <- sample(nrow(cell.params), nCells, replace = TRUE)
        cell.params <- cell.params[selected, ]

        cell.params$Phi <- (cell.params$Phi / sum(cell.params$Phi)) * nCells

        params <- setParam(params, "cell.params", cell.params)
    }

    thetas <- getParam(params, "theta")

    batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                      b = batch.cells)
    batches <- unlist(batches)

    if (verbose) {message("Simulating counts with BASiCS...")}
    if (has.spikes) {
        BASiCS.sim <- suppressMessages(
            BASiCS::BASiCS_Sim(Mu = mu, Mu_spikes = spike.mu, Delta = delta,
                               Phi = cell.params$Phi, S = cell.params$S,
                               Theta = thetas, BatchInfo = batches)
        )
        counts <- SingleCellExperiment::counts(BASiCS.sim)
        spike.counts <- SummarizedExperiment::assay(
            SingleCellExperiment::altExp(BASiCS.sim)
        )
        counts <- rbind(counts, spike.counts)
    } else {
        if (verbose) {message("Simulating without spikes")}
        BASiCS.sim <- suppressMessages(
            BASiCS::BASiCS_Sim(Mu = mu, Mu_spikes = NULL, Delta = delta,
                               Phi = NULL, S = cell.params$S,
                               Theta = thetas, BatchInfo = batches)
        )
        counts <- SingleCellExperiment::counts(BASiCS.sim)
    }

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))
    if (has.spikes) {
        gene.names <- c(gene.names, paste0("Spike", seq_len(nSpikes)))
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
                                metadata = list(Params = params))

    if (sparsify) {
        if (verbose) {message("Sparsifying assays...")}
        assays(sim) <- sparsifyMatrices(assays(sim), auto = TRUE,
                                        verbose = verbose)
    }

    if (verbose) {message("Done!")}
    return(sim)
}
