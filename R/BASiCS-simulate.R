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
#' @return SCESet containing simulated counts
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

    if (verbose) {message("Getting parameters...")}
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    nSpikes <- getParam(params, "nSpikes")

    gene.params <- getParam(params, "gene.params")
    # Sample gene.params if necessary
    if (nrow(gene.params) != nGenes) {
        warning("Number of gene.params not equal to nGenes, ",
                "gene.params will be sampled.")
        selected <- sample(nrow(gene.params), nGenes, replace = TRUE)
        gene.params <- gene.params[selected, ]
    }

    cell.params <- getParam(params, "cell.params")
    if (nrow(cell.params) != nCells) {
        warning("Number of cell.params not equal to nCells, ",
                "cell.params will be sampled.")
        selected <- sample(nrow(cell.params), nCells, replace = TRUE)
        cell.params <- cell.params[selected, ]
    }

    mu <- gene.params$Mean

    if (nSpikes > 0) {
        spike.mu <- getParam(params, "spike.means")
        if (length(spike.mu) != nSpikes) {
            warning("Number of spike-in means not equal to nSpikes, ",
                    "spike.means will be sampled.")
            selected <- sample(length(spike.mu), nSpikes, replace = TRUE)
            spike.mu <- spike.mu[selected]
        }
    }

    delta <- gene.params$Delta

    phi <- cell.params$Phi
    if (!(sum(phi) == nCells)) {
        warning("cell.params$Phi rescaled to sum to nCells")
        phi <- (phi / sum(phi)) * nCells
    }

    s <- cell.params$S
    theta <- getParam(params, "theta")

    if (verbose) {message("Simulating counts with BASiCS...")}
    BASiCS.sim <- suppressMessages(BASiCS::BASiCS_Sim(mu, spike.mu, delta, phi,
                                                      s, theta))

    counts <- assay(BASiCS)

    if (verbose) {message("Creating SCESet...")}
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))
    if (nSpikes > 0) {
        gene.names <- c(gene.names, paste0("Spike", seq_len(nSpikes)))
    }

    rownames(counts) <- gene.names
    colnames(counts) <- cell.names
    phenos <- new("AnnotatedDataFrame",
                  data = data.frame(Cell = cell.names,
                                    Phi = phi,
                                    S = s))
    rownames(phenos) <- cell.names
    features <- data.frame(Gene = gene.names,
                           Mean = c(mu, spike.mu),
                           Delta = c(delta, rep(NA, nSpikes)),
                           IsSpike = c(rep(FALSE, nGenes), rep(TRUE, nSpikes)))
    features <- new("AnnotatedDataFrame",
                    data = features)
    rownames(features) <- gene.names
    sim <- newSCESet(countData = counts, phenoData = phenos,
                     featureData = features)

    if (verbose) {message("Done!")}
    return(sim)
}
