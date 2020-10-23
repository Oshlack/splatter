#' PhenoPath simulation
#'
#' Simulate counts from a pseudotime trajectory using the PhenoPath method.
#'
#' @param params PhenoParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around
#' \code{\link[phenopath]{simulate_phenopath}} that takes a
#' \code{\link{PhenoParams}}, runs the simulation then converts the
#' output from log-expression to counts and returns a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. The original
#' simulated log-expression values are returned in the \code{LogExprs} assay.
#' See \code{\link[phenopath]{simulate_phenopath}} and the PhenoPath paper for
#' more details about how the simulation works.
#'
#' @return SingleCellExperiment containing simulated counts
#'
#' @references
#' Campbell K, Yau C. Uncovering genomic trajectories with heterogeneous genetic
#' and environmental backgrounds across single-cells and populations. bioRxiv
#' (2017).
#'
#' Paper: \url{10.1101/159913}
#'
#' Code: \url{https://github.com/kieranrcampbell/phenopath}
#'
#' @examples
#' if (requireNamespace("phenopath", quietly = TRUE)) {
#'     sim <- phenoSimulate()
#' }
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
phenoSimulate <- function(params = newPhenoParams(), sparsify = TRUE,
                          verbose = TRUE, ...) {

    checkmate::assertClass(params, "PhenoParams")
    params <- setParams(params, ...)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    n.de <- getParam(params, "n.de")
    n.pst <- getParam(params, "n.pst")
    n.pst.beta <- getParam(params, "n.pst.beta")
    n.de.pst.beta <- getParam(params, "n.de.pst.beta")

    if (verbose) {message("Simulating counts...")}
    pheno.sim <- phenopath::simulate_phenopath(N = nCells,
                                               G_de = n.de,
                                               G_pst = n.pst,
                                               G_pst_beta = n.pst.beta,
                                               G_de_pst_beta = n.de.pst.beta)

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    exprs <- t(pheno.sim$y)
    counts <- 2 ^ exprs - 1
    counts[counts < 0] <- 0
    counts <- round(counts)
    rownames(counts) <- gene.names
    colnames(counts) <- cell.names

    cells <- data.frame(Cell = cell.names,
                        Covariate = pheno.sim$x,
                        Pseudotime = pheno.sim$z)
    rownames(cells) <- cell.names

    features <- data.frame(Gene = gene.names,
                           Alpha = pheno.sim$parameters$alpha,
                           Lambda = pheno.sim$parameters$lambda,
                           Beta = pheno.sim$parameters$beta,
                           Regime = pheno.sim$parameters$regime)
    rownames(features) <- gene.names

    sim <- SingleCellExperiment(assays = list(counts = counts,
                                              LogExprs = exprs),
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
