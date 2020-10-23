#' SparseDC simulation
#'
#' Simulate counts from cluster in two conditions using the SparseDC method.
#'
#' @param params SparseDCParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around
#' \code{\link[SparseDC]{sim_data}} that takes a
#' \code{\link{SparseDCParams}}, runs the simulation then converts the
#' output from log-expression to counts and returns a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. The original
#' simulated log-expression values are returned in the \code{LogExprs} assay.
#' See \code{\link[SparseDC]{sim_data}} and the SparseDC paper for
#' more details about how the simulation works.
#'
#' @return SingleCellExperiment containing simulated counts
#'
#' @references
#' Campbell K, Yau C. Uncovering genomic trajectories with heterogeneous genetic
#' and environmental backgrounds across single-cells and populations. bioRxiv
#' (2017).
#'
#' Barron M, Zhang S, Li J. A sparse differential clustering algorithm for
#' tracing cell type changes via single-cell RNA-sequencing data.
#' Nucleic Acids Research (2017).
#'
#' Paper: \url{10.1093/nar/gkx1113}
#'
#' @examples
#' if (requireNamespace("SparseDC", quietly = TRUE)) {
#'     sim <- sparseDCSimulate()
#' }
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
sparseDCSimulate <- function(params = newSparseDCParams(),
                             sparsify = TRUE, verbose = TRUE, ...) {

    checkmate::assertClass(params, "SparseDCParams")
    params <- setParams(params, ...)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    markers.n <- getParam(params, "markers.n")
    markers.shared <- getParam(params, "markers.shared")
    markers.same <- getParam(params, "markers.same")
    clusts.c1 <- getParam(params, "clusts.c1")
    clusts.c2 <- getParam(params, "clusts.c2")
    mean.lower <- getParam(params, "mean.lower")
    mean.upper <- getParam(params, "mean.upper")

    if (verbose) {message("Simulating counts...")}
    sparsedc.sim <- SparseDC::sim_data(genes = nGenes,
                                       cells = nCells,
                                       sig.genes = markers.n,
                                       sig.genes.s = markers.shared,
                                       clus.t1 = clusts.c1,
                                       clus.t2 = clusts.c2,
                                       same.sig = markers.same,
                                       u.l = mean.lower,
                                       u.h = mean.upper)

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(2 * nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    exprs <- cbind(sparsedc.sim$dat.1, sparsedc.sim$dat.2)
    counts <- exp(exprs) - 1
    counts[counts < 0] <- 0
    counts <- round(counts)
    rownames(counts) <- gene.names
    colnames(counts) <- cell.names

    cells <- data.frame(Cell = cell.names,
                        Condition = factor(paste0("Condition",
                                                  rep(seq_len(2),
                                                      each = nCells))),
                        Cluster = factor(paste0("Cluster",
                                                c(sparsedc.sim$clusters1,
                                                  sparsedc.sim$clusters2))),
                        stringsAsFactors = FALSE)
    rownames(cells) <- cell.names

    features <- data.frame(Gene = gene.names,
                           BaseLogMean = sparsedc.sim$gene.means,
                           stringsAsFactors = FALSE)

    for (i in seq_len(ncol(sparsedc.sim$sig.gene.mat.1))) {
        col.name <- paste0("Condition1Cluster", i, "Marker")
        features[[col.name]] <- sparsedc.sim$sig.gene.mat.1[, i] == 1
    }

    for (i in seq_len(ncol(sparsedc.sim$sig.gene.mat.2))) {
        col.name <- paste0("Condition2Cluster", i, "Marker")
        features[[col.name]] <- sparsedc.sim$sig.gene.mat.2[, i] == 1
    }

    for (i in seq_len(ncol(sparsedc.sim$clus.gene.means))) {
        col.name <- paste0("Cluster", i, "LogMean")
        features[[col.name]] <- sparsedc.sim$clus.gene.means[, i]
    }

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
