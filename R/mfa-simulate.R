#' MFA simulation
#'
#' Simulate a bifurcating pseudotime path using the mfa method.
#'
#' @param params MFAParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose Logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' This function is just a wrapper around \code{\link[mfa]{create_synthetic}}
#' that takes a \code{\link{MFAParams}}, runs the simulation then converts the
#' output from log-expression to counts and returns a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. See
#' \code{\link[mfa]{create_synthetic}} and the mfa paper for more details about
#' how the simulation works.
#'
#' @return SingleCellExperiment containing simulated counts
#'
#' @references
#' Campbell KR, Yau C. Probabilistic modeling of bifurcations in single-cell
#' gene expression data using a Bayesian mixture of factor analyzers. Wellcome
#' Open Research (2017).
#'
#' Paper: \url{10.12688/wellcomeopenres.11087.1}
#'
#' Code: \url{https://github.com/kieranrcampbell/mfa}
#'
#' @examples
#' if (requireNamespace("mfa", quietly = TRUE)) {
#'     sim <- mfaSimulate()
#' }
#' @export
mfaSimulate <- function(params = newMFAParams(), sparsify = TRUE,
                        verbose = TRUE, ...) {

    checkmate::assertClass(params, "MFAParams")
    params <- setParams(params, ...)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    trans.prop <- getParam(params, "trans.prop")
    zero.neg <- getParam(params, "zero.neg")
    dropout.present <- getParam(params, "dropout.present")
    dropout.lambda <- getParam(params, "dropout.lambda")

    if (verbose) {message("Simulating counts...")}
    mfa.sim <- mfa::create_synthetic(C = nCells,
                                     G = nGenes,
                                     p_transient = trans.prop,
                                     zero_negative = zero.neg,
                                     model_dropout = dropout.present,
                                     lambda = dropout.lambda)

    if (verbose) {message("Creating final dataset...")}
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    exprs <- t(mfa.sim$X)
    rownames(exprs) <- gene.names
    colnames(exprs) <- cell.names

    counts <- 2 ^ exprs - 1
    counts[counts < 0] <- 0
    counts <- round(counts)

    cells <- data.frame(Cell = cell.names,
                        Branch = mfa.sim$branch,
                        Pseudotime = mfa.sim$pst)
    rownames(cells) <- cell.names

    features <- data.frame(Gene = gene.names,
                           KBranch1 = mfa.sim$k[, 1],
                           KBranch2 = mfa.sim$k[, 2],
                           PhiBranch1 = mfa.sim$phi[, 1],
                           PhiBranch2 = mfa.sim$phi[, 2],
                           DeltaBranch1 = mfa.sim$delta[, 1],
                           DeltaBranch2 = mfa.sim$delta[, 2])
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
