#' Lun simulation
#'
#' Simulate single-cell RNA-seq count data using the method described in Lun,
#' Bach and Marioni "Pooling across cells to normalize single-cell RNA
#' sequencing data with many zero counts".
#'
#' @param params LunParams object containing Lun simulation parameters.
#' @param verbose logical. Whether to print progress messages.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' The Lun simulation generates gene mean expression levels from a gamma
#' distribution with \code{shape = mean.shape} and \code{rate = mean.rate}.
#' Counts are then simulated from a negative binomial distribution with
#' \code{mu = means} and \code{size = 1 / bcv.common}. In addition each cell is
#' given a size factor (\code{2 ^ rnorm(nCells, mean = 0, sd = 0.5)}) and
#' differential expression can be simulated with fixed fold changes.
#'
#' See \code{\link{LunParams}} for details of the parameters.
#'
#' @return SingleCellExperiment object containing the simulated counts and
#' intermediate values.
#'
#' @references
#' Lun ATL, Bach K, Marioni JC. Pooling across cells to normalize single-cell
#' RNA sequencing data with many zero counts. Genome Biology (2016).
#'
#' Paper: \url{dx.doi.org/10.1186/s13059-016-0947-7}
#'
#' Code: \url{https://github.com/MarioniLab/Deconvolution2016}
#'
#' @examples
#' sim <- lunSimulate()
#'
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom stats rnorm rgamma rnbinom
#' @export
lunSimulate <- function(params = newLunParams(), sparsify = TRUE,
                        verbose = TRUE, ...) {

    checkmate::assertClass(params, "LunParams")

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- expandParams(params)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    nGroups <- getParam(params, "nGroups")
    groupCells <- getParam(params, "groupCells")
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    count.disp <- getParam(params, "count.disp")
    de.nGenes <- getParam(params, "de.nGenes")
    de.upProp <- getParam(params, "de.upProp")
    de.upFC <- getParam(params, "de.upFC")
    de.downFC <- getParam(params, "de.downFC")

    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    if (verbose) {message("Simulating means...")}
    gene.means <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

    if (verbose) {message("Simulating cell means...")}
    if (nGroups == 1) {
        cell.facs <- 2 ^ rnorm(nCells, sd = 0.5)
        cell.means <- outer(gene.means, cell.facs, "*")
    } else {
        groups <- list()
        cell.facs <- list()
        de.facs <- list()
        cell.means <- list()
        for (idx in seq_len(nGroups)) {
            groups[[idx]] <- rep(paste0("Group", idx), groupCells[idx])

            cell.facs.group <- 2 ^ rnorm(groupCells[idx], sd = 0.5)
            cell.facs[[idx]] <- cell.facs.group

            chosen <- de.nGenes[idx] * (idx - 1) + seq_len(de.nGenes[idx])
            is.up <- seq_len(de.nGenes[idx] * de.upProp[idx])
            de.up <- chosen[is.up]
            de.down <- chosen[-is.up]

            de.facs.group <- rep(1, nGenes)
            de.facs.group[de.up] <- de.upFC[idx]
            de.facs.group[de.down] <- de.downFC[idx]
            de.facs[[idx]] <- de.facs.group

            cell.means.group <- outer(gene.means, cell.facs.group)
            cell.means.group <- cell.means.group * de.facs.group
            cell.means[[idx]] <- cell.means.group
        }
        cell.means <- do.call(cbind, cell.means)
        cell.facs <- unlist(cell.facs)
        groups <- unlist(groups)
    }
    colnames(cell.means) <- cell.names
    rownames(cell.means) <- gene.names

    if (verbose) {message("Simulating counts...")}
    counts <- matrix(rnbinom(
            as.numeric(nGenes) * as.numeric(nCells),
            mu = cell.means, size = 1 / count.disp
        ),
    nrow = nGenes, ncol = nCells)

    if (verbose) {message("Creating final dataset...")}
    rownames(counts) <- gene.names
    colnames(counts) <- cell.names

    cells <- data.frame(Cell = cell.names, CellFac = cell.facs)
    rownames(cells) <- cell.names

    features <- data.frame(Gene = gene.names, GeneMean = gene.means)
    rownames(features) <- gene.names

    if (nGroups > 1) {
        cells$Group <- groups
        for (idx in seq_along(de.facs)) {
            features[[paste0("DEFacGroup", idx)]] <- de.facs[[idx]]
            features[[paste0("GeneMeanGroup", idx)]] <- gene.means *
                de.facs[[idx]]
        }
    }

    sim <- SingleCellExperiment(assays = list(counts = counts,
                                              CellMeans = cell.means),
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
