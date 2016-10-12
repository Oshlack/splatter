#' Lun simulation
#'
#' Simulate single-cell RNA-seq count data using the method described in Lun,
#' Bach and Marioni "Pooling across cells to normalize single-cell RNA
#' sequencing data with many zero counts".
#'
#' @param params splatParams object containing Lun simulation parameters.
#' @param verbose logical. Whether to print progress messages.
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
#' \code{lunSim} requires the following parameters:
#' \itemize{
#'   \item Common (see \code{\link{splatParams}})
#'     \itemize{
#'       \item nGenes
#'       \item nCells
#'       \item nGroups
#'       \item groupCells
#'       \item mean.shape
#'       \item mean.rate
#'       \item bcv.common
#'       \item de.prob
#'       \item de.downprob
#'     }
#'   \item Specific
#'     \itemize{
#'       \item [lun.upFC] - Fold change for up-regulated genes.
#'       \item [lun.downFC] - Fold change for down-regulated genes.
#'     }
#' }
#'
#' See \code{\link{defaultLunParams}} for example of these parameters.
#' Parameters are set in the tiered manner described in \code{\link{splat}}.
#'
#'
#'
#' @return SCESet object containing the simulated counts and intermediate values
#'
#' @references
#' Lun ATL, Bach K, Marioni JC. Pooling across cells to normalize single-cell
#' RNA sequencing data with many zero counts. Genome Biology (2016)
#'
#' Paper: \url{dx.doi.org/10.1186/s13059-016-0947-7}
#'
#' Code: \url{https://github.com/MarioniLab/Deconvolution2016}
#'
#' @examples
#' sim <- simLun()
#'
#' @importFrom Biobase fData fData<- pData pData<- assayData
#' @importFrom scater newSCESet
#' @importFrom stats rnorm rgamma rnbinom
#' @export
simLun <- function(params = defaultLunParams(), verbose = TRUE, ...) {

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- mergeParams(params, defaultLunParams())
    params <- expandParams(params)

    # Set random seed
    seed <- getParams(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nGenes <- getParams(params, "nGenes")
    nCells <- getParams(params, "nCells")
    nGroups <- getParams(params, "nGroups")
    groupCells <- getParams(params, "groupCells")
    mean.shape <- getParams(params, "mean.shape")
    mean.rate <- getParams(params, "mean.rate")
    dispersion <- getParams(params, "bcv.common")
    de.prob <- getParams(params, "de.prob")
    de.downProb <- getParams(params, "de.downProb")
    de.upFC <- getParams(params, "lun.upFC")
    de.downFC <- getParams(params, "lun.downFC")

    if (verbose) {message("Calculating alternative parameters...")}
    nGenes.de <- round(de.prob * nGenes)
    de.upProp <- 1 - de.downProb
    if (verbose) {
        message("Number of DE genes: ", paste(nGenes.de, collapse = ", "))
        message("DE up per group: ", paste(de.upProp, collapse = ", "))
    }

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
        for (idx in 1:nGroups) {
            groups[[idx]] <- rep(paste0("Group", idx), groupCells[idx])

            cell.facs.group <- 2 ^ rnorm(groupCells[idx], sd = 0.5)
            cell.facs[[idx]] <- cell.facs.group

            chosen <- nGenes.de[idx] * (idx - 1) + seq_len(nGenes.de[idx])
            is.up <- seq_len(nGenes.de[idx] * de.upProp[idx])
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

    if (verbose) {message("Simulating counts...")}
    counts <- matrix(rnbinom(nGenes * nCells, mu = cell.means,
                             size = 1 / dispersion),
                     nrow = nGenes, ncol = nCells)

    if (verbose) {message("Creating final SCESet...")}
    cell.names <- paste0("Cell", 1:nCells)
    gene.names <- paste0("Gene", 1:nGenes)
    rownames(counts) <- gene.names
    colnames(counts) <- cell.names

    phenos <- new("AnnotatedDataFrame",
                  data = data.frame(Cell = cell.names, CellFac = cell.facs))
    rownames(phenos) <- cell.names
    features <- new("AnnotatedDataFrame",
                    data = data.frame(Gene = gene.names, GeneMean = gene.means))
    rownames(features) <- gene.names
    sim <- newSCESet(countData = counts, phenoData = phenos,
                     featureData = features)

    rownames(cell.means) <- gene.names
    colnames(cell.means) <- cell.names
    assayData(sim)$CellMeans <- cell.means

    if (nGroups > 1) {
        pData(sim)$Group <- groups
        for (idx in seq_along(de.facs)) {
            fData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs[[idx]]
            fData(sim)[[paste0("GeneMeanGroup", idx)]] <- gene.means *
                                                            de.facs[[idx]]
        }
    }

    return(sim)
}
