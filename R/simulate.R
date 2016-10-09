# setup metadata
# Means
# Groups
# Paths
# Lib size
# Base Means
# BCV
# Means
# Counts
# Dropout
# Add metadata
#
# Add length
# Median outliers
# set.seed
# group DE



#' FUNCTION TITLE
#'
#' FUNCTION DESCRIPTION
#'
#' @param params DESCRIPTION.
#' @param method DESCRIPTION.
#' @param add.assay DESCRIPTION.
#' @param verbose DESCRIPTION.
#' @param ... DESCRIPTION.
#'
#' @return RETURN DESCRIPTION
#' @examples
#' # ADD EXAMPLES HERE
#' @importFrom Biobase pData fData
#' @importFrom scater newSCESet counts
splat <- function(params = defaultParams(), method = c("groups", "paths"),
                  add.assay = TRUE, verbose = TRUE, ...) {

    method <- match.arg(method)

    params <- setParams(params, ...)
    params <- mergeParams(params, defaultParams())

    # Get the parameters we are going to use
    n.cells <- getParams(params, "nCells")
    n.genes <- getParams(params, "nGenes")
    n.groups <- getParams(params, "nGroups")
    group.cells <- getParams(params, "groupCells")

    # Set up name vectors
    cell.names <- paste0("Cell", 1:n.cells)
    gene.names <- paste0("Gene", 1:n.genes)
    group.names <- paste0("Group", 1:n.groups)

    # Create SCESet with dummy counts to store simulation
    dummy.counts <- matrix(1, ncol = n.cells, nrow = n.genes)
    rownames(dummy.counts) <- gene.names
    colnames(dummy.counts) <- cell.names
    phenos <- new("AnnotatedDataFrame", data = data.frame(Cell = cell.names))
    rownames(phenos) <- cell.names
    features <- new("AnnotatedDataFrame", data = data.frame(Gene = gene.names))
    rownames(features) <- gene.names
    sim <- newSCESet(countData = dummy.counts, phenoData = phenos,
                     featureData = features)

    # Make groups vector which is the index of param$groupCells repeated
    # params$groupCells[index] times
    groups <- lapply(1:n.groups, function(i, g) {rep(i, g[i])},
                     g = group.cells)
    groups <- unlist(groups)
    pData(sim)$Group <- group.names[groups]

    sim <- simulateGeneMeans(sim, params)
    sim <- simulateDE(sim, params)

    # Create new SCESet to make sure values are calculated correctly
    sce <- newSCESet(countData = counts(sim),
                     phenoData = new("AnnotatedDataFrame", data = pData(sim)),
                     featureData = new("AnnotatedDataFrame", data = fData(sim)))

    return(sce)
}

simulateGeneMeans <- function(sim, params) {

    n.genes <- getParams(params, "nGenes")
    mean.shape <- getParams(params, "mean.shape")
    mean.rate <- getParams(params, "mean.rate")
    out.prob <- getParams(params, "out.prob")
    out.loProb <- getParams(params, "out.loProb")
    out.facLoc <- getParams(params, "out.facLoc")
    out.facScale <- getParams(params, "out.facScale")

    # Simulate base gene means
    base.means.gene <- rgamma(n.genes, shape = mean.shape, rate = mean.rate)

    # Add expression outliers
    outlier.facs <- getLNormFactors(n.genes, out.prob, out.loProb, out.facLoc,
                                    out.facScale)
    means.gene <- base.means.gene * outlier.facs

    fData(sim)$BaseGeneMean <- base.means.gene
    fData(sim)$OutlierFactor <- outlier.facs
    fData(sim)$GeneMean <- means.gene

    return(sim)
}

simulateDE <- function(sim, params) {

    n.genes <- getParams(params, "nGenes")
    de.prob <- getParams(params, "de.prob")
    de.downProb <- getParams(params, "de.downProb")
    de.facLoc <- getParams(params, "de.facLoc")
    de.facScale <- getParams(params, "de.facScale")
    means.gene <- fData(sim)$GeneMean
    group.names <- unique(pData(sim)$Group)

    for (group.name in group.names) {
        de.facs <- getLNormFactors(n.genes, de.prob, de.downProb, de.facLoc,
                                   de.facScale)
        group.means.gene <- means.gene * de.facs
        fData(sim)[[paste0("DEFac", group.name)]] <- de.facs
        fData(sim)[[paste0("GeneMean", group.name)]] <- group.means.gene
    }

    return(sim)
}

#' Get log-normal factors
#'
#' Randomly generate multiplication factors from a log-normal distribution.
#'
#' @param n.facs Number of factors to generate.
#' @param sel.prob Probability that a factor will be selected to be different
#'        from 1.
#' @param neg.prob Probability that a selected factor is less than one.
#' @param fac.loc Location parameter for the log-normal distribution.
#' @param fac.scale Scale factor for the log-normal distribution.
#'
#' @return Vector containing generated factors.
#' @examples
#' factors <- getLNormFactors(100, 0.5, 0.5, 4, 1)
getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale) {

    is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
    n.selected <- sum(is.selected)
    dir.selected <- (-1) ^ rbinom(n.selected, 1, neg.prob)
    facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
    # Reverse directions for factors that are less than one
    dir.selected[facs.selected < 1 & dir.selected == -1] <- 1
    dir.selected[facs.selected < 1 & dir.selected == 1] <- -1
    factors <- rep(1, n.facs)
    factors[is.selected] <- facs.selected ^ dir.selected

    return(factors)
}