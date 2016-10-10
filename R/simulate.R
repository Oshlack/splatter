# setup metadata x
# Means x
# Groups x
# Paths
# Lib size x
# Base Means X
# BCV X
# Means X
# Counts X
# Dropout *********
# Add metadata X
#
# Add length
# Median outliers
# set.seed
# group DE

# Min path length
# Path.from

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
#' @importFrom Biobase pData fData assayData
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
    if (method == "groups") {
        group.names <- paste0("Group", 1:n.groups)
    } else if (method == "paths") {
        group.names <- paste0("Path", 1:n.groups)
    }

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

    sim <- simLibSizes(sim, params)
    sim <- simGeneMeans(sim, params)
    if (method == "groups") {
        sim <- simGroupDE(sim, params)
        sim <- simGroupCellMeans(sim, params)
    } else {
        sim <- simPathDE(sim, params)
    }
    #sim <- simBCVMeans(sim, params)
    #sim <- simTrueCounts(sim, params)
    #im <- simDropout(sim, params)

    # Create new SCESet to make sure values are calculated correctly
    sce <- newSCESet(countData = counts(sim),
                     phenoData = new("AnnotatedDataFrame", data = pData(sim)),
                     featureData = new("AnnotatedDataFrame", data = fData(sim)))

    # Add intermediate matrices stored in assayData
    for (assay.name in names(assayData(sim))) {
        if (!(assay.name %in% names(assayData(sce)))) {
            assayData(sce)[[assay.name]] <- assayData(sim)[[assay.name]]
        }
    }

    return(sce)
}

simLibSizes <- function(sim, params) {

    n.cells <- getParams(params, "nCells")
    lib.loc <- getParams(params, "lib.loc")
    lib.scale <- getParams(params, "lib.scale")

    exp.lib.sizes <- rlnorm(n.cells, lib.loc, lib.scale)
    pData(sim)$ExpLibSize <- exp.lib.sizes

    return(sim)
}

simGeneMeans <- function(sim, params) {

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

simGroupDE <- function(sim, params) {

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

simPathDE <- function(sim, params) {

    n.genes <- getParams(params, "nGenes")
    de.prob <- getParams(params, "de.prob")
    de.downProb <- getParams(params, "de.downProb")
    de.facLoc <- getParams(params, "de.facLoc")
    de.facScale <- getParams(params, "de.facScale")
    path.from <- getParams(params, "path.from")
    means.gene <- fData(sim)$GeneMean
    path.names <- unique(pData(sim)$Group)

    path.order <- getPathOrder(path.from)
    for (path.name in path.names[path.order]) {
        de.facs <- getLNormFactors(n.genes, de.prob, de.downProb, de.facLoc,
                                   de.facScale)
        path.means.gene <- means.gene * de.facs
        fData(sim)[[paste0("DEFac", path.name)]] <- de.facs
        fData(sim)[[paste0("GeneMean", path.name)]] <- path.means.gene
    }

    return(sim)
}

simGroupCellMeans <- function(sim, params) {

    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    groups <- pData(sim)$Group
    group.names <- unique(groups)
    exp.lib.sizes <- pData(sim)$ExpLibSize

    group.means.gene <- fData(sim)[, paste0("GeneMean", group.names)]
    cell.means.gene <- as.matrix(group.means.gene[, factor(groups)])
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names

    assayData(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

simBCVMeans <- function(sim, params) {

    n.genes <- getParams(params, "nGenes")
    n.cells <- getParams(params, "nCells")
    bcv.common <- getParams(params, "bcv.common")
    bcv.DF <- getParams(params, "bcv.DF")
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    base.means.cell <- assayData(sim)$BaseCellMeans

    bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
        sqrt(bcv.DF / rchisq(n.genes, df = bcv.DF))

    means.cell <- matrix(rgamma(n.genes * n.cells, shape = 1 / (bcv ^ 2),
                                scale = base.means.cell * (bcv ^ 2)),
                         nrow = n.genes, ncol = n.cells)

    colnames(bcv) <- cell.names
    rownames(bcv) <- gene.names
    colnames(means.cell) <- cell.names
    rownames(means.cell) <- gene.names

    assayData(sim)$BCV <- bcv
    assayData(sim)$CellMeans <- means.cell

    return(sim)
}

simTrueCounts <- function(sim, params) {

    n.genes <- getParams(params, "nGenes")
    n.cells <- getParams(params, "nCells")
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    cell.means <- assayData(sim)$CellMeans

    true.counts <- matrix(rpois(n.genes * n.cells, lambda = cell.means),
                          nrow = n.genes, ncol = n.cells)

    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names

    assayData(sim)$TrueCounts <- true.counts

    return(sim)
}

simDropout <- function(sim, params) {

    dropout.present <- getParams(params, "dropout.present")
    true.counts <- assayData(sim)$TrueCounts

    if (dropout.present) {
        n.cells <- getParams(params, "nCells")
        n.genes <- getParams(params, "nGenes")
        dropout.mid <- getParams(params, "dropout.mid")
        dropout.shape <- getParams(params, "dropout.shape")
        cell.names <- pData(sim)$Cell
        gene.names <- fData(sim)$Gene
        cell.means <- assayData(sim)$CellMeans

        lib.sizes <- colSums(true.counts)
        cell.facs <- log(lib.sizes) / median(lib.sizes)
        drop.prob <- sapply(1:n.cells, function(idx) {
            eta <- cell.facs[idx] * (log(cell.means[, idx]))
            return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
        })

        keep <- matrix(rbinom(n.cells * n.genes, 1, 1 - drop.prob),
                       nrow = n.genes, ncol = n.cells)

        counts <- true.counts * keep

        colnames(drop.prob) <- cell.names
        rownames(drop.prob) <- gene.names
        colnames(keep) <- cell.names
        rownames(keep) <- gene.names

        assayData(sim)$DropProb <- drop.prob
        assayData(sim)$Dropout <- !keep

    } else {
        counts <- true.counts
    }

    counts(sim) <- counts

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

#' Get path order
#'
#' Identify the correct order to process paths so that preceding paths have
#' already been simulated.
#'
#' @param path.from Vector giving the path endpoints that each path orginates
#'        from.
#'
#' @return Vector giving the order to process paths in.
#' @examples
#' path.order <- getPathOrder(c(2, 0, 2))
getPathOrder <- function(path.from) {

    # Transform the vector into a list of (from, to) pairs
    path.pairs <- list()
    for (idx in 1:length(path.from)) {
        path.pairs[[idx]] <- c(path.from[idx], idx)
    }

    # Determine the processing order
    # If a path is in the "done" vector any path originating here can be
    # completed
    done <- 0
    while (length(path.pairs) > 0) {
        path.pair <- path.pairs[[1]]
        path.pairs <- path.pairs[-1]
        from <- path.pair[1]
        to <- path.pair[2]
        if (from %in% done) {
            done <- c(done, to)
        } else {
            path.pairs <- c(path.pairs, list(path.pair))
        }
    }

    # Remove the origin from the vector
    done <- done[-1]

    return(done)
}