# Add length
# Median outliers
# set.seed
# group DE

# Min path length
# Path.from + tests
# Bridge tests
# Verbose
# Dropout

#' Simulate scRNA-seq data
#'
#' Simulate count data from a fictional single-cell RNA-seq experiment.
#'
#' @param params splatParams object containing parameters for the simulation.
#'        See \code{\link{splatParams}} for details.
#' @param method which simulation method to use. Options are "single" which
#'        produces a single population, "groups" which produces distinct groups
#'        (eg. cell types) or "paths" which selects cells from continuous
#'        trajectories (eg. differentiation process).
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' Parameters can be set in a variety of ways. If no parameters are provided
#' the default parameters are used (see \code{\link{defaultParams}}). Any
#' parameters in \code{params} can be overridden by supplying additional
#' arguments through a call to \code{\link{setParams}}. Finally any parameters
#' the are still missing (\code{NA}) are replaced with the defaults through a
#' call to \code{\link{mergeParams}}. This design allows the user flexibility in
#' how they supply parameters and allows small adjustments without creating a
#' new \code{splatParams} object. See examples for a demonstration of how this
#' can be used.
#'
#' The simulation involves the following steps:
#' \enumerate{
#'   \item Set up simulation object
#'   \item Simulate library sizes
#'   \item Simulate gene means
#'   \item Simulate groups/paths
#'   \item Simulate BCV adjusted cell means
#'   \item Simulate true counts
#'   \item Simulate dropout
#'   \item Create final SCESet object
#' }
#'
#' The final output is an \code{\link{SCESet}} object that contains the
#' simulated counts but also the values for various intermediate steps. These
#' are stored in the \code{\link{phenoData}} (for cell specific information),
#' \code{\link{featureData}} (for gene specific information) or
#' \code{\link{assayData}} (for gene by cell matrices) slots. This additional
#' information includes:
#' \itemize{
#'   \item \code{phenoData}
#'     \itemize{
#'       \item Cell - unique cell identifier
#'       \item Group - the group or path the cell belongs to
#'       \item ExpLibSize - the expected library size for that cell
#'       \item Step (paths only) - how far along the path each cell is
#'     }
#'   \item \code{featureData}
#'     \itemize{
#'       \item Gene - unique gene identifier
#'       \item BaseGeneMean - the base expression level for that gene
#'       \item OutlierFactor - expression outlier factor for that gene. Values
#'             of 1 indicate the gene is not an expression outlier.
#'       \item GeneMean - expression level after applying outlier factors.
#'       \item DEFac[Group] - the differential expression factor for each gene
#'             in a particular group. Values of 1 indicate the gene is not
#'             differentially expressed.
#'       \item GeneMean[Group] - expression level of a gene in a particular
#'             group after applying differential expression factors.
#'     }
#'   \item \code{assayData}
#'     \itemize{
#'       \item BaseCellMeans - the expression of genes in each cell adjusted for
#'             expected library size.
#'       \item BCV - the Biological Coefficient of Variation for each gene in
#'             each cell.
#'       \item CellMeans - the expression level of genes in each cell adjusted
#'             for BCV.
#'       \item TrueCounts - the simulated counts before dropout.
#'       \item Dropout - logical matrix showing which values have been dropped
#'             in which cells.
#'     }
#'  }
#'
#'  Values that have been added by Splatter are named using \code{CamelCase} in
#'  order to differentiate them from the values added by Scater which uses
#'  \code{underscore_naming}.
#'
#' @return SCESet object containing the simulated counts and intermediate values
#' @examples
#' # Simulation with default parameters
#' sim <- splat()
#' # Simulation with different number of genes
#' sim <- splat(nGenes = 1000)
#' # Simulation with custom parameters
#' params <- splatParams(nGenes = 100, mean.rate = 0.5)
#' sim <- splat(params)
#' # Simulation with adjusted custom parameters
#' sim <- splat(params, mean.rate = 0.6, out.prob = 0.2)
#' # Simulate paths instead of groups
#' sim <- splat(method = "paths")
#' @importFrom Biobase fData pData pData<- assayData
#' @importFrom scater newSCESet counts
#' @export
splat <- function(params = defaultParams(),
                  method = c("single", "groups", "paths"),
                  verbose = TRUE, ...) {

    method <- match.arg(method)

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- mergeParams(params, defaultParams())
    params <- expandParams(params)

    # Set random seed
    seed <- getParams(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParams(params, "nCells")
    nGenes <- getParams(params, "nGenes")
    nGroups <- getParams(params, "nGroups")
    group.cells <- getParams(params, "groupCells")

    if (nGroups == 1 && method == "groups") {
        warning("nGroups is 1, switching to single mode")
        method <- "single"
    }

    if (verbose) {message("Creating simulation object...")}
    # Set up name vectors
    cell.names <- paste0("Cell", 1:nCells)
    gene.names <- paste0("Gene", 1:nGenes)
    if (method == "groups") {
        group.names <- paste0("Group", 1:nGroups)
    } else if (method == "paths") {
        group.names <- paste0("Path", 1:nGroups)
    }

    # Create SCESet with dummy counts to store simulation
    dummy.counts <- matrix(1, ncol = nCells, nrow = nGenes)
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
    if (method != "single") {
        groups <- lapply(1:nGroups, function(i, g) {rep(i, g[i])},
                         g = group.cells)
        groups <- unlist(groups)
        pData(sim)$Group <- group.names[groups]
    }

    if (verbose) {message("Simulating library sizes...")}
    sim <- simLibSizes(sim, params)
    if (verbose) {message("Simulating gene means...")}
    sim <- simGeneMeans(sim, params)
    if (method == "single") {
        sim <- simSingleCellMeans(sim, params)
    } else if (method == "groups") {
        if (verbose) {message("Simulating group DE...")}
        sim <- simGroupDE(sim, params)
        if (verbose) {message("Simulating cell means...")}
        sim <- simGroupCellMeans(sim, params)
    } else {
        if (verbose) {message("Simulating path endpoints...")}
        sim <- simPathDE(sim, params)
        if (verbose) {message("Simulating path steps...")}
        sim <- simPathCellMeans(sim, params)
    }
    if (verbose) {message("Simulating BCV...")}
    sim <- simBCVMeans(sim, params)
    if (verbose) {message("Simulating counts..")}
    sim <- simTrueCounts(sim, params)
    if (verbose) {message("Simulating dropout...")}
    sim <- simDropout(sim, params)

    if (verbose) {message("Creating final SCESet...")}
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

    if (verbose) {message("Done!")}
    return(sce)
}

#' @rdname splat
#' @export
splatSingle <- function(params = defaultParams(), verbose = TRUE, ...) {
    sim <- splat(params = params, method = "single", verbose = verbose, ...)
    return(sim)
}

#' @rdname splat
#' @export
splatGroups <- function(params = defaultParams(), verbose = TRUE, ...) {
    sim <- splat(params = params, method = "groups", verbose = verbose, ...)
    return(sim)
}

#' @rdname splat
#' @export
splatPaths <- function(params = defaultParams(), verbose = TRUE, ...) {
    sim <- splat(params = params, method = "paths", verbose = verbose, ...)
    return(sim)
}

#' Simulate library sizes
#'
#' Simulate expected library sizes from a log-normal distribution
#'
#' @param sim SCESet to add library size to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added library sizes.
#'
#' @importFrom Biobase pData pData<-
#' @importFrom stats rlnorm
simLibSizes <- function(sim, params) {

    nCells <- getParams(params, "nCells")
    lib.loc <- getParams(params, "lib.loc")
    lib.scale <- getParams(params, "lib.scale")

    exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    pData(sim)$ExpLibSize <- exp.lib.sizes

    return(sim)
}

#' Simulate gene means
#'
#' Simulate gene means from a gamma distribution. Also simulates outlier
#' expression factors. Genes with an outlier factor not equal to 1 are replaced
#' with the median mean expression multiplied by the outlier factor.
#'
#' @param sim SCESet to add gene means to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added gene means.
#'
#' @importFrom Biobase fData fData<-
#' @importFrom stats rgamma median
simGeneMeans <- function(sim, params) {

    nGenes <- getParams(params, "nGenes")
    mean.shape <- getParams(params, "mean.shape")
    mean.rate <- getParams(params, "mean.rate")
    out.prob <- getParams(params, "out.prob")
    out.loProb <- getParams(params, "out.loProb")
    out.facLoc <- getParams(params, "out.facLoc")
    out.facScale <- getParams(params, "out.facScale")

    # Simulate base gene means
    base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

    # Add expression outliers
    outlier.facs <- getLNormFactors(nGenes, out.prob, out.loProb, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    fData(sim)$BaseGeneMean <- base.means.gene
    fData(sim)$OutlierFactor <- outlier.facs
    fData(sim)$GeneMean <- means.gene

    return(sim)
}

#' Simulate group differential expression
#'
#' Simulate differential expression for groups. Differential expression
#' factors for each group are produced using \code{\link{getLNormFactors}} and
#' these are added along with updated means for each group.
#'
#' @param sim SCESet to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added differential expression.
#'
#' @importFrom Biobase fData
simGroupDE <- function(sim, params) {

    nGenes <- getParams(params, "nGenes")
    nGroups <- getParams(params, "nGroups")
    de.prob <- getParams(params, "de.prob")
    de.downProb <- getParams(params, "de.downProb")
    de.facLoc <- getParams(params, "de.facLoc")
    de.facScale <- getParams(params, "de.facScale")
    means.gene <- fData(sim)$GeneMean

    for (idx in 1:nGroups) {
        de.facs <- getLNormFactors(nGenes, de.prob[idx], de.downProb[idx],
                                   de.facLoc[idx], de.facScale[idx])
        group.means.gene <- means.gene * de.facs
        fData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs
        fData(sim)[[paste0("GeneMeanGroup", idx)]] <- group.means.gene
    }

    return(sim)
}

#' Simulate path differential expression
#'
#' Simulate differential expression for path. Similar to
#' \code{\link{simGroupDE}} but care has to be taken to make sure paths are
#' processed in the correct order.
#'
#' @param sim SCESet to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added differential expression.
#'
#' @importFrom Biobase fData
simPathDE <- function(sim, params) {

    nGenes <- getParams(params, "nGenes")
    de.prob <- getParams(params, "de.prob")
    de.downProb <- getParams(params, "de.downProb")
    de.facLoc <- getParams(params, "de.facLoc")
    de.facScale <- getParams(params, "de.facScale")
    path.from <- getParams(params, "path.from")

    path.order <- getPathOrder(path.from)
    for (path in path.order) {
        from <- path.from[path]
        if (from == 0) {
            means.gene <- fData(sim)$GeneMean
        } else {
            means.gene <- fData(sim)[[paste0("GeneMeanPath", from)]]
        }
        de.facs <- getLNormFactors(nGenes, de.prob[path], de.downProb[path],
                                   de.facLoc[path], de.facScale[path])
        path.means.gene <- means.gene * de.facs
        fData(sim)[[paste0("DEFacPath", path)]] <- de.facs
        fData(sim)[[paste0("GeneMeanPath", path)]] <- path.means.gene
    }

    return(sim)
}

#' Simulate single population cell means
#'
#' Simulate a gene by cell matrix giving the mean expression for each gene in
#' each cell.
#'
#' @param sim SCESet to add cell means to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added cell means.
#'
#' @importFrom Biobase fData pData assayData assayData<-
simSingleCellMeans <- function(sim, params) {

    nCells <- getParams(params, "nCells")
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    exp.lib.sizes <- pData(sim)$ExpLibSize

    cell.means.gene <- as.matrix(fData(sim)[, rep("GeneMean", nCells)])
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names

    assayData(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' Simulate group cell means
#'
#' Simulate a gene by cell matrix giving the mean expression for each gene in
#' each cell. Cells start with the mean expression for the group they belong to
#' which is adjusted for each cells expected library size.
#'
#' @param sim SCESet to add cell means to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added cell means.
#'
#' @importFrom Biobase fData pData assayData assayData<-
simGroupCellMeans <- function(sim, params) {

    nGroups <- getParams(params, "nGroups")
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

#' Simulate path cell means
#'
#' Simulate a gene by cell matrix giving the mean expression for each gene in
#' each cell. Cells are assigned assigned a random position on the appropriate
#' path. The mean at that position is then adjusted for each cells expected
#' library size.
#'
#' @param sim SCESet to add cell means to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added cell means.
#'
#' @importFrom Biobase fData pData assayData
#' @importFrom stats rbinom
simPathCellMeans <- function(sim, params) {

    nGenes <- getParams(params, "nGenes")
    nGroups <- getParams(params, "nGroups")
    group.cells <- getParams(params, "groupCells")
    path.from <- getParams(params, "path.from")
    path.length <- getParams(params, "path.length")
    path.skew <- getParams(params, "path.skew")
    path.nonlinearProb <- getParams(params, "path.nonlinearProb")
    path.sigmaFac <- getParams(params, "path.sigmaFac")
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    groups <- pData(sim)$Group
    group.names <- unique(groups)
    exp.lib.sizes <- pData(sim)$ExpLibSize

    # Generate paths. Each path is a matrix with path.length columns and
    # nGenes rows where the expression from each genes changes along the path.
    path.steps <- lapply(seq_along(path.from), function(idx) {
        from <- path.from[idx]
        # Find the means at the starting position
        if (from == 0) {
            means.start <- fData(sim)$GeneMean
        } else {
            means.start <- fData(sim)[[paste0("GeneMeanPath", from)]]
        }
        # Find the means at the end position
        means.end <- fData(sim)[[paste0("GeneMeanPath", idx)]]

        # Select genes to follow a non-linear path
        is.nonlinear <- as.logical(rbinom(nGenes, 1, path.nonlinearProb))
        sigma.facs <- rep(0, nGenes)
        sigma.facs[is.nonlinear] <- path.sigmaFac
        # Build Brownian bridges from start to end
        steps <- buildBridges(means.start, means.end, n = path.length[idx],
                              sigma.fac = sigma.facs)

        fData(sim)[[paste0("SigmaFacPath", idx)]] <- sigma.facs

        return(t(steps))
    })

    # Randomly assign a position in the appropriate path to each cell
    cell.steps <- lapply(1:nGroups, function(idx) {
        path.probs <- seq(path.skew[idx], 1 - path.skew[idx],
                          length = path.length[idx])
        path.probs <- path.probs / sum(path.probs)
        steps <- sort(sample(1:path.length[idx], group.cells[idx],
                             prob = path.probs, replace = TRUE))

        return(steps)
    })

    # Collect the underlying expression levels for each cell
    cell.means.gene <- lapply(1:nGroups, function(idx) {
        cell.means <- path.steps[[idx]][, cell.steps[[idx]]]
        return(cell.means)
    })
    cell.means.gene <- do.call(cbind, cell.means.gene)

    # Adjust expression based on library size
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names

    pData(sim)$Step <- unlist(cell.steps)
    assayData(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' Simulate BCV means
#'
#' Simulate means for each gene in each cell that are adjusted to follow a
#' mean-variance trend using Biological Coefficient of Variation taken from
#' and inverse gamma distribution.
#'
#' @param sim SCESet to add BCV means to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added BCV means.
#'
#' @importFrom Biobase fData pData assayData assayData<-
#' @importFrom stats rchisq rgamma
simBCVMeans <- function(sim, params) {

    nGenes <- getParams(params, "nGenes")
    nCells <- getParams(params, "nCells")
    bcv.common <- getParams(params, "bcv.common")
    bcv.DF <- getParams(params, "bcv.DF")
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    base.means.cell <- assayData(sim)$BaseCellMeans

    bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
        sqrt(bcv.DF / rchisq(nGenes, df = bcv.DF))

    means.cell <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
                                scale = base.means.cell * (bcv ^ 2)),
                         nrow = nGenes, ncol = nCells)

    colnames(bcv) <- cell.names
    rownames(bcv) <- gene.names
    colnames(means.cell) <- cell.names
    rownames(means.cell) <- gene.names

    assayData(sim)$BCV <- bcv
    assayData(sim)$CellMeans <- means.cell

    return(sim)
}

#' Simulate true counts
#'
#' Simulate a true counts matrix. Counts are simulated from a poisson
#' distribution where Each gene in each cell has it's own mean based on the
#' group (or path position), expected library size and BCV.
#'
#' @param sim SCESet to add true counts to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added true counts.
#'
#' @importFrom Biobase fData pData assayData
#' @importFrom stats rpois
simTrueCounts <- function(sim, params) {

    nGenes <- getParams(params, "nGenes")
    nCells <- getParams(params, "nCells")
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    cell.means <- assayData(sim)$CellMeans

    true.counts <- matrix(rpois(nGenes * nCells, lambda = cell.means),
                          nrow = nGenes, ncol = nCells)

    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names

    assayData(sim)$TrueCounts <- true.counts

    return(sim)
}

#' Simulate dropout
#'
#' A logistic function is used to form a relationshop between the expression
#' level of a gene and the probability of dropout, giving a probability for each
#' gene in each cell. These probabilities are used in a Bernoulli distribution
#' to decide which counts should be dropped.
#'
#' @param sim SCESet to add dropout to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with added dropout and observed counts.
#'
#' @importFrom Biobase fData pData assayData assayData<-
#' @importFrom stats rbinom
simDropout <- function(sim, params) {

    dropout.present <- getParams(params, "dropout.present")
    true.counts <- assayData(sim)$TrueCounts

    if (dropout.present) {
        nCells <- getParams(params, "nCells")
        nGenes <- getParams(params, "nGenes")
        dropout.mid <- getParams(params, "dropout.mid")
        dropout.shape <- getParams(params, "dropout.shape")
        cell.names <- pData(sim)$Cell
        gene.names <- fData(sim)$Gene
        cell.means <- assayData(sim)$CellMeans

        # Generate probabilites based on expression
        lib.sizes <- colSums(true.counts)
        cell.facs <- log(lib.sizes) / median(lib.sizes)
        drop.prob <- sapply(1:nCells, function(idx) {
            eta <- cell.facs[idx] * (log(cell.means[, idx]))
            return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
        })

        # Decide which counts to keep
        keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob),
                       nrow = nGenes, ncol = nCells)

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

    scater::counts(sim) <- counts

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
#'
#' @importFrom stats rbinom rlnorm
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
#' @param path.from vector giving the path endpoints that each path orginates
#'        from.
#'
#' @return Vector giving the order to process paths in.
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

#' Brownian bridge
#'
#' Calculate a smoothed Brownian bridge between two points. A Brownian bridge is
#' a random walk with fixed end points.
#'
#' @param x starting value.
#' @param y end value.
#' @param N number of steps in random walk.
#' @param n number of points in smoothed bridge.
#' @param sigma.fac multiplier specifying how extreme each step can be.
#'
#' @return Vector of length n following a path from x to y.
#'
#' @importFrom stats runif rnorm
bridge <- function (x = 0, y = 0, N = 5, n = 100, sigma.fac = 0.8) {

    dt <- 1 / (N - 1)
    t <- seq(0, 1, length = N)
    sigma2 <- runif(1, 0, sigma.fac * mean(c(x, y)))
    X <- c(0, cumsum(rnorm(N - 1, sd = sigma2) * sqrt(dt)))
    BB <- x + X - t * (X[N] - y + x)
    BB <- akima::aspline(BB, n = n)$y
    BB[BB < 0] <- 0

    return(BB)
}
buildBridges <- Vectorize(bridge, vectorize.args = c("x", "y", "sigma.fac"))