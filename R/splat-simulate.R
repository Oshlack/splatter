#' Splat simulation
#'
#' Simulate count data from a fictional single-cell RNA-seq experiment using
#' the Splat method.
#'
#' @param params SplatParams object containing parameters for the simulation.
#'        See \code{\link{SplatParams}} for details.
#' @param method which simulation method to use. Options are "single" which
#'        produces a single population, "groups" which produces distinct groups
#'        (eg. cell types) or "paths" which selects cells from continuous
#'        trajectories (eg. differentiation processes).
#' @param verbose logical. Whether to print progress messages.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' Parameters can be set in a variety of ways. If no parameters are provided
#' the default parameters are used. Any parameters in \code{params} can be
#' overridden by supplying additional arguments through a call to
#' \code{\link{setParams}}. This design allows the user flexibility in
#' how they supply parameters and allows small adjustments without creating a
#' new \code{SplatParams} object. See examples for a demonstration of how this
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
#' The final output is an \code{\link[scater]{SCESet}} object that contains the
#' simulated counts but also the values for various intermediate steps. These
#' are stored in the \code{\link[Biobase]{phenoData}} (for cell specific
#' information), \code{\link[Biobase]{featureData}} (for gene specific
#' information) or \code{\link[Biobase]{assayData}} (for gene by cell matrices)
#' slots. This additional information includes:
#' \describe{
#'   \item{\code{phenoData}}{
#'     \describe{
#'       \item{Cell}{Unique cell identifier.}
#'       \item{Group}{The group or path the cell belongs to.}
#'       \item{ExpLibSize}{The expected library size for that cell.}
#'       \item{Step (paths only)}{how far along the path each cell is.}
#'     }
#'   }
#'   \item{\code{featureData}}{
#'     \describe{
#'       \item{Gene}{Unique gene identifier.}
#'       \item{BaseGeneMean}{The base expression level for that gene.}
#'       \item{OutlierFactor}{Expression outlier factor for that gene. Values
#'       of 1 indicate the gene is not an expression outlier.}
#'       \item{GeneMean}{Expression level after applying outlier factors.}
#'       \item{DEFac[Group]}{The differential expression factor for each gene
#'       in a particular group. Values of 1 indicate the gene is not
#'       differentially expressed.}
#'       \item{GeneMean[Group]}{Expression level of a gene in a particular
#'       group after applying differential expression factors.}
#'     }
#'   }
#'   \item{\code{assayData}}{
#'     \describe{
#'       \item{BaseCellMeans}{The expression of genes in each cell adjusted for
#'       expected library size.}
#'       \item{BCV}{The Biological Coefficient of Variation for each gene in
#'       each cell.}
#'       \item{CellMeans}{The expression level of genes in each cell adjusted
#'       for BCV.}
#'       \item{TrueCounts}{The simulated counts before dropout.}
#'       \item{Dropout}{Logical matrix showing which values have been dropped
#'       in which cells.}
#'     }
#'   }
#' }
#'
#' Values that have been added by Splatter are named using \code{CamelCase} in
#' order to differentiate them from the values added by Scater which uses
#' \code{underscore_naming}.
#'
#' @return SCESet object containing the simulated counts and intermediate
#' values.
#'
#' @seealso
#' \code{\link{splatSimLibSizes}}, \code{\link{splatSimGeneMeans}},
#' \code{\link{splatSimDE}}, \code{\link{splatSimCellMeans}},
#' \code{\link{splatSimBCVMeans}}, \code{\link{splatSimTrueCounts}},
#' \code{\link{splatSimDropout}}
#'
#' @examples
#' # Simulation with default parameters
#' sim <- splatSimulate()
#' # Simulation with different number of genes
#' sim <- splatSimulate(nGenes = 1000)
#' # Simulation with custom parameters
#' params <- newSplatParams(nGenes = 100, mean.rate = 0.5)
#' sim <- splatSimulate(params)
#' # Simulation with adjusted custom parameters
#' sim <- splatSimulate(params, mean.rate = 0.6, out.prob = 0.2)
#' # Simulate groups
#' sim <- splatSimulate(method = "groups")
#' # Simulate paths
#' sim <- splatSimulate(method = "paths")
#' @importFrom Biobase fData pData pData<- assayData
#' @importFrom methods validObject
#' @importFrom scater newSCESet counts set_exprs<- get_exprs
#' @export
splatSimulate <- function(params = newSplatParams(),
                          method = c("single", "groups", "paths"),
                          verbose = TRUE, ...) {

    checkmate::assertClass(params, "SplatParams")

    method <- match.arg(method)

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    nGroups <- getParam(params, "nGroups")
    group.cells <- getParam(params, "groupCells")

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
    sim <- splatSimLibSizes(sim, params)
    if (verbose) {message("Simulating gene means...")}
    sim <- splatSimGeneMeans(sim, params)
    if (method == "single") {
        sim <- splatSimSingleCellMeans(sim, params)
    } else if (method == "groups") {
        if (verbose) {message("Simulating group DE...")}
        sim <- splatSimGroupDE(sim, params)
        if (verbose) {message("Simulating cell means...")}
        sim <- splatSimGroupCellMeans(sim, params)
    } else {
        if (verbose) {message("Simulating path endpoints...")}
        sim <- splatSimPathDE(sim, params)
        if (verbose) {message("Simulating path steps...")}
        sim <- splatSimPathCellMeans(sim, params)
    }
    if (verbose) {message("Simulating BCV...")}
    sim <- splatSimBCVMeans(sim, params)
    if (verbose) {message("Simulating counts..")}
    sim <- splatSimTrueCounts(sim, params)
    if (verbose) {message("Simulating dropout...")}
    sim <- splatSimDropout(sim, params)

    if (verbose) {message("Creating final SCESet...")}
    # Create new SCESet to make sure values are calculated correctly
    yo <- counts(sim)
    rownames(yo) <- featureNames(sim)
    colnames(yo) <- sampleNames(sim)
    sce <- newSCESet(countData = yo,
                     phenoData = new("AnnotatedDataFrame", data = pData(sim)),
                     featureData = new("AnnotatedDataFrame", data = fData(sim)))

    # Add intermediate matrices stored in assayData
    for (assay.name in names(assayData(sim))) {
        if (!(assay.name %in% names(assayData(sce)))) {
            set_exprs(sce, assay.name) <- get_exprs(sim, assay.name)
        }
    }

    if (verbose) {message("Done!")}
    return(sce)
}

#' @rdname splatSimulate
#' @export
splatSimulateSingle <- function(params = newSplatParams(),
                                verbose = TRUE, ...) {
    sim <- splatSimulate(params = params, method = "single",
                         verbose = verbose, ...)
    return(sim)
}

#' @rdname splatSimulate
#' @export
splatSimulateGroups <- function(params = newSplatParams(),
                                verbose = TRUE, ...) {
    sim <- splatSimulate(params = params, method = "groups",
                         verbose = verbose, ...)
    return(sim)
}

#' @rdname splatSimulate
#' @export
splatSimulatePaths <- function(params = newSplatParams(), verbose = TRUE, ...) {
    sim <- splatSimulate(params = params, method = "paths",
                         verbose = verbose, ...)
    return(sim)
}

#' Simulate library sizes
#'
#' Simulate expected library sizes from a log-normal distribution
#'
#' @param sim SCESet to add library size to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SCESet with simulated library sizes.
#'
#' @importFrom Biobase pData pData<-
#' @importFrom stats rlnorm
splatSimLibSizes <- function(sim, params) {

    nCells <- getParam(params, "nCells")
    lib.loc <- getParam(params, "lib.loc")
    lib.scale <- getParam(params, "lib.scale")

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
#' @param params SplatParams object with simulation parameters.
#'
#' @return SCESet with simulated gene means.
#'
#' @importFrom Biobase fData fData<-
#' @importFrom stats rgamma median
splatSimGeneMeans <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    out.prob <- getParam(params, "out.prob")
    out.loProb <- getParam(params, "out.loProb")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")

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
#' Simulate differential expression. Differential expression factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taked to make sure
#' they are simualated in the correct order.
#'
#' @param sim SCESet to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SCESet with simulated differential expression.
#'
#' @name splatSimDE
NULL

#' @rdname splatSimDE
#' @importFrom Biobase fData
splatSimGroupDE <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nGroups <- getParam(params, "nGroups")
    de.prob <- getParam(params, "de.prob")
    de.downProb <- getParam(params, "de.downProb")
    de.facLoc <- getParam(params, "de.facLoc")
    de.facScale <- getParam(params, "de.facScale")
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

#' @rdname splatSimDE
#' @importFrom Biobase fData
splatSimPathDE <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    de.prob <- getParam(params, "de.prob")
    de.downProb <- getParam(params, "de.downProb")
    de.facLoc <- getParam(params, "de.facLoc")
    de.facScale <- getParam(params, "de.facScale")
    path.from <- getParam(params, "path.from")

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

#' Simulate cell means
#'
#' Simulate a gene by cell matrix giving the mean expression for each gene in
#' each cell. Cells start with the mean expression for the group they belong to
#' (when simulating groups) or cells are assigned the mean expression from a
#' random position on the appropriate path (when simulating paths). The selected
#' means are adjusted for each cell's expected library size.
#'
#' @param sim SCESet to add cell means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SCESet with added cell means.
#'
#' @name splatSimCellMeans
NULL

#' @rdname splatSimCellMeans
#' @importFrom Biobase fData pData
#' @importFrom scater set_exprs<-
splatSimSingleCellMeans <- function(sim, params) {

    nCells <- getParam(params, "nCells")
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    exp.lib.sizes <- pData(sim)$ExpLibSize

    cell.means.gene <- as.matrix(fData(sim)[, rep("GeneMean", nCells)])
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    set_exprs(sim, "BaseCellMeans") <- base.means.cell

    return(sim)
}

#' @rdname splatSimCellMeans
#' @importFrom Biobase fData pData
#' @importFrom scater set_exprs<-
splatSimGroupCellMeans <- function(sim, params) {

    nGroups <- getParam(params, "nGroups")
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
    set_exprs(sim, "BaseCellMeans") <- base.means.cell

    return(sim)
}

#' @rdname splatSimCellMeans
#' @importFrom Biobase fData pData
#' @importFrom scater set_exprs<-
#' @importFrom stats rbinom
splatSimPathCellMeans <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nGroups <- getParam(params, "nGroups")
    group.cells <- getParam(params, "groupCells")
    path.from <- getParam(params, "path.from")
    path.length <- getParam(params, "path.length")
    path.skew <- getParam(params, "path.skew")
    path.nonlinearProb <- getParam(params, "path.nonlinearProb")
    path.sigmaFac <- getParam(params, "path.sigmaFac")
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

    pData(sim)$Step <- unlist(cell.steps)
    set_exprs(sim, "BaseCellMeans") <- base.means.cell

    return(sim)
}

#' Simulate BCV means
#'
#' Simulate means for each gene in each cell that are adjusted to follow a
#' mean-variance trend using Biological Coefficient of Variation taken from
#' and inverse gamma distribution.
#'
#' @param sim SCESet to add BCV means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SCESet with simulated BCV means.
#'
#' @importFrom Biobase fData pData
#' @importFrom scater get_exprs set_exprs<-
#' @importFrom stats rchisq rgamma
splatSimBCVMeans <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    bcv.common <- getParam(params, "bcv.common")
    bcv.df <- getParam(params, "bcv.df")
    base.means.cell <- get_exprs(sim, "BaseCellMeans")

    bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
        sqrt(bcv.df / rchisq(nGenes, df = bcv.df))

    means.cell <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
                                scale = base.means.cell * (bcv ^ 2)),
                         nrow = nGenes, ncol = nCells)

    set_exprs(sim, "BCV") <- bcv
    set_exprs(sim, "CellMeans") <- means.cell

    return(sim)
}

#' Simulate true counts
#'
#' Simulate a true counts matrix. Counts are simulated from a poisson
#' distribution where Each gene in each cell has it's own mean based on the
#' group (or path position), expected library size and BCV.
#'
#' @param sim SCESet to add true counts to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SCESet with simulated true counts.
#'
#' @importFrom Biobase fData pData
#' @importFrom scater get_exprs set_exprs<-
#' @importFrom stats rpois
splatSimTrueCounts <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    cell.means <- get_exprs(sim, "CellMeans")

    true.counts <- matrix(rpois(nGenes * nCells, lambda = cell.means),
                          nrow = nGenes, ncol = nCells)

    set_exprs(sim, "TrueCounts") <- true.counts

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
#' @param params SplatParams object with simulation parameters.
#'
#' @return SCESet with simulated dropout and observed counts.
#'
#' @importFrom Biobase fData pData
#' @importFrom scater get_exprs set_exprs<-
#' @importFrom stats rbinom
splatSimDropout <- function(sim, params) {

    dropout.present <- getParam(params, "dropout.present")
    true.counts <- get_exprs(sim, "TrueCounts")

    if (dropout.present) {
        nCells <- getParam(params, "nCells")
        nGenes <- getParam(params, "nGenes")
        dropout.mid <- getParam(params, "dropout.mid")
        dropout.shape <- getParam(params, "dropout.shape")
        cell.means <- get_exprs(sim, "CellMeans")

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

        set_exprs(sim, "DropProb") <- drop.prob
        set_exprs(sim, "Dropout") <- !keep
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
    BB[BB < 0] <- 1e-6

    return(BB)
}
buildBridges <- Vectorize(bridge, vectorize.args = c("x", "y", "sigma.fac"))