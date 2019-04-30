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
#'     \item Set up simulation object
#'     \item Simulate library sizes
#'     \item Simulate gene means
#'     \item Simulate groups/paths
#'     \item Simulate BCV adjusted cell means
#'     \item Simulate true counts
#'     \item Simulate dropout
#'     \item Create final dataset
#' }
#'
#' The final output is a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object that
#' contains the simulated counts but also the values for various intermediate
#' steps. These are stored in the \code{\link{colData}} (for cell specific
#' information), \code{\link{rowData}} (for gene specific information) or
#' \code{\link{assays}} (for gene by cell matrices) slots. This additional
#' information includes:
#' \describe{
#'     \item{\code{colData}}{
#'         \describe{
#'             \item{Cell}{Unique cell identifier.}
#'             \item{Group}{The group or path the cell belongs to.}
#'             \item{ExpLibSize}{The expected library size for that cell.}
#'             \item{Step (paths only)}{how far along the path each cell is.}
#'         }
#'     }
#'     \item{\code{rowData}}{
#'         \describe{
#'             \item{Gene}{Unique gene identifier.}
#'             \item{BaseGeneMean}{The base expression level for that gene.}
#'             \item{OutlierFactor}{Expression outlier factor for that gene.
#'             Values of 1 indicate the gene is not an expression outlier.}
#'             \item{GeneMean}{Expression level after applying outlier factors.}
#'             \item{BatchFac[Batch]}{The batch effects factor for each gene for
#'             a particular batch.}
#'             \item{DEFac[Group]}{The differential expression factor for each
#'             gene in a particular group. Values of 1 indicate the gene is not
#'             differentially expressed.}
#'             \item{SigmaFac[Path]}{Factor applied to genes that have
#'             non-linear changes in expression along a path.}
#'         }
#'     }
#'     \item{\code{assays}}{
#'         \describe{
#'             \item{BatchCellMeans}{The mean expression of genes in each cell
#'             after adding batch effects.}
#'             \item{BaseCellMeans}{The mean expression of genes in each cell
#'             after any differential expression and adjusted for expected
#'             library size.}
#'             \item{BCV}{The Biological Coefficient of Variation for each gene
#'             in each cell.}
#'             \item{CellMeans}{The mean expression level of genes in each cell
#'             adjusted for BCV.}
#'             \item{TrueCounts}{The simulated counts before dropout.}
#'             \item{Dropout}{Logical matrix showing which values have been
#'             dropped in which cells.}
#'         }
#'     }
#' }
#'
#' Values that have been added by Splatter are named using \code{UpperCamelCase}
#' in order to differentiate them from the values added by analysis packages
#' which typically use \code{underscore_naming}.
#'
#' @return SingleCellExperiment object containing the simulated counts and
#' intermediate values.
#'
#' @references
#' Zappia L, Phipson B, Oshlack A. Splatter: simulation of single-cell RNA
#' sequencing data. Genome Biology (2017).
#'
#' Paper: \url{10.1186/s13059-017-1305-0}
#'
#' Code: \url{https://github.com/Oshlack/splatter}
#'
#' @seealso
#' \code{\link{splatSimLibSizes}}, \code{\link{splatSimGeneMeans}},
#' \code{\link{splatSimBatchEffects}}, \code{\link{splatSimBatchCellMeans}},
#' \code{\link{splatSimDE}}, \code{\link{splatSimCellMeans}},
#' \code{\link{splatSimBCVMeans}}, \code{\link{splatSimTrueCounts}},
#' \code{\link{splatSimDropout}}
#'
#' @examples
#' # Simulation with default parameters
#' sim <- splatSimulate()
#' \dontrun{
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
#' }
#' @importFrom SummarizedExperiment rowData colData colData<- assays
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods validObject
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
    nBatches <- getParam(params, "nBatches")
    batch.cells <- getParam(params, "batchCells")
    nGroups <- getParam(params, "nGroups")
    group.prob <- getParam(params, "group.prob")

    if (nGroups == 1 && method == "groups") {
        warning("nGroups is 1, switching to single mode")
        method <- "single"
    }

    if (verbose) {message("Creating simulation object...")}
    # Set up name vectors
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))
    batch.names <- paste0("Batch", seq_len(nBatches))
    if (method == "groups") {
        group.names <- paste0("Group", seq_len(nGroups))
    } else if (method == "paths") {
        group.names <- paste0("Path", seq_len(nGroups))
    }

    # Create SingleCellExperiment to store simulation
    cells <-  data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))

    # Make batches vector which is the index of param$batchCells repeated
    # params$batchCells[index] times
    batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                      b = batch.cells)
    batches <- unlist(batches)
    colData(sim)$Batch <- batch.names[batches]

    if (method != "single") {
        groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                         replace = TRUE)
        colData(sim)$Group <- factor(group.names[groups], levels = group.names)
    }

    if (verbose) {message("Simulating library sizes...")}
    sim <- splatSimLibSizes(sim, params)
    if (verbose) {message("Simulating gene means...")}
    sim <- splatSimGeneMeans(sim, params)
    if (nBatches > 1) {
        if (verbose) {message("Simulating batch effects...")}
        sim <- splatSimBatchEffects(sim, params)
    }
    sim <- splatSimBatchCellMeans(sim, params)
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
    if (verbose) {message("Simulating counts...")}
    sim <- splatSimTrueCounts(sim, params)
    if (verbose) {message("Simulating dropout (if needed)...")}
    sim <- splatSimDropout(sim, params)

    if (verbose) {message("Done!")}
    return(sim)
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
#' Simulate expected library sizes. Typically a log-normal distribution is used
#' but there is also the option to use a normal distribution. In this case any
#' negative values are set to half the minimum non-zero value.
#'
#' @param sim SingleCellExperiment to add library size to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated library sizes.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats rlnorm rnorm
splatSimLibSizes <- function(sim, params) {

    nCells <- getParam(params, "nCells")
    lib.loc <- getParam(params, "lib.loc")
    lib.scale <- getParam(params, "lib.scale")
    lib.norm <- getParam(params, "lib.norm")

    if (lib.norm) {
        exp.lib.sizes <- rnorm(nCells, lib.loc, lib.scale)
        min.lib <- min(exp.lib.sizes[exp.lib.sizes > 0])
        exp.lib.sizes[exp.lib.sizes < 0] <- min.lib / 2
    } else {
        exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    }

    colData(sim)$ExpLibSize <- exp.lib.sizes

    return(sim)
}

#' Simulate gene means
#'
#' Simulate gene means from a gamma distribution. Also simulates outlier
#' expression factors. Genes with an outlier factor not equal to 1 are replaced
#' with the median mean expression multiplied by the outlier factor.
#'
#' @param sim SingleCellExperiment to add gene means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated gene means.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom stats rgamma median
splatSimGeneMeans <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    mean.shape <- getParam(params, "mean.shape")
    mean.rate <- getParam(params, "mean.rate")
    out.prob <- getParam(params, "out.prob")
    out.facLoc <- getParam(params, "out.facLoc")
    out.facScale <- getParam(params, "out.facScale")

    # Simulate base gene means
    base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

    # Add expression outliers
    outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                    out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]

    rowData(sim)$BaseGeneMean <- base.means.gene
    rowData(sim)$OutlierFactor <- outlier.facs
    rowData(sim)$GeneMean <- means.gene

    return(sim)
}

#' Simulate batch effects
#'
#' Simulate batch effects. Batch effect factors for each batch are produced
#' using \code{\link{getLNormFactors}} and these are added along with updated
#' means for each batch.
#'
#' @param sim SingleCellExperiment to add batch effects to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated batch effects.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
splatSimBatchEffects <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nBatches <- getParam(params, "nBatches")
    batch.facLoc <- getParam(params, "batch.facLoc")
    batch.facScale <- getParam(params, "batch.facScale")
    batch.rmEffect <- getParam(params, "batch.rmEffect")
    means.gene <- rowData(sim)$GeneMean

    for (idx in seq_len(nBatches)) {
        batch.facs <- getLNormFactors(nGenes, 1, 0.5, batch.facLoc[idx],
                                        batch.facScale[idx])
        if (batch.rmEffect) batch.facs <- rep(1, length(batch.facs))
        # batch.means.gene <- means.gene * batch.facs
        rowData(sim)[[paste0("BatchFacBatch", idx)]] <- batch.facs
    }

    return(sim)
}

#' Simulate batch means
#'
#' Simulate a mean for each gene in each cell incorporating batch effect
#' factors.
#'
#' @param sim SingleCellExperiment to add batch means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated batch means.
#'
#' @importFrom SummarizedExperiment rowData rowData<- colData
splatSimBatchCellMeans <- function(sim, params) {

    nBatches <- getParam(params, "nBatches")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    gene.means <- rowData(sim)$GeneMean

    if (nBatches > 1) {
        batches <- colData(sim)$Batch
        batch.names <- unique(batches)

        batch.facs.gene <- rowData(sim)[, paste0("BatchFac", batch.names)]
        batch.facs.cell <- as.matrix(batch.facs.gene[,
                                                  as.numeric(factor(batches))])
    } else {
        nCells <- getParam(params, "nCells")
        nGenes <- getParam(params, "nGenes")

        batch.facs.cell <- matrix(1, ncol = nCells, nrow = nGenes)
    }

    batch.means.cell <- batch.facs.cell * gene.means

    colnames(batch.means.cell) <- cell.names
    rownames(batch.means.cell) <- gene.names
    assays(sim)$BatchCellMeans <- batch.means.cell

    return(sim)
}

#' Simulate group differential expression
#'
#' Simulate differential expression. Differential expression factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taken to make sure
#' they are simulated in the correct order.
#'
#' @param sim SingleCellExperiment to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated differential expression.
#'
#' @name splatSimDE
NULL

#' @rdname splatSimDE
#' @importFrom SummarizedExperiment rowData
splatSimGroupDE <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nGroups <- getParam(params, "nGroups")
    de.prob <- getParam(params, "de.prob")
    de.downProb <- getParam(params, "de.downProb")
    de.facLoc <- getParam(params, "de.facLoc")
    de.facScale <- getParam(params, "de.facScale")
    means.gene <- rowData(sim)$GeneMean

    for (idx in seq_len(nGroups)) {
        de.facs <- getLNormFactors(nGenes, de.prob[idx], de.downProb[idx],
                                   de.facLoc[idx], de.facScale[idx])
        group.means.gene <- means.gene * de.facs
        rowData(sim)[[paste0("DEFacGroup", idx)]] <- de.facs
    }

    return(sim)
}

#' @rdname splatSimDE
#' @importFrom SummarizedExperiment rowData
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
            means.gene <- rowData(sim)$GeneMean
        } else {
            means.gene <- rowData(sim)[[paste0("GeneMeanPath", from)]]
        }
        de.facs <- getLNormFactors(nGenes, de.prob[path], de.downProb[path],
                                   de.facLoc[path], de.facScale[path])
        path.means.gene <- means.gene * de.facs
        rowData(sim)[[paste0("DEFacPath", path)]] <- de.facs
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
#' @param sim SingleCellExperiment to add cell means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with added cell means.
#'
#' @name splatSimCellMeans
NULL

#' @rdname splatSimCellMeans
#' @importFrom SummarizedExperiment rowData colData assays assays<-
splatSimSingleCellMeans <- function(sim, params) {

    nCells <- getParam(params, "nCells")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$BatchCellMeans

    cell.means.gene <- batch.means.cell
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' @rdname splatSimCellMeans
#' @importFrom SummarizedExperiment rowData colData assays assays<-
splatSimGroupCellMeans <- function(sim, params) {

    nGroups <- getParam(params, "nGroups")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    groups <- colData(sim)$Group
    group.names <- levels(groups)
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$BatchCellMeans

    group.facs.gene <- rowData(sim)[, paste0("DEFac", group.names)]
    cell.facs.gene <- as.matrix(group.facs.gene[, paste0("DEFac", groups)])
    cell.means.gene <- batch.means.cell * cell.facs.gene
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' @rdname splatSimCellMeans
#' @importFrom SummarizedExperiment rowData colData colData<- assays assays<-
#' @importFrom stats rbinom
splatSimPathCellMeans <- function(sim, params) {

    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    nGroups <- getParam(params, "nGroups")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    path.from <- getParam(params, "path.from")
    path.nSteps <- getParam(params, "path.nSteps")
    path.skew <- getParam(params, "path.skew")
    path.nonlinearProb <- getParam(params, "path.nonlinearProb")
    path.sigmaFac <- getParam(params, "path.sigmaFac")
    groups <- colData(sim)$Group
    exp.lib.sizes <- colData(sim)$ExpLibSize
    batch.means.cell <- assays(sim)$BatchCellMeans

    group.sizes <- table(groups)

    # Generate non-linear path factors
    for (idx in seq_along(path.from)) {
        # Select genes to follow a non-linear path
        is.nonlinear <- as.logical(rbinom(nGenes, 1, path.nonlinearProb))
        sigma.facs <- rep(0, nGenes)
        sigma.facs[is.nonlinear] <- path.sigmaFac
        rowData(sim)[[paste0("SigmaFacPath", idx)]] <- sigma.facs
    }

    # Generate non-linear path factors
    for (idx in seq_along(path.from)) {
        # Select genes to follow a non-linear path
        is.nonlinear <- as.logical(rbinom(nGenes, 1, path.nonlinearProb))
        sigma.facs <- rep(0, nGenes)
        sigma.facs[is.nonlinear] <- path.sigmaFac
        rowData(sim)[[paste0("SigmaFacPath", idx)]] <- sigma.facs
    }

    # Generate paths. Each path is a matrix with path.nSteps columns and
    # nGenes rows where the expression from each genes changes along the path.
    path.steps <- lapply(seq_along(path.from), function(idx) {
        from <- path.from[idx]
        # Find the factors at the starting position
        if (from == 0) {
            facs.start <- rep(1, nGenes)
        } else {
            facs.start <- rowData(sim)[[paste0("DEFacPath", from)]]
        }
        # Find the factors at the end position
        facs.end <- rowData(sim)[[paste0("DEFacPath", idx)]]

        # Get the non-linear factors
        sigma.facs <- rowData(sim)[[paste0("SigmaFacPath", idx)]]

        # Build Brownian bridges from start to end
        steps <- buildBridges(facs.start, facs.end, n = path.nSteps[idx],
                              sigma.fac = sigma.facs)

        return(t(steps))
    })

    # Randomly assign a position in the appropriate path to each cell
    path.probs <- lapply(seq_len(nGroups), function(idx) {
        probs <- seq(path.skew[idx], 1 - path.skew[idx],
                          length = path.nSteps[idx])
        probs <- probs / sum(probs)
        return(probs)
    })

    steps <- sapply(factor(groups), function(path) {
        step <- sample(seq_len(path.nSteps[path]), 1, prob = path.probs[[path]])
    })

    # Collect the underlying expression levels for each cell
    cell.facs.gene <- lapply(seq_len(nCells), function(idx) {
        path <- factor(groups)[idx]
        step <- steps[idx]
        cell.means <- path.steps[[path]][, step]
    })
    cell.facs.gene <- do.call(cbind, cell.facs.gene)

    # Adjust expression based on library size
    cell.means.gene <- batch.means.cell * cell.facs.gene
    cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names

    colData(sim)$Step <- steps
    assays(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' Simulate BCV means
#'
#' Simulate means for each gene in each cell that are adjusted to follow a
#' mean-variance trend using Biological Coefficient of Variation taken from
#' and inverse gamma distribution.
#'
#' @param sim SingleCellExperiment to add BCV means to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated BCV means.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<-
#' @importFrom stats rchisq rgamma
splatSimBCVMeans <- function(sim, params) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    bcv.common <- getParam(params, "bcv.common")
    bcv.df <- getParam(params, "bcv.df")
    base.means.cell <- assays(sim)$BaseCellMeans

    if (is.finite(bcv.df)) {
        bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
            sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    } else {
        warning("'bcv.df' is infinite. This parameter will be ignored.")
        bcv <- (bcv.common + (1 / sqrt(base.means.cell)))
    }

    means.cell <- matrix(rgamma(
        as.numeric(nGenes) * as.numeric(nCells),
        shape = 1 / (bcv ^ 2), scale = base.means.cell * (bcv ^ 2)),
    nrow = nGenes, ncol = nCells)

    colnames(means.cell) <- cell.names
    rownames(means.cell) <- gene.names

    assays(sim)$BCV <- bcv
    assays(sim)$CellMeans <- means.cell

    return(sim)
}

#' Simulate true counts
#'
#' Simulate a true counts matrix. Counts are simulated from a poisson
#' distribution where Each gene in each cell has it's own mean based on the
#' group (or path position), expected library size and BCV.
#'
#' @param sim SingleCellExperiment to add true counts to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated true counts.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<-
#' @importFrom stats rpois
splatSimTrueCounts <- function(sim, params) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    cell.means <- assays(sim)$CellMeans

    true.counts <- matrix(rpois(
        as.numeric(nGenes) * as.numeric(nCells),
        lambda = cell.means),
    nrow = nGenes, ncol = nCells)

    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names

    assays(sim)$TrueCounts <- true.counts

    return(sim)
}

#' Simulate dropout
#'
#' A logistic function is used to form a relationship between the expression
#' level of a gene and the probability of dropout, giving a probability for each
#' gene in each cell. These probabilities are used in a Bernoulli distribution
#' to decide which counts should be dropped.
#'
#' @param sim SingleCellExperiment to add dropout to.
#' @param params SplatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated dropout and observed counts.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<-
#' @importFrom stats rbinom
splatSimDropout <- function(sim, params) {

    dropout.type <- getParam(params, "dropout.type")
    true.counts <- assays(sim)$TrueCounts
    dropout.mid <- getParam(params, "dropout.mid")
    dropout.shape <- getParam(params, "dropout.shape")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    nBatches <- getParam(params, "nBatches")
    nGroups <- getParam(params, "nGroups")
    cell.means <- assays(sim)$CellMeans

    switch(dropout.type,
           experiment = {
               if ((length(dropout.mid) != 1) || length(dropout.shape) != 1) {
                   stop("dropout.type is set to 'experiment' but dropout.mid ",
                        "and dropout.shape aren't length 1")
               }

               dropout.mid <- rep(dropout.mid, nCells)
               dropout.shape <- rep(dropout.shape, nCells)
           },
           batch = {
               if ((length(dropout.mid) != nBatches) ||
                   length(dropout.shape) != nBatches) {
                   stop("dropout.type is set to 'batch' but dropout.mid ",
                        "and dropout.shape aren't length equal to nBatches ",
                        "(", nBatches, ")")
               }

               batches <- as.numeric(factor(colData(sim)$Batch))
               dropout.mid <- dropout.mid[batches]
               dropout.shape <- dropout.shape[batches]
           },
           group = {
               if ((length(dropout.mid) != nGroups) ||
                   length(dropout.shape) != nGroups) {
                   stop("dropout.type is set to 'group' but dropout.mid ",
                        "and dropout.shape aren't length equal to nGroups ",
                        "(", nGroups, ")")
               }

               if ("Group" %in% colnames(colData(sim))) {
                   groups <- as.numeric(colData(sim)$Group)
               } else {
                   stop("dropout.type is set to 'group' but groups have not ",
                        "been simulated")
               }

               dropout.mid <- dropout.mid[groups]
               dropout.shape <- dropout.shape[groups]
           },
           cell = {
               if ((length(dropout.mid) != nCells) ||
                   length(dropout.shape) != nCells) {
                   stop("dropout.type is set to 'cell' but dropout.mid ",
                        "and dropout.shape aren't length equal to nCells ",
                        "(", nCells, ")")
               }
           })

    if (dropout.type != "none") {

        # Generate probabilities based on expression
        drop.prob <- sapply(seq_len(nCells), function(idx) {
            eta <- log(cell.means[, idx])
            return(logistic(eta, x0 = dropout.mid[idx], k = dropout.shape[idx]))
        })

        # Decide which counts to keep
        keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob),
                       nrow = nGenes, ncol = nCells)

        counts <- true.counts * keep

        colnames(drop.prob) <- cell.names
        rownames(drop.prob) <- gene.names
        colnames(keep) <- cell.names
        rownames(keep) <- gene.names

        assays(sim)$DropProb <- drop.prob
        assays(sim)$Dropout <- !keep
    } else {
        counts <- true.counts
    }

    BiocGenerics::counts(sim) <- counts

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
    dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
    factors <- rep(1, n.facs)
    factors[is.selected] <- facs.selected ^ dir.selected

    return(factors)
}

#' Get path order
#'
#' Identify the correct order to process paths so that preceding paths have
#' already been simulated.
#'
#' @param path.from vector giving the path endpoints that each path originates
#'        from.
#'
#' @return Vector giving the order to process paths in.
getPathOrder <- function(path.from) {

    # Transform the vector into a list of (from, to) pairs
    path.pairs <- list()
    for (idx in seq_along(path.from)) {
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
