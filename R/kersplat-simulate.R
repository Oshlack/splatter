#' Kersplat simulation
#'
#' Simulate scRNA-seq count data using the Kersplat model
#'
#' @param params KersplatParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#'
#' This functions is for simulating data in a single step. It consists of a
#' call to \code{\link{kersplatSetup}} followed by a call to
#' \code{\link{kersplatSample}}. Please see the documentation for those
#' functions for more details of the individual steps.
#'
#' @seealso
#' \code{\link{kersplatSetup}}, \code{\link{kersplatSample}}
#'
#' @return SingleCellExperiment containing simulated counts and intermediate
#' values
#'
#' @examples
#'
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'     sim <- kersplatSimulate
#' }
#'
#' @export
kersplatSimulate <- function(params = newKersplatParams(), sparsify = TRUE,
                             verbose = TRUE, ...) {

    params <- kersplatSetup(params, verbose, ...)
    sim <- kersplatSample(params, sparsify, verbose)

    return(sim)
}

#' Kersplat setup
#'
#' Setup the parameters required for the Kersplat simulation
#'
#' @param params KersplatParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' The first stage is a two-step Kersplat simulation is to generate some of the
#' intermediate parameters. The resulting parameters allow multiple simulated
#' datasets to be generated from the same biological structure (using
#' \code{\link{kersplatSample}}). As with all the other parameters these values
#' can be manually overwritten if desired.
#'
#' The setup involves the following steps:
#' \enumerate{
#'     \item Generate a gene network (if not already present)
#'     \item Select regulator genes (if not already present)
#'     \item Simulate gene means (if not already present)
#'     \item Simulate cell paths
#' }
#'
#' The resulting \code{\link{KersplatParams}} object will have the following
#' parameters set (if they weren't already).
#'
#' \itemize{
#'     \item \code{mean.values}
#'     \item \code{network.graph}
#'     \item \code{network.regsSet}
#'     \item \code{paths.means}
#' }
#'
#' See \code{\link{KersplatParams}} for more details about these parameters and
#' the functions for the individual steps for more details about the process.
#'
#' @seealso
#' \code{\link{kersplatGenNetwork}}, \code{\link{kersplatSelectRegs}},
#' \code{\link{kersplatSimGeneMeans}}, \code{\link{kersplatSimPaths}},
#' \code{\link{KersplatParams}}
#'
#' @return A complete KersplatParams object
#' @export
#'
#' @examples
#'
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'     params <- kersplatSetup()
#' }
kersplatSetup <- function(params = newKersplatParams(), verbose = TRUE, ...) {

    checkmate::assertClass(params, "KersplatParams")
    params <- setParams(params, ...)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    if (verbose) {message("Setting up parameters...")}
    params <- kersplatGenNetwork(params, verbose)
    params <- kersplatSelectRegs(params, verbose)
    params <- kersplatSimGeneMeans(params, verbose)
    params <- kersplatSimPaths(params, verbose)

    return(params)
}

#' Kersplat sample
#'
#' Sample cells for the Kersplat simulation
#'
#' @param params KersplatParams object containing simulation parameters.
#' @param sparsify logical. Whether to automatically convert assays to sparse
#'        matrices if there will be a size reduction.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' The second stage is a two-step Kersplat simulation is to generate cells based
#' on a complete \code{\link{KersplatParams}} object.
#' intermediate parameters.
#'
#' The sampling process involves the following steps:
#' \enumerate{
#'     \item Simulate library sizes for each cell
#'     \item Simulate means for each cell
#'     \item Simulate endogenous counts for each cell
#'     \item Simulate ambient counts for each cell
#'     \item Simulate final counts for each cell
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
#'             \item{Type}{Whether the cell is a Cell, Doublet or Empty.}
#'             \item{CellLibSize}{The expected number of endogenous counts for
#'             that cell.}
#'             \item{AmbientLibSize}{The expected number of ambient counts for
#'             that cell.}
#'             \item{Path}{The path the cell belongs to.}
#'             \item{Step}{How far along the path each cell is.}
#'             \item{Path1}{For doublets the path of the first partner in the
#'             doublet (otherwise \code{NA}).}
#'             \item{Step1}{For doublets the step of the first partner in the
#'             doublet (otherwise \code{NA}).}
#'             \item{Path2}{For doublets the path of the second partner in the
#'             doublet (otherwise \code{NA}).}
#'             \item{Step2}{For doublets the step of the second partner in the
#'             doublet (otherwise \code{NA}).}
#'         }
#'     }
#'     \item{\code{rowData}}{
#'         \describe{
#'             \item{Gene}{Unique gene identifier.}
#'             \item{BaseMean}{The base expression level for that gene.}
#'             \item{AmbientMean}{The ambient expression level for that gene.}
#'         }
#'     }
#'     \item{\code{assays}}{
#'         \describe{
#'             \item{CellMeans}{The mean expression of genes in each cell
#'             after any differential expression and adjusted for expected
#'             library size.}
#'             \item{CellCounts}{Endogenous count matrix.}
#'             \item{AmbientCounts}{Ambient count matrix.}
#'             \item{counts}{Final count matrix.}
#'         }
#'     }
#' }
#'
#' Values that have been added by Splatter are named using \code{UpperCamelCase}
#' in order to differentiate them from the values added by analysis packages
#' which typically use \code{underscore_naming}.
#'
#' @seealso
#' \code{\link{kersplatSimLibSizes}}, \code{\link{kersplatSimCellMeans}},
#' \code{\link{kersplatSimCellCounts}}, \code{\link{kersplatSimAmbientCounts}},
#' \code{\link{kersplatSimCounts}}
#'
#' @return SingleCellExperiment object containing the simulated counts and
#' intermediate values.
#' @export
#'
#' @examples
#'
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'     params <- kersplatSetup()
#'     sim <- kersplatSample(params)
#' }
kersplatSample <- function(params, sparsify = TRUE, verbose = TRUE) {

    # Check that parameters are set up
    checkmate::assertClass(params, "KersplatParams")
    network.graph <- getParam(params, "network.graph")
    if (is.null(network.graph)) {
        stop("'network.graph' not set, run kersplatSetup first")
    }
    network.regsSet <- getParam(params, "network.regsSet")
    if (!network.regsSet) {
        stop("network regulators not set, run kersplatSetup first")
    }
    mean.values <- getParam(params, "mean.values")
    if (length(mean.values) == 0) {
        stop("'mean.values' not set, run kersplatSetup first")
    }
    paths.means <- getParam(params, "paths.means")
    if (length(mean.values) == 0) {
        stop("'paths.means' not set, run kersplatSetup first")
    }

    if (verbose) {message("Creating simulation object...")}
    nGenes <- getParam(params, "nGenes")
    gene.names <- paste0("Gene", seq_len(nGenes))

    nCells <- getParam(params, "nCells")
    doublet.prop <- getParam(params, "doublet.prop")
    nDoublets <- floor(nCells * doublet.prop)
    if (doublet.prop > 0) {
        nCells <- nCells - nDoublets
        cell.names <- c(paste0("Cell", seq_len(nCells)),
                        paste0("Doublet", seq_len(nDoublets)))
    } else {
        cell.names <- paste0("Cell", seq_len(nCells))
    }

    nEmpty <- getParam(params, "ambient.nEmpty")
    if (nEmpty > 0) {
        empty.names <- paste0("Empty", seq_len(nEmpty))
        cell.names <- c(cell.names, empty.names)
    }

    cells <-  data.frame(Cell = cell.names,
                         Type = rep(c("Cell", "Doublet", "Empty"),
                                    c(nCells, nDoublets, nEmpty)),
                         row.names = cell.names)
    features <- data.frame(Gene = gene.names,
                           BaseMean = getParam(params, "mean.values"),
                           row.names = gene.names)
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))

    sim <- kersplatSimLibSizes(sim, params, verbose)
    sim <- kersplatSimCellMeans(sim, params, verbose)
    sim <- kersplatSimCellCounts(sim, params, verbose)
    sim <- kersplatSimAmbientCounts(sim, params, verbose)
    sim <- kersplatSimCounts(sim, params, verbose)

    if (sparsify) {
        if (verbose) {message("Sparsifying assays...")}
        assays(sim) <- sparsifyMatrices(assays(sim), auto = TRUE,
                                        verbose = verbose)
    }

    return(sim)

}

#' Generate Kersplat gene network
#'
#' Generate a gene network for the Kersplat simulation
#'
#' @param params KersplatParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' Currently a very simple approach is used which needs to be improved. A
#' network is generated using the \code{\link[igraph]{sample_forestfire}}
#' function and edge weights are sampled from a standard normal distribution.
#'
#' @return KersplatParams object with gene network
kersplatGenNetwork <- function(params, verbose) {

    nGenes <- getParam(params, "nGenes")
    network.graph <- getParam(params, "network.graph")

    if (!is.null(network.graph)) {
        if (verbose) {message("Using provided gene network...")}
        return(params)
    }

    if (verbose) {message("Generating gene network...")}

    graph.raw <- igraph::sample_forestfire(nGenes, 0.1)
    graph.data <- igraph::get.data.frame(graph.raw)
    graph.data <- graph.data[, c("from", "to")]
    graph.data$weight <- rnorm(nrow(graph.data))
    graph <- igraph::graph.data.frame(graph.data)

    params <- setParam(params, "network.graph", graph)

    return(params)
}

#' Select Kersplat regulators
#'
#' Select regulator genes in the gene network for a Kersplat simulation
#'
#' @param params KersplatParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' Regulators are randomly selected, weighted according to the difference
#' between their out degree and in degree. This is an arbitrary weighting and
#' may be improved or replace in the future.
#'
#' @return KersplatParams object with gene regulators
kersplatSelectRegs <- function(params, verbose) {

    network.regsSet <- getParam(params, "network.regsSet")

    if (network.regsSet) {
        if (verbose) {message("Using selected regulators...")}
        return(params)
    }

    if (verbose) {message("Selecting regulators...")}
    network.nRegs <- getParam(params, "network.nRegs")
    network.graph <- getParam(params, "network.graph")

    out.degree <- igraph::degree(network.graph, mode = "out")
    in.degree <- igraph::degree(network.graph, mode = "in")
    reg.prob <-  out.degree - in.degree
    reg.prob <- reg.prob + rnorm(length(reg.prob))
    reg.prob[reg.prob <= 0] <- 1e-10
    reg.prob <- reg.prob / sum(reg.prob)
    reg.nodes <- names(rev(sort(reg.prob))[seq_len(network.nRegs)])
    is.reg <- igraph::V(network.graph)$name %in% reg.nodes
    network.graph <- igraph::set_vertex_attr(network.graph, "IsReg",
                                             value = is.reg)

    params <- setParam(params, "network.graph", network.graph)

    return(params)
}

#' Simulate Kersplat gene means
#'
#' @param params KersplatParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' Gene means are simulated in one of two ways depending on the value of the
#' \code{mean.method} parameter.
#'
#' If \code{mean.method} is "fit" (default) then means are sampled from a Gamma
#' distribution with shape equals \code{mean.shape} and rate equals
#' \code{mean.rate}. Expression outliers are then added by replacing some
#' values with the median multiplied by a factor from a log-normal distribution.
#' This is the same process used for the Splat simulation.
#'
#' If \code{mean.method} is "density" then means are sampled from the
#' density object in the \code{mean.density} parameter using a rejection
#' sampling method. This approach is more flexible but may violate some
#' statistical assumptions.
#'
#' @return KersplatParams object with gene means
kersplatSimGeneMeans <- function(params, verbose) {

    mean.values <- getParam(params, "mean.values")

    # Generate means
    if (length(mean.values) != 0) {
        if (verbose) {message("Using defined means...")}
        return(params)
    }

    if (verbose) {message("Simulating means...")}
    nGenes <- getParam(params, "nGenes")
    mean.method <- getParam(params, "mean.method")

    if (mean.method == "fit") {
        if (verbose) {message("Sampling from gamma distribution...")}
        mean.shape <- getParam(params, "mean.shape")
        mean.rate <- getParam(params, "mean.rate")
        mean.outProb <- getParam(params, "mean.outProb")
        mean.outLoc <- getParam(params, "mean.outLoc")
        mean.outScale <- getParam(params, "mean.outScale")

        mean.values <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

        outlier.facs <- getLNormFactors(nGenes, mean.outProb, 0, mean.outLoc,
                                        mean.outScale)
        median.means.gene <- median(mean.values)
        outlier.means <- median.means.gene * outlier.facs
        is.outlier <- outlier.facs != 1
        mean.values[is.outlier] <- outlier.means[is.outlier]
    } else if (mean.method == "density") {
        if (verbose) {message("Sampling from density object...")}
        mean.dens <- getParam(params, "mean.dens")

        mean.values <- exp(sampleDensity(nGenes, mean.dens, lower = -Inf))
    }

    params <- setParam(params, "mean.values", mean.values)

    return(params)
}

#' Simulate Kersplat paths
#'
#' Simulate gene means for each step along each path of a Kersplat simulation
#'
#' @param params KersplatParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' The method of simulating paths is inspired by the method used in the PROSSTT
#' simulation. Changes in expression are controlled by \code{paths.nPrograms}
#' regulatory programs. Each of the regulatory genes in the gene network has
#' some association with each program. This is analogous to there being changes
#' in the environment (the programs) which are sensed by receptors (regulatory
#' genes) and cause changes in expression downstream. For each path a random
#' walk is generated for each program and the changes passed on to the
#' regulatory genes. At each step the changes propagate through the network
#' according to the weights on edges between genes. This algorithm is fairly
#' simple but should result in correlation relationships between genes. However
#' it is likely to be improved and adjusted in the future.
#'
#' The path structure itself is specified by the \code{paths.design} parameter.
#' This is a \code{data.frame} with three columns: "Path", "From", and "Steps".
#' The Path field is an ID for each path while the Steps field controls the
#' length of each path. Increasing the number of steps will increase the
#' difference in expression between the ends of the paths. The From field sets
#' the originating point of each path. For example a From of \code{0, 0, 0}
#' would indicate three paths from the origin while a From of \code{0, 1, 1}
#' would give a branching structure with Path 1 beginning at the origin and
#' Path 2 and Path 3 beginning at the end of Path 1.
#'
#' @references
#'
#' Papadopoulos N, Parra RG, Söding J. PROSSTT: probabilistic simulation of
#' single-cell RNA-seq data for complex differentiation processes.
#' Bioinformatics (2019). \url{https://doi.org/10.1093/bioinformatics/btz078}.
#'
#' @return KersplatParams object with path means
kersplatSimPaths <- function(params, verbose) {

    paths.means <- getParam(params, "paths.means")

    if (length(paths.means) != 0) {
        if (verbose) {message("Using defined path means...")}
        return(params)
    }

    if (verbose) {message("Simulating paths...")}
    nGenes <- getParam(params, "nGenes")
    paths.design <- getParam(params, "paths.design")
    network.graph <- getParam(params, "network.graph")
    network.weights <- igraph::as_adjacency_matrix(network.graph,
                                                   attr = "weight")
    network.nRegs <- getParam(params, "network.nRegs")
    network.isReg <- igraph::vertex_attr(network.graph, "IsReg")
    paths.nPrograms <- getParam(params, "paths.nPrograms")

    programs.weights <- matrix(rnorm(network.nRegs * paths.nPrograms),
                               nrow = network.nRegs, ncol = paths.nPrograms)
    paths.changes <- vector("list", nrow(paths.design))
    paths.factors <- vector("list", nrow(paths.design))

    paths.graph <- igraph::graph_from_data_frame(paths.design)
    paths.order <- names(igraph::topo_sort(paths.graph, mode = "in"))
    paths.order <- as.numeric(paths.order)
    # Remove the origin because it is not a path
    paths.order <- paths.order[paths.order != 0]

    for (path in paths.order) {
        if (verbose) {message("Simulating path ", path, "...")}
        nSteps <- paths.design$Steps[path]
        from <- paths.design$From[path]
        changes <- matrix(0, nrow = nGenes, ncol = nSteps + 1)

        if (from != 0) {
            from.changes <- paths.changes[[from]]
            changes[, 1] <- from.changes[, ncol(from.changes)]
        }

        for (step in seq_len(nSteps) + 1) {
            programs.changes <- rnorm(paths.nPrograms, sd = 0.01)
            reg.changes <- as.vector(programs.weights %*% programs.changes)
            changes[network.isReg, step] <- reg.changes
            change <- as.vector(changes[, step - 1] %*% network.weights)
            changes[, step] <- changes[, step] + change
        }

        if (from == 0) {
            changes <- changes[, seq_len(nSteps)]
            factors <- matrixStats::rowCumsums(changes)
        } else {
            changes <- changes[, seq_len(nSteps) + 1]
            from.factors <- paths.factors[[from]][, ncol(paths.factors[[from]])]
            factors <- matrixStats::rowCumsums(changes) + from.factors
        }
        paths.changes[[path]] <- changes
        paths.factors[[path]] <- factors
    }

    mean.values <- getParam(params, "mean.values")
    paths.means <- lapply(paths.factors, function(x) {
        (2 ^ x) * mean.values
    })

    names(paths.means) <- paste0("Path", paths.design$Path)
    params <- setParam(params, "paths.means", paths.means)

    return(params)
}

#' Simulate Kersplat library sizes
#'
#' Generate library sizes for cells in the Kersplat simulation
#'
#' @param sim SingleCellExperiment containing simulation.
#' @param params KersplatParams object with simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' Library sizes are simulated in one of two ways depending on the value of the
#' \code{lib.method} parameter.
#'
#' If \code{lib.method} is "fit" (default) then means are sampled from a
#' log-normal distribution with meanlog equals \code{lib.loc} and sdlog equals
#' \code{lib.scale}.
#'
#' If \code{mean.method} is "density" then library sizes are sampled from the
#' density object in the \code{lib.density} parameter using a rejection
#' sampling method. This approach is more flexible but may violate some
#' statistical assumptions.
#'
#' Ambient library sizes are also generated from a log-normal distribution based
#' on the parameters for the cell library size and adjusted using the
#' \code{ambient.scale} parameter.
#'
#' @return SingleCellExperiment with library sizes
kersplatSimLibSizes <- function(sim, params, verbose) {

    if (verbose) {message("Simulating library sizes...")}
    nCells <- getParam(params, "nCells")
    nEmpty <- getParam(params, "ambient.nEmpty")
    is.doublet <- colData(sim)$Type == "Doublet"
    lib.method <- getParam(params, "lib.method")

    if (lib.method == "fit") {
        if (verbose) {message("Sampling from log-normal distribution...")}
        lib.loc <- getParam(params, "lib.loc")
        lib.scale <- getParam(params, "lib.scale")

        cell.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    } else if (lib.method == "density") {
        if (verbose) {message("Sampling from density object...")}
        lib.dens <- getParam(params, "lib.dens")

        cell.lib.sizes <- sampleDensity(nCells, lib.dens)
    }

    cell.lib.sizes <- c(cell.lib.sizes, rep(0, nEmpty))
    cell.lib.sizes[is.doublet] <- 1.5 * cell.lib.sizes[is.doublet]
    colData(sim)$CellLibSize <- cell.lib.sizes

    ambient.scale <- getParam(params, "ambient.scale")
    if (ambient.scale > 0) {
        ambient.loc <- log(exp(lib.loc) * ambient.scale)

        ambient.lib.sizes <- rlnorm(nCells + nEmpty, ambient.loc, 0.3)
        colData(sim)$AmbientLibSize <- ambient.lib.sizes
    }

    return(sim)
}

#' Simulate Kersplat cell means
#'
#' Simulate endogenous counts for each cell in a Kersplat simulation
#'
#' @param sim SingleCellExperiment containing simulation.
#' @param params KersplatParams object with simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' Cells are first assigned to a path and a step along that path. This is
#' controlled by the \code{cells.design} parameter which is a \code{data.frame}
#' with the columns "Path", "Probability", "Alpha" and "Beta". The Path field
#' is an ID for each path and the Probability field is the probability that a
#' cell will come from that path (must sum to 1). The Alpha and Beta parameters
#' control the density of cells along the path. After they are assigned to paths
#' the step for each cell is sampled from a Beta distribution with parameters
#' shape1 equals Alpha and shape2 equals beta. This approach is very flexible
#' and allows almost any distribution of cells along a path. The distribution
#' can be viewed using \code{hist(rbeta(10000, Alpha, Beta), breaks = 100)}.
#' Some useful combinations of parameters are:
#'
#' \describe{
#'     \item{\code{Alpha = 1}, \code{Beta = 1}}{Uniform distribution along the
#'     path}
#'     \item{\code{Alpha = 0}, \code{Beta = 1}}{All cells at the start of the
#'     path.}
#'     \item{\code{Alpha = 1}, \code{Beta = 0}}{All cells at the end of the
#'     path.}
#'     \item{\code{Alpha = 0}, \code{Beta = 0}}{Cells only at each end of the
#'     path.}
#'     \item{\code{Alpha = 1}, \code{Beta = 2}}{Linear skew towards the start
#'     of the path}
#'     \item{\code{Alpha = 0.5}, \code{Beta = 1}}{Curved skew towards the start
#'     of the path}
#'     \item{\code{Alpha = 2}, \code{Beta = 1}}{Linear skew towards the end
#'     of the path}
#'     \item{\code{Alpha = 1}, \code{Beta = 0.5}}{Curved skew towards the end
#'     of the path}
#'     \item{\code{Alpha = 0.5}, \code{Beta = 0.5}}{Curved skew towards both
#'     ends of the path}
#'     \item{\code{Alpha = 0.5}, \code{Beta = 0.5}}{Curved skew away from both
#'     ends of the path}
#' }
#'
#' Once cells are assigned to paths and steps the correct means are extracted
#' from the \code{paths.means} parameter and adjusted based on each cell's
#' library size. An adjustment for BCV is then applied. Doublets are also
#' simulated at this stage by selecting two path/step combinations and averaging
#' the means.
#'
#' @return SingleCellExperiment with cell means
kersplatSimCellMeans <- function(sim, params, verbose) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nDoublets <- sum(colData(sim)$Type == "Doublet")
    nCells <- getParam(params, "nCells") - nDoublets
    cells.design <- getParam(params, "cells.design")
    paths.design <- getParam(params, "paths.design")
    paths.means <- getParam(params, "paths.means")
    cell.lib.sizes <- colData(sim)$CellLibSize
    nEmpty <- getParam(params, "ambient.nEmpty")
    not.empty <- colData(sim)$Type != "Empty"

    if (verbose) {message("Assigning cells to paths...")}
    cells.paths <- sample(cells.design$Path, nCells, replace = TRUE,
                          prob = cells.design$Probability)

    if (verbose) {message("Assigning cells to steps...")}
    paths.cells.design <- merge(paths.design, cells.design)
    steps.probs <- apply(paths.cells.design, 1, function(path) {
        steps <- path["Steps"]
        probs <- getBetaStepProbs(path["Steps"], path["Alpha"], path["Beta"])

        # Return a list to avoid getting a matrix if all path lengths are equal
        return(list(probs))
    })
    # Remove unnecessary list level
    steps.probs <- lapply(steps.probs, "[[", 1)
    names(steps.probs) <- paths.cells.design$Path

    cells.steps <- vapply(cells.paths, function(path) {
        probs <- steps.probs[[path]]
        step <- sample(seq_len(length(probs)), 1, prob = probs)
        return(step)
    }, c(Step = 0))

    if (verbose) {message("Simulating cell means...")}
    cells.means <- vapply(seq_len(nCells), function(cell) {
        path <- cells.paths[cell]
        step <- cells.steps[cell]
        means <- paths.means[[path]][, step]
        return(means)
    }, as.numeric(seq_len(nGenes)))

    if (nDoublets > 0) {
        if (verbose) {message("Assigning doublets...")}
        doublet.paths1 <- sample(cells.design$Path, nDoublets, replace = TRUE,
                                 prob = cells.design$Probability)
        doublet.paths2 <- sample(cells.design$Path, nDoublets, replace = TRUE,
                                 prob = cells.design$Probability)

        doublet.steps1 <- vapply(doublet.paths1, function(path) {
            probs <- steps.probs[[path]]
            step <- sample(seq_len(length(probs)), 1, prob = probs)
            return(step)
        }, c(Step1 = 0))
        doublet.steps2 <- vapply(doublet.paths2, function(path) {
            probs <- steps.probs[[path]]
            step <- sample(seq_len(length(probs)), 1, prob = probs)
            return(step)
        }, c(Step2 = 0))

        if (verbose) {message("Simulating doublet means...")}
        doublet.means1 <- vapply(seq_len(nDoublets), function(doublet) {
            path <- doublet.paths1[doublet]
            step <- doublet.steps1[doublet]
            means <- paths.means[[path]][, step]
            return(means)
        }, as.numeric(seq_len(nGenes)))
        doublet.means2 <- vapply(seq_len(nDoublets), function(doublet) {
            path <- doublet.paths2[doublet]
            step <- doublet.steps2[doublet]
            means <- paths.means[[path]][, step]
            return(means)
        }, as.numeric(seq_len(nGenes)))
        doublet.means <- (doublet.means1 + doublet.means2) * 0.5

        cells.means <- cbind(cells.means, doublet.means)
    }

    # Adjust mean based on library size
    cells.props <- t(t(cells.means) / colSums(cells.means))
    cells.means <- t(t(cells.props) * cell.lib.sizes[not.empty])

    if (verbose) {message("Applying BCV adjustment...")}
    nGenes <- getParam(params, "nGenes")
    bcv.common <- getParam(params, "bcv.common")
    bcv.df <- getParam(params, "bcv.df")

    if (is.finite(bcv.df)) {
        bcv <- (bcv.common + (1 / sqrt(cells.means))) *
            sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    } else {
        warning("'bcv.df' is infinite. This parameter will be ignored.")
        bcv <- (bcv.common + (1 / sqrt(cells.means)))
    }

    cells.means <- matrix(rgamma(
        as.numeric(nGenes) * as.numeric(nCells + nDoublets),
        shape = 1 / (bcv ^ 2), scale = cells.means * (bcv ^ 2)),
        nrow = nGenes, ncol = nCells + nDoublets)

    empty.means <- matrix(0, nrow = nGenes, ncol = nEmpty)
    cells.means <- cbind(cells.means, empty.means)

    colnames(cells.means) <- cell.names
    rownames(cells.means) <- gene.names

    colData(sim)$Path <- factor(c(cells.paths, rep(NA, nDoublets),
                                  rep(NA, nEmpty)))
    colData(sim)$Step <- c(cells.steps, rep(NA, nDoublets), rep(NA, nEmpty))

    if (nDoublets > 0) {
        colData(sim)$Path1 <- factor(c(rep(NA, nCells), doublet.paths1,
                                       rep(NA, nEmpty)))
        colData(sim)$Step1 <- c(rep(NA, nCells), doublet.steps1,
                                rep(NA, nEmpty))
        colData(sim)$Path2 <- factor(c(rep(NA, nCells), doublet.paths2,
                                       rep(NA, nEmpty)))
        colData(sim)$Step2 <- c(rep(NA, nCells), doublet.steps2,
                                rep(NA, nEmpty))
    }

    assays(sim)$CellMeans <- cells.means

    return(sim)
}

#' Simulate Kersplat cell counts
#'
#' Simulate cell counts for the Kersplat simulation
#'
#' @param sim SingleCellExperiment containing simulation.
#' @param params KersplatParams object with simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' Counts are sampled from a Poisson distribution with lambda equal to the
#' cell means matrix.
#'
#' @return SingleCellExperiment with cell counts
kersplatSimCellCounts <- function(sim, params, verbose) {

    if (verbose) {message("Simulating cell counts...")}
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    nEmpty <- getParam(params, "ambient.nEmpty")
    cells.means <- assays(sim)$CellMeans

    cell.counts <- matrix(rpois(
        as.numeric(nGenes) * as.numeric(nCells + nEmpty),
        lambda = cells.means),
        nrow = nGenes, ncol = nCells + nEmpty)

    colnames(cell.counts) <- cell.names
    rownames(cell.counts) <- gene.names
    assays(sim)$CellCounts <- cell.counts

    return(sim)
}

#' Simulate Kersplat ambient counts
#'
#' @param sim SingleCellExperiment containing simulation.
#' @param params KersplatParams object with simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' The overall expression profile to calculated by averaging the cell counts
#' of the (non-empty) cells. This is then multiplied by the ambient library
#' sizes to get a mean for each cell. Counts are then sampled from a Poisson
#' distribution using these means.
#'
#' @return SingleCellExperiment with ambient counts
kersplatSimAmbientCounts <- function(sim, params, verbose) {

    if (verbose) {message("Simulating ambient counts...")}
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    nEmpty <- getParam(params, "ambient.nEmpty")
    cell.counts <- assays(sim)$CellCounts
    not.empty <- colData(sim)$Type != "Empty"
    ambient.lib.sizes <- colData(sim)$AmbientLibSize

    not.empty.means <- rowMeans(cell.counts[, not.empty])
    ambient.props <- not.empty.means / sum(not.empty.means)

    ambient.means <- ambient.props %*% t(ambient.lib.sizes)

    ambient.counts <- matrix(rpois(
        as.numeric(nGenes) * as.numeric(nCells + nEmpty),
        lambda = ambient.means),
        nrow = nGenes, ncol = nCells + nEmpty)

    colnames(ambient.counts) <- cell.names
    rownames(ambient.counts) <- gene.names
    assays(sim)$AmbientCounts <- ambient.counts
    rowData(sim)$AmbientMean <- not.empty.means

    return(sim)
}

#' Simulate Kersplat final counts
#'
#' Simulate the final counts matrix for a Kersplat simulation
#'
#' @param sim SingleCellExperiment containing simulation.
#' @param params KersplatParams object with simulation parameters.
#' @param verbose logical. Whether to print progress messages
#'
#' @details
#' The cell counts matrix and ambient counts matrix are added together. The
#' result is then downsampled to the cell library size (for cells and doublets)
#' or the ambient library size (for empty cells) using the
#' \code{\link[scuttle]{downsampleMatrix}} function.
#'
#' @seealso \code{\link[scuttle]{downsampleMatrix}}
#'
#' @return SingleCellExperiment with counts matrix
kersplatSimCounts <- function(sim, params, verbose) {

    if (verbose) {message("Simulating final counts...")}
    cell.lib.sizes <- colData(sim)$CellLibSize
    ambient.lib.sizes <- colData(sim)$AmbientLibSize
    empty <- colData(sim)$Type == "Empty"
    cell.counts <- assays(sim)$CellCounts
    ambient.counts <- assays(sim)$AmbientCounts

    lib.sizes <- cell.lib.sizes
    lib.sizes[empty] <- ambient.lib.sizes[empty]

    counts <- cell.counts + ambient.counts

    down.prop <- lib.sizes / colSums(counts)
    # Avoid proportion creeping over 1 for empty cells
    down.prop <- min(down.prop, 1)

    counts <- scuttle::downsampleMatrix(counts, down.prop)

    assays(sim)$counts <- counts

    return(sim)
}

#' Get Beta step probabilities
#'
#' Use a Beta distribution for set probabilities along a path
#'
#' @param steps Number of steps
#' @param alpha Alpha parameter
#' @param beta Beta parameter
#'
#' @details
#' The density is sampled from a Beta distribution between 0 and 1. Infinite
#' densities at edges are adjusted and then the values are scaled to give
#' probabilities.
#'
#' @return Vector of probabilities
#'
#' @importFrom stats dbeta
getBetaStepProbs <- function(steps, alpha, beta) {
    dens <- dbeta(seq(0, 1, length.out = steps), alpha, beta)

    # Adjust for infinite values at edge of distribution
    dens.inf <- !is.finite(dens)
    if (any(dens.inf) && all(dens[!dens.inf] == 0)) {
        dens[dens.inf] <- 1
    }
    if (!is.finite(dens[1])) {
        dens[1] <- 1.1 * dens[2]
    }
    if (!is.finite(dens[steps])) {
        dens[steps] <- 1.1 * dens[steps - 1]
    }

    probs <- dens / sum(dens)

    return(probs)
}

#' Sample density
#'
#' Sample from a density object using rejection sampling
#'
#' @param n Number of values to sample
#' @param dens Density object to sample from
#' @param lower Lower x-axis bound on sampled values
#'
#' @details
#' Random points (x and y) are generated inside the range of the density object.
#' If they value is less than the density for that x value (and x is greater
#' than \code{lower}) then that x value is retained. Ten thousand points are
#' generated at a time until enough valid values have been sampled.
#'
#' @return Vector of sampled values
#'
#' @importFrom stats approxfun
sampleDensity <- function(n, dens, lower = 0) {

    xmin <- min(dens$x)
    xmax <- max(dens$x)
    ymin <- min(dens$y)
    ymax <- max(dens$y)

    boundary <- approxfun(dens$x, dens$y)

    values <- c()
    nsel <- 0

    while(nsel < n) {
        x <- runif(1e4, xmin, xmax)
        y <- runif(1e4, ymin, ymax)
        sel <- y < boundary(x) & x > lower

        nsel <- nsel + sum(sel)
        values <- c(values, x[sel])
    }

    values <- values[seq_len(n)]

    return(values)
}
