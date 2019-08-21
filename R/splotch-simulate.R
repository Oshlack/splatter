#' Splotch simulation
#'
#' Simulate counts from...
#'
#' @param params SplotchParams object containing simulation parameters.
#' @param verbose logical. Whether to print progress messages
#' @param ... any additional parameter settings to override what is provided in
#'        \code{params}.
#'
#' @details
#' Details...
#'
#' @return SingleCellExperiment containing simulated counts
#' @examples
#' sim <- splotchSimulate()
#'
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
splotchSimulate <- function(params = newSplotchParams(), verbose = TRUE, ...) {

    params <- splotchSetup(params, verbose, ...)
    sim <- splotchSample(params, verbose)

    return(sim)
}

splotchSetup <- function(params = newSplotchParams(), verbose = TRUE, ...) {

    checkmate::assertClass(params, "SplotchParams")
    params <- setParams(params, ...)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    if (verbose) {message("Setting up parameters...")}
    params <- splotchGenNetwork(params, verbose)
    params <- splotchSelectRegs(params, verbose)
    params <- splotchSimGeneMeans(params, verbose)
    params <- splotchSimPaths(params, verbose)

    return(params)
}

splotchSample <- function(params, verbose = TRUE) {

    # Check that parameters are set up
    checkmate::assertClass(params, "SplotchParams")
    network.graph <- getParam(params, "network.graph")
    if (is.null(network.graph)) {
        stop("'network.graph' not set, run splotchSetup first")
    }
    network.regsSet <- getParam(params, "network.regsSet")
    if (!network.regsSet) {
        stop("network regulators not set, run splotchSetup first")
    }
    mean.values <- getParam(params, "mean.values")
    if (length(mean.values) == 0) {
        stop("'mean.values' not set, run splotchSetup first")
    }
    paths.means <- getParam(params, "paths.means")
    if (length(mean.values) == 0) {
        stop("'paths.means' not set, run splotchSetup first")
    }

    if (verbose) {message("Creating simulation object...")}
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))

    cells <-  data.frame(Cell = cell.names, row.names = cell.names)
    features <- data.frame(Gene = gene.names, row.names = gene.names)
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))

    sim <- splotchSimLibSizes(sim, params, verbose)
    sim <- splotchSimCellMeans(sim, params, verbose)
    sim <- splotchSimCounts(sim, params, verbose)

    return(sim)

}

splotchGenNetwork <- function(params, verbose) {

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

splotchSelectRegs <- function(params, verbose) {

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
    reg.nodes <- names(rev(sort(reg.prob))[1:network.nRegs])
    is.reg <- igraph::V(network.graph)$name %in% reg.nodes
    network.graph <- igraph::set_vertex_attr(network.graph, "IsReg",
                                             value = is.reg)

    params <- setParam(params, "network.graph", network.graph)

    return(params)
}

splotchSimGeneMeans <- function(params, verbose) {

    mean.values <- getParam(params, "mean.values")

    # Generate means
    if (length(mean.values) != 0) {
        if (verbose) {message("Using defined means...")}
        return(params)
    }

    if (verbose) {message("Simulating means...")}
    nGenes <- getParam(params, "nGenes")
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

    params <- setParam(params, "mean.values", mean.values)

    return(params)
}

splotchSimPaths <- function(params, verbose) {

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
            changes <- changes[, 1:nSteps]
            factors <- matrixStats::rowCumsums(changes)
        } else {
            changes <- changes[, 2:(nSteps + 1)]
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

splotchSimLibSizes <- function(sim, params, verbose) {

    if (verbose) {message("Simulating library sizes...")}
    nCells <- getParam(params, "nCells")
    # lib.loc <- getParam(params, "lib.loc")
    # lib.scale <- getParam(params, "lib.scale")
    lib.dens <- getParam(params, "lib.dens")

    # exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    exp.lib.sizes <- rejectionSample(nCells, lib.dens)

    colData(sim)$ExpLibSize <- exp.lib.sizes

    return(sim)
}

#' @importFrom stats dbeta
splotchSimCellMeans <- function(sim, params, verbose) {

    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nCells <- getParam(params, "nCells")
    cells.design <- getParam(params, "cells.design")
    paths.design <- getParam(params, "paths.design")
    paths.means <- getParam(params, "paths.means")
    exp.lib.sizes <- colData(sim)$ExpLibSize

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

    cells.steps <- sapply(cells.paths, function(path) {
        probs <- steps.probs[[path]]
        step <- sample(1:length(probs), 1, prob = probs)
        return(step)
    })

    if (verbose) {message("Simulating cell means...")}
    cells.means <- sapply(seq_len(nCells), function(cell) {
        path <- cells.paths[cell]
        step <- cells.steps[cell]
        means <- paths.means[[path]][, step]
        return(means)
    })

    # Adjust mean based on library size
    cells.props <- t(t(cells.means) / colSums(cells.means))
    cells.means <- t(t(cells.props) * exp.lib.sizes)

    colnames(cells.means) <- cell.names
    rownames(cells.means) <- gene.names

    colData(sim)$Path <- cells.paths
    colData(sim)$Step <- cells.steps
    assays(sim)$CellMeans <- cells.means

    return(sim)
}

splotchSimCounts <- function(sim, params, verbose) {

    if (verbose) {message("Simulating counts...")}
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    nGenes <- getParam(params, "nGenes")
    nCells <- getParam(params, "nCells")
    cells.means <- assays(sim)$CellMeans

    true.counts <- matrix(rpois(
        as.numeric(nGenes) * as.numeric(nCells),
        lambda = cells.means),
        nrow = nGenes, ncol = nCells)

    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names
    assays(sim)$counts <- true.counts

    return(sim)
}

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

#' @importFrom stats approxfun
rejectionSample <- function(n, dens, lower = 0) {

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

    values <- values[1:n]

    return(values)
}
