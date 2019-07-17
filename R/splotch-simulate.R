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
#' @importFrom stats dbeta
splotchSimulate <- function(params = newSplotchParams(), verbose = TRUE, ...) {

    checkmate::assertClass(params, "SplotchParams")
    params <- setParams(params, ...)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nGenes <- getParam(params, "nGenes")
    network.graph <- getParam(params, "network.graph")

    # Generate network
    if (is.null(network.graph)) {
        network.graph <- generateNetwork(nGenes, verbose)
        params <- setParam(params, "network.graph", network.graph)
    }

    # Select regulators
    if (!getParam(params, "network.regsSet")) {
        network.nRegs <- getParam(params, "network.nRegs")
        network.graph <- selectRegulators(network.graph, network.nRegs,
                                          verbose)
        params <- setParam(params, "network.graph", network.graph)
    } else {
        if (verbose) {message("Using selected regulators...")}
    }

    # Generate means
    if (length(getParam(params, "mean.values")) == 0) {
        if (verbose) {message("Simulating means...")}
        mean.shape <- getParam(params, "mean.shape")
        mean.rate <- getParam(params, "mean.rate")
        mean.values <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)
        params <- setParam(params, "mean.values", mean.values)
    } else {
        if (verbose) {message("Using defined means...")}
    }

    # Generate paths
    if (verbose) {message("Simulating paths...")}
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
    # The origin is not a path
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
            programs.changes <- rnorm(paths.nPrograms)
            reg.changes <- as.vector(programs.weights %*% programs.changes)
            changes[network.isReg, step] <- reg.changes
            change <- as.vector(changes[, step - 1] %*% network.weights)
            changes[, step] <- changes[, step] + change
        }

        if (from == 0) {
            changes <- changes[, 1:nSteps]
            factors <- matrixStats::rowCumsums(changes)
        } else {
            changes <- changes[, 2:nSteps + 1]
            from.factors <- paths.factors[[from]][, ncol(paths.factors[[from]])]
            factors <- matrixStats::rowCumsums(changes) + from.factors
        }
        paths.changes[[path]] <- changes
        paths.factors[[path]] <- factors
    }
    means.values <- getParam(params, "mean.values")
    paths.means <- lapply(paths.factors, function(x) {
        2 ^ x * means.values
    })
    names(paths.means) <- paste0("Path", paths.design$Path)
    params <- setParam(params, "paths.means", paths.means)

    if (verbose) {message("Simulating library sizes...")}
    nCells <- getParam(params, "nCells")
    lib.loc <- getParam(params, "lib.loc")
    lib.scale <- getParam(params, "lib.scale")
    lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)

    if (verbose) {message("Simulating cell means...")}
    nCells <- getParam(params, "nCells")
    cells.design <- getParam(params, "cells.design")
    cells.paths <- sample(cells.design$Path, nCells, replace = TRUE,
                          prob = cells.design$Probability)
    paths.design <- getParam(params, "paths.design")
    paths.cells.design <- merge(paths.design, cells.design)
    steps.probs <- apply(paths.cells.design, 1, function(path) {
        steps <- path["Steps"]
        dens <- dbeta(seq(0, 1, length.out = steps),
                      path["Alpha"], path["Beta"])
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
    cells.means <- sapply(seq_len(nCells), function(cell) {
        path <- cells.paths[cell]
        step <- cells.steps[cell]
        means <- paths.means[[path]][, step]
        return(means)
    })

    # sim <- SingleCellExperiment(assays = list(counts = counts),
    #                             rowData = features,
    #                             colData = cells,
    #                             metadata = list(params = params))
    #
    # return(sim)

    return(params)
}

generateNetwork <- function(n.nodes, verbose) {

    if (verbose) {message("Generating gene network...")}

    graph.raw <- igraph::sample_forestfire(n.nodes, 0.1)
    graph.data <- igraph::get.data.frame(graph.raw)
    graph.data <- graph.data[, c("from", "to")]
    graph.data$weight <- rnorm(nrow(graph.data))
    graph <- igraph::graph.data.frame(graph.data)
    graph <- igraph::set_vertex_attr(graph, "mean",
                                     value = rnorm(igraph::gorder(graph)))

    return(graph)
}

selectRegulators <- function(graph, nReg, verbose) {

    if (verbose) {message("Selecting regulators...")}

    out.degree <- igraph::degree(graph, mode = "out")
    in.degree <- igraph::degree(graph, mode = "in")
    reg.prob <-  out.degree - in.degree
    reg.prob <- reg.prob + rnorm(length(reg.prob))
    reg.prob[reg.prob <= 0] <- 1e-10
    reg.prob <- reg.prob / sum(reg.prob)
    reg.nodes <- names(rev(sort(reg.prob))[1:nReg])
    is.reg <- igraph::V(graph)$name %in% reg.nodes
    graph <- igraph::set_vertex_attr(graph, "IsReg", value = is.reg)

    return(graph)
}
