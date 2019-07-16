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

    if (length(getParam(params, "mean.values")) == 0) {
        if (verbose) {message("Simulating means...")}
        mean.shape <- getParam(params, "mean.shape")
        mean.rate <- getParam(params, "mean.rate")
        mean.values <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)
        params <- setParam(params, "mean.values", mean.values)
    } else {
        if (verbose) {message("Using defined means...")}
    }

    # Simulate base gene means

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
