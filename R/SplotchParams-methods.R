#' @rdname newParams
#' @importFrom methods new
#' @export
newSplotchParams <- function(...) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("The Splotch simulation requires the 'igraph' package.")
    }

    params <- new("SplotchParams")
    params <- setParams(params, ...)

    return(params)
}

setValidity("SplotchParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                seed = checkmate::checkInt(v$seed, lower = 0),
                network.graph = checkmate::checkClass(v$network, "igraph",
                                                      null.ok = TRUE),
                network.nRegs = checkmate::checkInt(v$network.nRegs,
                                                    lower = 0))

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

setMethod("show", "SplotchParams", function(object) {

    pp.network <- list("Network:" = c("[Graph]" = "network.graph",
                                      "[nRegs]" = "network.nRegs"))

    callNextMethod()

    network.graph <- getParam(object, "network.graph")
    if (is.null(network.graph)) {
        showPP(object, pp.network)
    } else {
        cat(crayon::bold("Network:"), "\n")
        cat(crayon::bold(crayon::blue("[GRAPH]\n")))
        cat(crayon::bold(crayon::green(paste(
            "Graph with", igraph::gorder(network.graph), "nodes and",
            igraph::gsize(network.graph), "edges\n"
        ))))
        show(network.graph)
        network.nRegs <- getParam(object, "network.nRegs")
        names(network.nRegs) <- c("[nRegs]")
        if (network.nRegs == 100) {
            showValues(network.nRegs, FALSE)
        } else {
            showValues(network.nRegs, TRUE)
        }
    }

    # showPP(object, pp)
})
