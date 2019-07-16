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

    checks <- c(#nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                seed = checkmate::checkInt(v$seed, lower = 0),
                mean.rate = checkmate::checkNumber(v$mean.rate, lower = 0),
                mean.shape = checkmate::checkNumber(v$mean.shape, lower = 0),
                network.graph = checkmate::checkClass(v$network, "igraph",
                                                      null.ok = TRUE),
                network.nRegs = checkmate::checkInt(v$network.nRegs,
                                                    lower = 0),
                network.regsSet = checkmate::checkFlag(v$network.regsSet))

    if (checkmate::testNumeric(v$mean.values, len = 0)) {
        checks <- c(checks, mean.values = TRUE)
    } else {
        checks <- c(checks,
                    mean.values = checkmate::checkNumeric(v$mean.values,
                                                          lower = 0,
                                                          finite = TRUE,
                                                          any.missing = FALSE,
                                                          len = v$nGenes))
    }

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

#' @importFrom methods show
setMethod("show", "SplotchParams", function(object) {

    pp.top <- list("Mean:" = c("(Rate)"   = "mean.rate",
                               "(Shape)"  = "mean.shape",
                               "[Values]" = "mean.values"))

    pp.network <- list("Network:" = c("[Graph]"   = "network.graph",
                                      "[nRegs]"   = "network.nRegs",
                                      "[regsSet]" = "network.regsSet"))

    callNextMethod()

    showPP(object, pp.top)

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
        network.values <- list("[nRegs]" = getParam(object, "network.nRegs"),
                               "[regsSet]" = getParam(object,
                                                      "network.regsSet"))
        network.default <- c(network.values$`[nRegs]` != 100,
                             network.values$`[regsSet]` != FALSE)
        showValues(network.values, network.default)
    }
})

#' @rdname setParam
setMethod("setParam", "SplotchParams", function(object, name, value) {
    checkmate::assertString(name)

    if (name == "nGenes") {
        if (!is.null(getParam(object, "network.graph"))) {
            stop("nGenes can not be changed once network.graph is set")
        }
        if (!length(getParam(object, "mean.values")) == 0) {
            stop("nGenes can not be changed once mean.values is set")
        }
    }

    if (name == "mean.values") {
        if (!is.null(getParam(object, "network.graph"))) {
            if (length(value) != getParam(object, "nGenes")) {
                stop("new mean.values does not match the number of genes in ",
                     "the network")
            }
        } else {
            object <- setParam(object, "nGenes", length(value))
        }
    }

    if (name == "network.graph") {
        checkmate::assertClass(value, "igraph")
        object <- setParamUnchecked(object, "nGenes", igraph::gorder(value))
        if (!(length(getParam(object, "mean.values")) == 0)) {
            warning("changing network.graph resets mean.values")
            object <- setParam(object, "mean.values", numeric())
        }

        if ("IsReg" %in% igraph::vertex_attr_names(value)) {
            isReg <- igraph::vertex_attr(value, "IsReg")
            checkmate::assertLogical(isReg, any.missing = FALSE)
            object <- setParamUnchecked(object, "network.nRegs", sum(isReg))
            object <- setParamUnchecked(object, "network.regsSet", TRUE)
        } else {
            object <- setParamUnchecked(object, "network.regsSet", FALSE)
        }
    }

    if (name == "network.nRegs" && getParam(object, "network.regsSet")) {
        stop("network.nRegs can not be changed after regulators are set")
    }

    if (name == "network.regsSet") {
        stop("network.regsSet can not be set manually")
    }

    object <- callNextMethod()

    return(object)
})

#' @rdname setParams
setMethod("setParams", "SplotchParams", function(object, update = NULL, ...) {

    checkmate::assertList(update, null.ok = TRUE)

    update <- c(update, list(...))

    update <- bringItemsForward(update, c("network.graph"))

    object <- callNextMethod(object, update)

    return(object)
})
