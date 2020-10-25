#' @rdname newParams
#' @importFrom methods new
#' @export
newKersplatParams <- function(...) {

    checkDependencies("kersplat")

    if (getOption("splatter.warn.kersplat", TRUE)) {
        warning("The Kersplat simulation is still experimental and may ",
                "produce unreliable results. Please try it and report any ",
                "issues to ",
                "https://github.com/Oshlack/splatter/issues. The development ",
                "version may have improved features.")
        options(splatter.warn.kersplat = FALSE)
    }

    params <- new("KersplatParams")
    params <- setParams(params, ...)

    return(params)
}

setValidity("KersplatParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                seed = checkmate::checkInt(v$seed, lower = 0),
                mean.rate = checkmate::checkNumber(v$mean.rate, lower = 0),
                mean.shape = checkmate::checkNumber(v$mean.shape, lower = 0),
                mean.outProb = checkNumber(v$mean.outProb, lower = 0,
                                           upper = 1),
                mean.outLoc = checkNumber(v$mean.outLoc),
                mean.outScale = checkNumber(v$mean.outScale, lower = 0),
                mean.dens = checkmate::checkClass(v$mean.dens, "density"),
                mean.method = checkmate::checkChoice(v$mean.method,
                                                     c("fit", "density")),
                bcv.common = checkNumber(v$bcv.common, lower = 0),
                bcv.df = checkNumber(v$bcv.df, lower = 0),
                network.graph = checkmate::checkClass(v$network, "igraph",
                                                      null.ok = TRUE),
                network.nRegs = checkmate::checkInt(v$network.nRegs,
                                                    lower = 0),
                network.regsSet = checkmate::checkFlag(v$network.regsSet),
                paths.nPrograms = checkmate::checkInt(v$paths.nPrograms,
                                                      lower = 1),
                paths.design = checkmate::checkDataFrame(v$paths.design,
                                                         types = "numeric",
                                                         any.missing = FALSE,
                                                         min.rows = 1,
                                                         ncols = 3),
                lib.loc = checkmate::checkNumber(v$lib.loc),
                lib.scale = checkmate::checkNumber(v$lib.scale, lower = 0),
                lib.dens = checkmate::checkClass(v$lib.dens, "density"),
                lib.method = checkmate::checkChoice(v$lib.method,
                                                    c("fit", "density")),
                cells.design = checkmate::checkDataFrame(v$cells.design,
                                                         types = "numeric",
                                                         any.missing = FALSE,
                                                         nrows = nrow(
                                                             v$paths.design),
                                                         ncols = 4),
                doublet.prop = checkmate::check_number(v$doublet.prop,
                                                       lower = 0,
                                                       upper = 1),
                ambient.scale = checkmate::check_number(v$ambient.scale,
                                                        lower = 0,
                                                        upper = 1),
                ambient.nEmpty = checkmate::check_number(v$ambient.nEmpty,
                                                         lower = 0,
                                                         finite = TRUE))

    if (!checkmate::testNumeric(v$mean.values, len = 0)) {
        checks <- c(checks,
                    mean.values = checkmate::checkNumeric(v$mean.values,
                                                          lower = 0,
                                                          finite = TRUE,
                                                          any.missing = FALSE,
                                                          len = v$nGenes))
    }

    if (!all(colnames(v$paths.design) == c("Path", "From", "Steps"))) {
        checks <- c(checks, paths.design = "Incorrect column names")
    } else {
        if (!(0 %in% v$paths.design$From)) {
            checks <- c(checks, paths.design = paste("origin must be specified",
                      "in paths.design"))
        }

        paths.graph <- igraph::graph_from_data_frame(v$paths.design)
        if (!igraph::is_simple(paths.graph)) {
            checks <- c(checks, paths.design = "graph is not simple")
        }
        if (!igraph::is_connected(paths.graph)) {
            checks <- c(checks, paths.design = "graph is not connected")
        }
        if (!igraph::is_dag(paths.graph)) {
            checks <- c(checks, paths.design = "graph is not a DAG")
        }
    }

    if (!checkmate::testList(v$paths.means, len = 0)) {
        checks <- c(checks,
                    paths.means = checkmate::checkList(v$paths.means,
                                                       types = "matrix",
                                                       any.missing = FALSE,
                                                       len = nrow(v$paths.design),
                                                       names = "unique"))
    }

    if (!all(colnames(v$cells.design) == c("Path", "Probability", "Alpha",
                                           "Beta"))) {
        checks <- c(checks, cells.design = "Incorrect column names")
    } else {
        if (!all(sort(v$cells.design$Path) == sort(v$paths.design$Path))) {
            checks <- c(checks,
                        cells.design = "Path names don't match paths.design")
        }
        if (sum(v$cells.design$Probability) != 1) {
            checks <- c(checks, cells.design = "Probability must sum to 1")
        }
        checks <- c(checks,
                    "cells.design$Alpha" = checkmate::checkNumeric(
                        v$cells.design$Alpha,
                        lower = 0, any.missing = FALSE))
        checks <- c(checks,
                    "cells.design$Beta" = checkmate::checkNumeric(
                        v$cells.design$Beta,
                        lower = 0, any.missing = FALSE))
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
setMethod("show", "KersplatParams", function(object) {

    pp.top <- list("Mean:" = c("(Rate)"         = "mean.rate",
                               "(Shape)"        = "mean.shape",
                               "(Out Prob)"     = "mean.outProb",
                               "(Out Location)" = "mean.outLoc",
                               "(Out Scale)"    = "mean.outScale",
                               "(Density)"      = "mean.dens",
                               "[Method]"       = "mean.method",
                               "[Values]*"      = "mean.values"),
                   "BCV:"  = c("(Common Disp)"  = "bcv.common",
                               "[DoF]"          = "bcv.df"))

    pp.network <- list("Network:" = c("[Graph]*"   = "network.graph",
                                      "[nRegs]"    = "network.nRegs",
                                      "[regsSet]*" = "network.regsSet"))

    pp.paths <- list("Paths:" = c("[nPrograms]" = "paths.nPrograms",
                                  "[Design]"  = "paths.design"))

    pp.bot <- list("Library size:" = c("(Location)" = "lib.loc",
                                       "(Scale)"    = "lib.scale",
                                       "(Density)"  = "lib.dens",
                                       "[Method]"   = "lib.method"),
                   "Cells:"        = c("[Design]"   = "cells.design"),
                   "Doublets:"     = c("[Prop]"     = "doublet.prop"),
                   "Ambient:"      = c("[Scale]"    = "ambient.scale",
                                       "[Empty]"    = "ambient.nEmpty"))

    paths.means <- getParam(object, "paths.means")
    if (length(paths.means) == 0) {
        pp.paths[[1]] <- c(pp.paths[[1]], "[Means]*" = "paths.means")
    }

    callNextMethod()

    showPP(object, pp.top)

    network.graph <- getParam(object, "network.graph")
    if (is.null(network.graph)) {
        showPP(object, pp.network)
    } else {
        cat(crayon::bold("Network:"), "\n")
        cat(crayon::bold(crayon::bgYellow(crayon::blue("[GRAPH]\n"))))
        cat(crayon::bold(crayon::green(paste(
            "Graph with", igraph::gorder(network.graph), "nodes and",
            igraph::gsize(network.graph), "edges\n"
        ))))
        show(network.graph)
        network.values <- list("[nRegs]" = getParam(object, "network.nRegs"),
                               "[regsSet]*" = getParam(object,
                                                       "network.regsSet"))
        network.default <- c(network.values$`[nRegs]` != 100,
                             network.values$`[regsSet]` != FALSE)
        showValues(network.values, network.default)
    }

    showPP(object, pp.paths)

    if (length(paths.means) != 0) {
        cat(crayon::bgYellow(crayon::bold(crayon::blue("[MEANS]\n"))))
        cat(crayon::bold(crayon::green(paste(
            "List of", length(paths.means), "matrices with names:",
            paste(names(paths.means), collapse = ", "), "\n\n"
        ))))
    }

    showPP(object, pp.bot)
})

#' @rdname setParam
setMethod("setParam", "KersplatParams", function(object, name, value) {
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
        if (getParam(object, "nGenes") != igraph::gorder(value)) {
            if (!(length(getParam(object, "mean.values")) == 0)) {
                warning("changing network.graph resets mean.values")
                object <- setParam(object, "mean.values", numeric())
            }
            object <- setParamUnchecked(object, "nGenes", igraph::gorder(value))
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

    if (name == "paths.design" &&
        (nrow(value) != nrow(getParam(object, "cells.design")))) {
        warning("cells.design reset to match paths.design")
        cells.design <- data.frame(Path = value$Path,
                                   Probability = 1 / nrow(value),
                                   Alpha = 1, Beta = 0)
        object <- setParamUnchecked(object, "cells.design", cells.design)
    }

    object <- callNextMethod()

    return(object)
})

#' @rdname setParams
setMethod("setParams", "KersplatParams", function(object, update = NULL, ...) {

    checkmate::assertList(update, null.ok = TRUE)

    update <- c(update, list(...))

    # If both cells.design and paths.design are given set cells.design first
    # to avoid reset warning
    if ("cells.design" %in% names(update) &&
        "paths.design" %in% names(update)) {
        object <- setParamUnchecked(object, "cells.design",
                                    update$cells.design)
        update$cells.design <- NULL
    }

    update <- bringItemsForward(update, c("network.graph", "paths.design"))

    object <- callNextMethod(object, update)

    return(object)
})
