#' @rdname newParams
#' @importFrom methods new
#' @export
newSparseDCParams <- function(...) {

    checkDependencies("sparseDC")

    params <- new("SparseDCParams")
    params <- setParams(params, ...)

    return(params)
}

setValidity("SparseDCParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                markers.n = checkmate::checkInt(v$markers.n, lower = 0,
                                                upper = v$nGenes),
                markers.shared = checkmate::checkInt(v$markers.shared,
                                                     lower = 0,
                                                     upper = v$markers.n),
                markers.same = checkmate::checkFlag(v$markers.same),
                clusts.c1 = checkmate::checkIntegerish(v$clusts.c1, lower = 1,
                                                       any.missing = FALSE,
                                                       min.len = 1),
                clusts.c2 = checkmate::checkIntegerish(v$clusts.c2, lower = 1,
                                                       any.missing = FALSE,
                                                       min.len = 1),
                mean.lower = checkmate::checkNumber(v$mean.lower,
                                                    finite = TRUE),
                mean.upper = checkmate::checkNumber(v$mean.upper,
                                                    finite = TRUE),
                seed = checkmate::checkInt(v$seed, lower = 0))

    if (length(v$clusts.c1) == 1) {
        if (v$clusts.c1 > 1) {
            checks <- c(checks,
                        clusts.c1 = "If clusts.c1 is length 1 it must equal 1")
        }
    }

    if (length(v$clusts.c2) == 1) {
        if (v$clusts.c2 > 1) {
            checks <- c(checks,
                        clusts.c2 = "If clusts.c2 is length 1 it must equal 1")
        }
    }

    clusts.combined <- c(v$clusts.c1, v$clusts.c2)
    if (!all(seq_len(max(clusts.combined)) %in% clusts.combined)) {
        checks <- c(checks, clusts = "Cluster labels must be sequential")
    }

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

setMethod("show", "SparseDCParams", function(object) {

    pp <- list("Markers:"  = c("(Number)"      = "markers.n",
                               "(Shared)"      = "markers.shared",
                               "[Same]"        = "markers.same"),
               "Clusters:" = c("(Condition 1)" = "clusts.c1",
                               "(Condition 2)" = "clusts.c2"),
               "Means:"    = c("[Lower]"       = "mean.lower",
                               "[Upper]"       = "mean.upper"))

    callNextMethod()
    showPP(object, pp)
})
