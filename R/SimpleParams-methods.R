#' @rdname newParams
#' @importFrom methods new
#' @export
newSimpleParams <- function(...) {

    params <- new("SimpleParams")
    params <- setParams(params, ...)

    return(params)
}

setValidity("SimpleParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                mean.rate = checkmate::checkNumber(v$mean.rate, lower = 0),
                mean.shape = checkmate::checkNumber(v$mean.shape, lower = 0),
                count.disp = checkmate::checkNumber(v$count.disp, lower = 0),
                seed = checkmate::checkInt(v$seed, lower = 0))

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

#' @importFrom methods is
setMethod("show", "SimpleParams", function(object) {

    pp <- list("Mean:"   = c("(Rate)"       = "mean.rate",
                             "(Shape)"      = "mean.shape"),
               "Counts:" = c("[Dispersion]" = "count.disp"))

    # Mean parameters aren't estimated for the LunParams object which
    # inherits from SimpleParams
    if (is(object, "LunParams")) {
        pp[["Mean:"]] <- c("[Rate]" = "mean.rate", "[Shape]" = "mean.shape")
    }

    callNextMethod()
    showPP(object, pp)
})
