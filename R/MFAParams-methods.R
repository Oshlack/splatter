#' @rdname newParams
#' @importFrom methods new
#' @export
newMFAParams <- function(...) {

    checkDependencies("mfa")

    params <- new("MFAParams")
    params <- setParams(params, ...)

    return(params)
}

setValidity("MFAParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                trans.prop = checkmate::checkNumber(v$trans.prop, lower = 0,
                                                    upper = 1),
                zero.neg = checkmate::checkLogical(v$zero.neg,
                                                   any.missing = FALSE,
                                                   len = 1),
                dropout.present = checkmate::checkLogical(v$dropout.present,
                                                          any.missing = FALSE,
                                                          len = 1),
                dropout.lambda = checkmate::checkNumber(v$dropout.lambda),
                seed = checkmate::checkInt(v$seed, lower = 0))

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

setMethod("show", "MFAParams", function(object) {

    pp <- list("Transient:" = c("[Proportion]" = "trans.prop"),
               "Negative:"  = c("[Zero]"       = "zero.neg"),
               "Dropout:"   = c("[Present]"    = "dropout.present",
                              "(Lambda)"       = "dropout.lambda"))

    callNextMethod()
    showPP(object, pp)
})
