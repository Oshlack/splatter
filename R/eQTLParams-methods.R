#' @rdname newParams
#' @importFrom methods new
#' @export
neweQTLParams <- function(...) {
    
    params <- new("eQTLParams")
    params <- setParams(params, ...)

    return(params)
}


#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric
#' checkFlag
setValidity("eQTLParams", function(object) {
    
    v <- getParams(object, c(slotNames(object)))
    
    checks <- c(eqtl.n = checkInt(v$eqtl.n, lower = 1),
                eqtl.dist = checkInt(v$eqtl.dist, lower = 1),
                eqtl.maf = checkNumber(v$eqtl.maf, lower = 0, upper = 1),
                eqtl.mafd = checkNumber(v$eqtl.mafd, lower = 0, upper = 1),
                eqtlES.shape = checkNumber(v$eqtlES.shape, lower = 0),
                eqtlES.rate = checkNumber(v$eqtlES.rate, lower = 0),
                bulkmean.shape = checkNumber(v$bulkmean.shape, lower = 0),
                bulkmean.rate = checkNumber(v$bulkmean.rate, lower = 0),
                bulkcv.bins = checkInt(v$bulkcv.bins, lower=1),
                bulkcv.param = checkDataFrame(v$bulkcv.param))
    
    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }
    
    return(valid)
})


#' @importFrom methods callNextMethod
setMethod("show", "eQTLParams", function(object) {
    
    pp <- list("eQTL.General:"      = c("[eQTL.N]"    = "eqtl.n",
                                       "[Distance]"   = "eqtl.dist",
                                       "[MAF]"        = "eqtl.maf",
                                       "[MAF dev]"    = "eqtl.mafd"),
               "eQTL.Effect Size:"= c("(Rate)"        = "eqtlES.rate",
                                      "(Shape)"       = "eqtlES.shape"),
               "eQTL.Mean:"      = c("(Rate)"         = "bulkmean.rate",
                                     "(Shape)"        = "bulkmean.shape",
                                     "[CV bins]"      = "bulkcv.bins",
                                     "(CV params)"    = "bulkcv.param"))
    
    callNextMethod()
    showPP(object, pp)
})


#' @rdname setParam
setMethod("setParam", "eQTLParams", function(object, name, value) {
    checkmate::assertString(name)
    
    if (name == "nCells") {
        warning(name, " parameter does not impact eQTL simulation, set nCells
             in the sc simulation Params object instead.")
    }
    
    if (name == "bulkcv.param") {
        if (getParam(object, "bulkcv.bins") != nrow(value)) {
            stop("Need to set bulkcv.bins to length of bulkcv.param")
        }
    }
    
    object <- callNextMethod()
    
    return(object)
})

#' @rdname setParams
setMethod("setParams", "eQTLParams", function(object, update = NULL, ...) {
    
    checkmate::assertClass(object, classes = "eQTLParams")
    checkmate::assertList(update, null.ok = TRUE)
    
    update <- c(update, list(...))
    
    object <- callNextMethod(object, update)
    
    return(object)
})