#' @rdname neweQTLParams
#' @importFrom methods new
#' @export
neweQTLParams <- function(...) {
    
    eQTLparams <- new("eQTLParams")
    eQTLparams <- setParams(eQTLparams, ...)

    return(eQTLparams)
}


#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric
#' checkFlag
setValidity("eQTLParams", function(object) {
    
    object <- expandParams(object)
    v <- getParams(object, c(slotNames(object)))
    
    checks <- c(esnp.n = checkInt(v$esnp.n, lower = 1),
                eqtl.dist = checkInt(v$eqtl.dist, lower = 1),
                eqtl.maf = checkNumber(v$eqtl.maf, lower = 0, upper = 1),
                eqtl.mafd = checkNumber(v$eqtl.mafd, lower = 0, upper = 1),
                eqtlES.shape = checkNumber(v$eqtlES.shape, lower = 0),
                eqtlES.rate = checkNumber(v$eqtlES.rate, lower = 0),
                bulkmean.shape = checkNumber(v$bulkmean.shape, lower = 0),
                bulkmean.rate = checkNumber(v$bulkmean.rate, lower = 0),
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
    
    pp <- list("eQTL.Params:"      = c("[pair Number]"= "esnp.n",
                                       "[pair Distance]"= "eqtl.dist",
                                       "[eSNP MAF]"= "eqtl.maf",
                                       "[eSNP MAF dev]"= "eqtl.mafd"),
               "eQTL.Effect Size:"= c("(Rate)"        = "eqtlES.rate",
                                      "(Shape)"        = "eqtlES.shape"),
               "eQTL.Mean:"      = c("(Rate)"         = "bulkmean.rate",
                                     "(Shape)"        = "bulkmean.shape",
                                     "(CV params)"    = "bulkcv.param"))
    
    callNextMethod()
    showPP(object, pp)
})


#' @rdname setParam
setMethod("setParam", "eQTLParams", function(object, name, value) {
    checkmate::assertString(name)
    
    # Function to check that user defined changes to eQTLParams don't break rules
    
    # For example:
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
    
    

    object <- callNextMethod()
    
    return(object)
})

