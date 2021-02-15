#' @rdname newParams
#' @importFrom methods new
#' @export
newSplatPopParams <- function(...) {

    #checkDependencies("splatPop")

    params <- new("SplatPopParams")
    params <- setParams(params, ...)

    return(params)
}


#' @importFrom checkmate checkInt checkIntegerish checkNumber checkNumeric
#' checkFlag
setValidity("SplatPopParams", function(object) {
    
    object <- expandParams(object)
    v <- getParams(object, c(slotNames(object)))
    
    nConditions <- v$nConditions
    
    checks <- c(eqtl.n = checkNumber(v$eqtl.n, lower = 0),
                eqtl.dist = checkInt(v$eqtl.dist, lower = 1),
                eqtl.maf.min = checkNumber(v$eqtl.maf.min, lower = 0,
                                           upper = 0.5),
                eqtl.maf.max = checkNumber(v$eqtl.maf.max, lower = 0,
                                           upper = 0.5),
                eqtl.ES.shape = checkNumber(v$eqtl.ES.shape, lower = 0),
                eqtl.ES.rate = checkNumber(v$eqtl.ES.rate, lower = 0),
                eqtl.group.specific = checkNumber(v$eqtl.group.specific,
                                                 lower = 0, upper = 1),
                eqtl.condition.specific = checkNumber(v$eqtl.condition.specific,
                                                  lower = 0, upper = 1),
                pop.mean.shape = checkNumber(v$pop.mean.shape, lower = 0),
                pop.mean.rate = checkNumber(v$pop.mean.rate, lower = 0),
                pop.quant.norm = checkFlag(v$pop.quant.norm),
                pop.cv.bins = checkInt(v$pop.cv.bins, lower = 1),
                pop.cv.param = checkDataFrame(v$pop.cv.param),
                similarity.scale = checkNumber(v$similarity.scale, lower = 0),
                batch.size = checkInt(v$batch.size, lower = 1),
                nConditions = checkInt(v$nConditions, lower = 1),
                condition.prob = checkNumeric(v$condition.prob, lower = 0, 
                                              upper = 1, len = nConditions),
                cde.prob = checkNumeric(v$cde.prob, lower = 0, upper = 1,
                                       len = nConditions),
                cde.downProb = checkNumeric(v$cde.downProb, lower = 0, 
                                            upper = 1, len = nConditions),
                cde.facLoc = checkNumeric(v$cde.facLoc, len = nConditions),
                cde.facScale = checkNumeric(v$cde.facScale, lower = 0,
                                           len = nConditions))

    # Check group.prob sums to 1
    if (sum(v$condition.prob) != 1) {
        checks <- c(checks, "condition.probs must sum to 1")
    }
    
    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})



#' @rdname setParam
setMethod("setParam", "SplatPopParams", function(object, name, value) {

    checkmate::assertString(name)

    # splatPopParam checks
    if (name == "pop.cv.param") {
        if (getParam(object, "pop.cv.bins") != nrow(value)) {
            stop("Need to set pop.cv.bins to length of pop.cv.param")
        }
    }

    if (name == "eqtl.maf.min") {
        if (getParam(object, "eqtl.maf.min") >= getParam(object, "eqtl.maf.max")) {
            stop("Range of acceptable Minor Allele Frequencies is too small...
                 Be sure eqtl.maf.min < eqtl.maf.max.")
        }
    }

    if (name == "nConditions" ) {
        stop(name, " cannot be set directly, set condition.prob instead")
    }
    
    if (name == "condition.prob") {
        object <- setParamUnchecked(object, "nConditions", length(value))
    }
    
    object <- callNextMethod()
    
    return(object)
})


#' @importFrom methods callNextMethod
setMethod("show", "SplatPopParams", function(object) {
    
    pp <- list("Population params:" = c("(mean.shape)" = "pop.mean.shape",
                                        "(mean.rate)" = "pop.mean.rate",
                                        "[pop.quant.norm]" = "pop.quant.norm",
                                        "[similarity.scale]" = "similarity.scale",
                                        "[batch.size]" = "batch.size",
                                        "[cv.bins]" = "pop.cv.bins",
                                        "(cv.params)" = "pop.cv.param"),
               "eQTL params:" = c("[eqtl.n]" = "eqtl.n",
                                  "[eqtl.dist]" = "eqtl.dist",
                                  "[eqtl.maf.min]" = "eqtl.maf.min",
                                  "[eqtl.maf.max]" = "eqtl.maf.max",
                                  "[eqtl.group.specific]" = 
                                      "eqtl.group.specific",
                                  "[eqtl.condition.specific]" = 
                                      "eqtl.condition.specific",
                                  "(eqtl.ES.shape)" = "eqtl.ES.shape",
                                  "(eqtl.ES.rate)" = "eqtl.ES.rate"),
               "Condition params:" = c("[nConditions]"= "nConditions",
                                       "[condition.prob]" = "condition.prob",
                                       "[cde.prob]"  = "cde.prob",
                                       "[cde.downProb]" = "cde.downProb",
                                       "[cde.facLoc]" = "cde.facLoc",
                                       "[cde.facScale]" = "cde.facScale"))
    
    callNextMethod()
    showPP(object, pp)
})


#' @rdname expandParams
setMethod("expandParams", "SplatPopParams", function(object) {

    n <- getParam(object, "nConditions")
    
    vectors <- c("cde.prob", "cde.downProb", "cde.facLoc", "cde.facScale")
    
    # I should be able to call expandParams(object, vectors, n) or 
    # callNextMethod() with those variables (like for expanding LunParams or 
    # BasicParams), but it gives me an error that I can't figure out: 
    # Error in .local(object, ...) : 
    #      unused arguments (c("cde.prob", "cde.downProb", "cde.facLoc", "cde.facScale"), 2)
    
    # Remove below after bug above fixed! 
    cde.prob <- getParam(object, "cde.prob")
    cde.downProb <- getParam(object, "cde.downProb")
    cde.facLoc <- getParam(object, "cde.facLoc")
    cde.facScale <- getParam(object, "cde.facScale")

    object <- setParamsUnchecked(object, list("cde.prob" = rep(cde.prob, n),
                                           "cde.downProb" = rep(cde.downProb, n),
                                           "cde.facLoc" = rep(cde.facLoc, n),
                                           "cde.facScale" = rep(cde.facScale, n)))
    # 
    #object <- expandParams(object, vectors, n)

    return(object)
})
