#' @rdname getParam
#' @importFrom methods slot
setMethod("getParam", "Params", function(object, name) {
    slot(object, name)
})

#' @rdname setParam
#' @importFrom methods slot<- validObject
setMethod("setParam", "Params", function(object, name, value) {
    checkmate::assertString(name)
    slot(object, name) <- value
    validObject(object)
    return(object)
})

#' @rdname setParamUnchecked
#' @importFrom methods slot<-
setMethod("setParamUnchecked", "Params", function(object, name, value) {
    checkmate::assertString(name)
    slot(object, name) <- value
    return(object)
})

#' @importFrom methods slotNames
setMethod("show", "Params", function(object) {

    pp <- list("Global:" = c("(Genes)" = "nGenes",
                             "(Cells)" = "nCells",
                             "[Seed]"  = "seed"))

    cat("A", crayon::bold("Params"), "object of class",
        crayon::bold(class(object)), "\n")
    cat("Parameters can be (estimable) or", crayon::blue("[not estimable],"),
        "'Default' or ", crayon::bold(crayon::green("'NOT DEFAULT'")), "\n\n")
    showPP(object, pp)
    cat(length(slotNames(object)) - 3, "additional parameters", "\n\n")
})

setMethod("expandParams", "Params", function(object, vectors, n) {

    update <- list()
    for (parameter in vectors) {
        value <- getParam(object, parameter)
        if (length(value) == 1) {
            update[[parameter]] <- rep(value, n)
        }
    }

    object <- setParamsUnchecked(object, update)

    return(object)
})
