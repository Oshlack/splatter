#' @rdname getParam
setMethod("getParam", "Params", function(object, name) {
    slot(object, name)
})

#' @rdname setParam
setMethod("setParam", "Params", function(object, name, value) {
    checkmate::assertString(name)
    slot(object, name) <- value
    validObject(object)
    return(object)
})

#' @rdname setParamUnchecked
setMethod("setParamUnchecked", "Params", function(object, name, value) {
              checkmate::assertString(name)
              slot(object, name) <- value
              return(object)
          })

setMethod("show", "Params", function(object) {

    pp <- list("Global:" = c("(Genes)" = "nGenes",
                             "(Cells)" = "nCells",
                             "[Seed]"  = "seed"))

    cat("A Params object of class", class(object), "\n")
    cat("Parameters can be (estimable) or [not estimable],",
        "'Default' or 'NOT DEFAULT'.", "\n\n")
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