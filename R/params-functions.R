getParams <- function(params, names) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertCharacter(names, min.len = 1, any.missing = FALSE)

    sapply(names, getParam, object = params, simplify = FALSE)
}

setParams <- function(params, update = NULL, ...) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(update, null.ok = TRUE)

    update <- c(update, list(...))

    if (length(update) > 0) {
        for (name in names(update)) {
            value <- update[[name]]
            params <- setParam(params, name, value)
        }
    }

    return(params)
}

showPP <- function(params, pp) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(pp, types = "character", min.len = 1)

    default <- new(class(params))
    for (category in names(pp)) {
        parameters <- pp[[category]]
        values <- getParams(params, parameters)
        values <- sapply(values, paste, collapse = ", ")
        default.values <- getParams(default, parameters)
        default.values <- sapply(default.values, paste, collapse = ", ")
        not.default <- values != default.values
        names(values)[not.default] <- toupper(names(values[not.default]))
        cat(category, "\n")
        print(noquote(values), print.gap = 2)
        cat("\n")
    }
}

# mergeParams <- function(params1, params2) {
#
#     if (class(params1) != class(params2)) {
#         stop("params1 and params2 must be of the same Params class")
#     }
#
#     default <- new(class(params1))
#
#     update <- list()
#     for (parameter in slotNames(params1)) {
#         value1 <- getParam(params1, parameter)
#         default.value <- getParam(default, parameter)
#         if (value1 == default.value) {
#             value2 <- getParam(params2, parameter)
#             update[[parameter]] <- value2
#         } else {
#             update[[parameter]] <- value1
#         }
#     }
#
#     merged <- setParams(default, update)
#
#     return(merged)
# }