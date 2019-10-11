#' Get parameters
#'
#' Get multiple parameter values from a Params object.
#'
#' @param params Params object to get values from.
#' @param names vector of names of the parameters to get.
#'
#' @return List with the values of the selected parameters.
#' @examples
#' params <- newSimpleParams()
#' getParams(params, c("nGenes", "nCells", "mean.rate"))
#' @export
getParams <- function(params, names) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertCharacter(names, min.len = 1, any.missing = FALSE)

    params.list <- lapply(names, getParam, object = params)
    names(params.list) <- names

    return(params.list)
}

#' Set parameters UNCHECKED
#'
#' Set multiple parameters in a Params object.
#'
#' @param params Params object to set parameters in.
#' @param update list of parameters to set where \code{names(update)} are the
#'        names of the parameters to set and the items in the list are values.
#' @param ... additional parameters to set. These are combined with any
#'        parameters specified in \code{update}.
#'
#' @details
#' Each parameter is set by a call to \code{\link{setParam}}. If the same
#' parameter is specified multiple times it will be set multiple times.
#' Parameters can be specified using a list via \code{update} (useful when
#' collecting parameter values in some way) or individually (useful when setting
#' them manually), see examples. THE FINAL OBJECT IS NOT CHECKED FOR VALIDITY!
#'
#' @return Params object with updated values.
setParamsUnchecked <- function(params, update = NULL, ...) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(update, null.ok = TRUE)

    update <- c(update, list(...))

    if (length(update) > 0) {
        for (name in names(update)) {
            value <- update[[name]]
            params <- setParamUnchecked(params, name, value)
        }
    }

    return(params)
}

#' Show pretty print
#'
#' Function used for pretty printing params object.
#'
#' @param params object to show.
#' @param pp list specifying how the object should be displayed.
#'
#' @return Print params object to console
showPP <- function(params, pp) {

    checkmate::assertClass(params, classes = "Params")
    checkmate::assertList(pp, types = "character", min.len = 1)

    default <- new(class(params))
    for (category in names(pp)) {
        parameters <- pp[[category]]
        values <- getParams(params, parameters)
        is.df <- vapply(values, is.data.frame, FALSE)

        default.values <- getParams(default, parameters)
        not.default <- vapply(seq_along(values), function(i) {
            !identical(values[i], default.values[i])
        }, FALSE)
        empty.values <- vapply(values, function(x) {
            is.null(x) || length(x) == 0
        }, FALSE)
        values[empty.values] <- "Not set"

        names(values) <- names(parameters)
        cat(crayon::bold(category), "\n")
        if (sum(!is.df) > 0) {
            showValues(values[!is.df], not.default[!is.df])
        }
        if (sum(is.df) > 0) {
            showDFs(values[is.df], not.default[is.df])
        }
        cat("\n")
    }
}

#' Show values
#'
#' Function used for pretty printing scalar or vector parameters.
#'
#' @param values list of values to show.
#' @param not.default logical vector giving which have changed from the default.
#'
#' @return Print values
#'
#' @importFrom utils head
showValues <- function(values, not.default) {

    checkmate::check_list(values, any.missing = FALSE, min.len = 1)
    checkmate::check_logical(not.default, any.missing = FALSE,
                             len = length(values))

    short.values <- vapply(values, function(x) {
        if (is.list(x)) {
            classes <- class(x)
            if (length(classes) == 1 && classes == "list") {
                paste("List with", length(x), "items")
            } else {
                paste("Object of class", paste(classes, collapse = ", "))
            }
        } else {
            if (length(x) > 4) {
                paste0(paste(head(x, n = 4), collapse = ", "), ",...")
            } else {
                paste(x, collapse = ", ")
            }
        }
    }, c(Value = "None"))

    names(short.values)[not.default] <- toupper(names(values[not.default]))

    max.len <- max(nchar(short.values), nchar(names(short.values)))

    screen.width <- options("width")$width
    items.per.line <- floor(screen.width / (max.len + 2))

    short.names <- names(short.values)
    not.est <- !grepl("\\(", short.names)
    secondary <- grepl("\\*", short.names)
    short.names <- gsub("\\*", "", short.names)

    short.names[not.est] <- crayon::blue(short.names[not.est])
    short.names[secondary] <- crayon::bgYellow(short.names[secondary])
    short.names[not.default] <- crayon::bold(short.names[not.default])
    short.values[not.default] <- crayon::green(short.values[not.default])
    short.values[not.default] <- crayon::bold(short.values[not.default])

    short.values <- crayon::col_align(short.values, max.len, "right")
    short.names <- crayon::col_align(short.names, max.len, "right")

    names(short.values) <- short.names

    values.list <- split(short.values,
                         ceiling(seq_along(short.values) / items.per.line))

    for (line in values.list) {
        cat(paste(names(line), collapse = "  "), "\n")
        cat(paste(unname(line), collapse = "  "), "\n")
    }
}

#' Show data.frame
#'
#' Function used for pretty printing data.frame parameters.
#'
#' @param dfs list of data.frames to show.
#' @param not.default logical vector giving which have changed from the default.
#'
#' @return Print data.frame parameters
#'
#' @importFrom utils head
showDFs <- function(dfs, not.default) {

    checkmate::check_list(dfs, types = "data.frame", any.missing = FALSE,
                          min.len = 1)
    checkmate::check_logical(not.default, any.missing = FALSE,
                             len = length(dfs))

    names(dfs)[not.default] <- toupper(names(dfs)[not.default])

    not.est <- !grepl("\\(", names(dfs))
    names(dfs)[not.est] <- crayon::blue(names(dfs)[not.est])
    names(dfs)[not.default] <- crayon::bold(names(dfs)[not.default])

    for (i in seq_along(dfs)) {
        df <- dfs[[i]]
        name <- names(dfs)[i]

        msg <- paste0("data.frame (", nrow(df), " x ", ncol(df),
                      ") with columns: ", paste(colnames(df), collapse = ", "))

        if (not.default[i]) {
            msg <- crayon::bold(crayon::green(msg))
        }

        cat(paste0("\n", name, "\n"))
        cat(msg, "\n")
        print(head(df, n = 4))
        if (nrow(df) > 4) {
            cat("# ... with", nrow(df) - 4, "more rows\n")
        }
    }
}
