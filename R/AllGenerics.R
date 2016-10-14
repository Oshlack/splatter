#' New Params
#'
#' Create a new Params object. Functions exist for each of the different
#' Params subtypes.
#'
#' @param ... additional parameters passed to \code{\link{setParams}}.
#'
#' @return New Params object.
#' @examples
#' params <- newSimpleParams()
#' params <- newSimpleParams(nGenes = 200, nCells = 10)
#'
#' @name newParams
NULL

#' Get a parameter
#'
#' Accessor function for getting parameter values.
#'
#' @param object object to get parameter from.
#' @param name name of the parameter to get.
#'
#' @return The extracted parameter value
#'
#' @rdname getParam
#' @export
setGeneric("getParam", function(object, name) {standardGeneric("getParam")})

#' Set a parameter
#'
#' Function for setting parameter values.
#'
#' @param object object to set parameter in.
#' @param name name of the parameter to set.
#' @param value value to set the paramter to.
#'
#' @return Object with new parameter value.
#'
#' @rdname setParam
#' @export
setGeneric("setParam",
           function(object, name, value) {
               standardGeneric("setParam")
})

#' Set a parameter UNCHECKED
#'
#' Function for setting parameter values. THE OUTPUT IS NOT CHECKED FOR
#' VALIDITY!
#'
#' @param object object to set parameter in.
#' @param name name of the parameter to set.
#' @param value value to set the paramter to.
#'
#' @return Object with new parameter value.
#'
#' @rdname setParamUnchecked
setGeneric("setParamUnchecked",
           function(object, name, value) {
               standardGeneric("setParamUnchecked")
           })

#' Expand parameters
#'
#' Expand the parameters that can be vectors so that they are the same length as
#' the number of groups.
#'
#' @param object object to expand.
#' @param ... additional arguments.
#'
#' @return Expanded object.
#'
#' @rdname expandParams
setGeneric("expandParams", function(object, ...) {
    standardGeneric("expandParams")
})