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
           function(object, name, value) {standardGeneric("setParam")})

#' #@name arrange
#' #@rdname arrange
#' #@docType methods
#' #@export
#setGeneric("expandParams", function(object) {standardGeneric("expandParams")})