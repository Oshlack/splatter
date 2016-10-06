#' Logistic function
#'
#' Implementation of the logistic function
#'
#' @param x value to apply the function to.
#' @param x0 midpoint parameter. Gives the centre of the function.
#' @param k shape parameter. Gives the slope of the function.
#'
#' @return RETURN DESCRIPTION
#' @examples
#' logistic(0, 1, 1)
logistic <- function(x, x0, k) {
    1 / (1 + exp(-k * (x - x0)))
}