#' @importFrom utils globalVariables
utils::globalVariables(c("pval_nominal", 'perc', 'gene_id',  'type', 'eQTL',
                         'chr', 'gene_mid', "maf", "data", "vcf", "gff", "dna"))

#' Logistic function
#'
#' Implementation of the logistic function
#'
#' @param x value to apply the function to.
#' @param x0 midpoint parameter. Gives the centre of the function.
#' @param k shape parameter. Gives the slope of the function.
#'
#' @return Value of logistic function with given parameters
logistic <- function(x, x0, k) {
    1 / (1 + exp(-k * (x - x0)))
}

#' Bind rows (matched)
#'
#' Bind the rows of two data frames, keeping only the columns that are common
#' to both.
#'
#' @param df1 first data.frame to bind.
#' @param df2 second data.frame to bind.
#'
#' @return data.frame containing rows from \code{df1} and \code{df2} but only
#'         common columns.
rbindMatched <- function(df1, df2) {
    common.names <- intersect(colnames(df1), colnames(df2))
    if (length(common.names) < 2) {
        stop("There must be at least two columns in common")
    }
    combined <- rbind(df1[, common.names], df2[, common.names])

    return(combined)
}

#' Bring items forward
#'
#' Move selected items to the start of a list.
#'
#' @param ll list to adjust item order.
#' @param items vector of items to bring to the front. Any not in the list will
#'        be ignored.
#'
#' @return list with selected items first
bringItemsForward <- function(ll, items) {

    checkmate::check_list(ll, min.len = 1, names = "unique")
    checkmate::check_character(items, any.missing = FALSE, min.len = 1,
                               unique = TRUE)

    items <- items[items %in% names(ll)]

    if (length(items) > 0) {
        ll.front <- ll[items]
        ll.back <- ll[!(names(ll) %in% items)]

        ll <- c(ll.front, ll.back)
    }

    return(ll)
}

#' Winsorize vector
#'
#' Set outliers in a numeric vector to a specified percentile.
#'
#' @param x Numeric vector to winsorize
#' @param q Percentile to set from each end
#'
#' @return Winsorized numeric vector
winsorize <- function(x, q) {

    checkmate::check_numeric(x, any.missing = FALSE)
    checkmate::check_number(q, lower = 0, upper = 1)

    lohi <- stats::quantile(x, c(q, 1 - q), na.rm = TRUE)

    if (diff(lohi) < 0) { lohi <- rev(lohi) }

    x[!is.na(x) & x < lohi[1]] <- lohi[1]
    x[!is.na(x) & x > lohi[2]] <- lohi[2]

    return(x)
}

#' Calculate coefficient of variation
#'
#' Implementation of the coefficient of variation
#'
#' @param x vector of values.
#'
#' @return Value of coefficient of variation for vector
#' @importFrom stats sd
co.var <- function(x) {
    stats::sd(x)/mean(x) 
} 


