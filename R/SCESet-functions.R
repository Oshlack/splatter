#' Add feature statistics
#'
#' Add additional feature statistics to an SCESet object
#'
#' @param sce SCESet to add feature statistics to.
#' @param value the expression value to calculate statistics for. Options are
#'        "counts", "cpm", "tpm" or "fpkm". The values need to exist in the
#'        given SCESet.
#' @param log logical. Whether to take log2 before calculating statistics.
#' @param offset offset to add to avoid taking log of zero.
#' @param no.zeros logical. Whether to remove all zeros from each feature before
#'        calculating statistics.
#'
#' @details
#' Currently adds the following statistics: mean, variance, coefficient of
#' variation, median and median absolute deviation. Statistics are added to
#' the \code{fData} slot and are named \code{stat_[log]_value_[no0]} where
#' \code{log} and \code{no0} are added if those arguments are true.
#'
#' @return SCESet with additional feature statistics
#'
#' @importFrom Biobase fData fData<-
addFeatureStats <- function(sce, value = c("counts", "cpm", "tpm", "fpkm"),
                            log = FALSE, offset = 1, no.zeros = FALSE) {

    value <- match.arg(value)

    switch(value,
           counts = {
               values = scater::counts(sce)
           },
           cpm = {
               values = scater::cpm(sce)
           },
           tpm = {
               values = scater::tpm(sce)
           },
           fpkm = {
               values = scater::fpkm(sce)
           }
    )

    suffix <- value

    if (no.zeros) {
        values[values == 0] <- NA
        suffix = paste0(suffix, "_no0")
    }

    if (log) {
        values = log2(values + offset)
        suffix = paste0("log_", suffix)
    }

    mean.str <- paste0("mean_", suffix)
    var.str  <- paste0("var_",  suffix)
    cv.str   <- paste0("cv_",   suffix)
    med.str  <- paste0("med_",  suffix)
    mad.str  <- paste0("mad_",  suffix)

    fData(sce)[, mean.str] <- rowMeans(values, na.rm = TRUE)
    fData(sce)[, var.str]  <- matrixStats::rowVars(values, na.rm = TRUE)
    fData(sce)[, cv.str]   <- sqrt(fData(sce)[, var.str]) /
        fData(sce)[, mean.str]
    fData(sce)[, med.str]  <- matrixStats::rowMedians(values, na.rm = TRUE)
    fData(sce)[, mad.str]  <- matrixStats::rowMads(values, na.rm = TRUE)

    return(sce)
}