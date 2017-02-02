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

#' Add gene lengths
#'
#' Add gene lengths to an SCESet object
#'
#' @param sce SCESet to add gene lengths to.
#' @param method Method to use for creating lengths.
#' @param loc Location parameter for the generate method.
#' @param scale Scale parameter for the generate method.
#' @param lengths Vector of lengths for the sample method.
#'
#' @details
#' This function adds simulated gene lengths to the \code{fData} slot of an
#' \code{SCESet} object that can be used for calculating length normalised
#' expression values such as TPM or FPKM. The \code{generate} simulates lengths
#' using a (rounded) log-normal distribution, with the default \code{loc} and
#' \code{scale} parameters based on human coding genes. Alternatively the
#' \code{sample} method can be used which randomly samples lengths (with
#' replacement) from a supplied vector.
#'
#' @return SCESet with added gene lengths
#' @examples
#' # Default generate method
#' sce <- simpleSimulate()
#' sce <- addGeneLengths(sce)
#' head(fData(sce))
#' # Sample method (human coding genes)
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(GenomicFeatures)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' tx.lens <- transcriptLengths(txdb, with.cds_len = TRUE)
#' tx.lens <- tx.lens[tx.lens$cds_len > 0, ]
#' gene.lens <- max(splitAsList(tx.lens$tx_len, tx.lens$gene_id))
#' sce <- addGeneLengths(sce, method = "sample", lengths = gene.lens)
#' }
#' @export
#' @importFrom stats rlnorm
addGeneLengths <- function(sce, method = c("generate", "sample"), loc = 7.9,
                           scale = 0.7, lengths = NULL) {

    method <- match.arg(method)
    checkmate::assertClass(sce, "SCESet")
    checkmate::assertNumber(loc)
    checkmate::assertNumber(scale, lower = 0)
    checkmate::assertNumeric(lengths, lower = 0, null.ok = TRUE)

    switch(method,
           generate = {
               sim.lengths <- rlnorm(nrow(sce), meanlog = loc, sdlog = scale)
               sim.lengths <- round(sim.lengths)
           },
           sample = {
               if (is.null(lengths)) {
                   stop("Lengths must be supplied to use the sample method.")
               } else {
                   sim.lengths <- sample(lengths, nrow(sce), replace = TRUE)
               }
           }
    )

    fData(sce)$Length <- sim.lengths

    return(sce)
}

#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#tx.lens <- GenomicFeatures::transcriptLengths(txdb, with.cds_len = TRUE)
#tx.lens <- tx.lens[tx.lens$cds_len > 0, ]
#gene.lens <- max(IRanges::splitAsList(tx.lens$tx_len, tx.lens$gene_id))
#lens <- rlnorm(length(gene.lens), meanlog = 7.9, sdlog = 0.7)
