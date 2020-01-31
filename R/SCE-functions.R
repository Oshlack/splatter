#' Add feature statistics
#'
#' Add additional feature statistics to a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment to add feature statistics to.
#' @param value the expression value to calculate statistics for. Options are
#'        "counts", "cpm", "tpm" or "fpkm". The values need to exist in the
#'        given SingleCellExperiment.
#' @param log logical. Whether to take log2 before calculating statistics.
#' @param offset offset to add to avoid taking log of zero.
#' @param no.zeros logical. Whether to remove all zeros from each feature before
#'        calculating statistics.
#'
#' @details
#' Currently adds the following statistics: mean, variance, coefficient of
#' variation, median and median absolute deviation. Statistics are added to
#' the \code{\link{rowData}} slot and are named \code{Stat[Log]Value[No0]} where
#' \code{Log} and \code{No0} are added if those arguments are true.
#' UpperCamelCase is used to differentiate these columns from those added by
#' analysis packages.
#'
#' @return SingleCellExperiment with additional feature statistics
#'
#' @importFrom SummarizedExperiment rowData rowData<-
addFeatureStats <- function(sce, value = c("counts", "cpm", "tpm", "fpkm"),
                            log = FALSE, offset = 1, no.zeros = FALSE) {

    checkmate::assertClass(sce, "SingleCellExperiment")
    checkmate::assertLogical(log)
    checkmate::assertNumber(offset, lower = 0)
    checkmate::assertLogical(no.zeros)
    value <- match.arg(value)

    switch(value,
           counts = {
               values = BiocGenerics::counts(sce)
               suffix <- "Counts"
           },
           cpm = {
               values = SingleCellExperiment::cpm(sce)
               suffix <- "CPM"
           },
           tpm = {
               values = SingleCellExperiment::tpm(sce)
               suffix <- "TPM"
           },
           fpkm = {
               values = SummarizedExperiment::assays(sce)$fpkm
               suffix <- "FPKM"
           }
    )

    if (no.zeros) {
        values[values == 0] <- NA
        suffix = paste0(suffix, "No0")
    }

    if (log) {
        values = log2(values + offset)
        suffix = paste0("Log", suffix)
    }

    mean.str <- paste0("Mean", suffix)
    var.str  <- paste0("Var",  suffix)
    cv.str   <- paste0("CV",   suffix)
    med.str  <- paste0("Med",  suffix)
    mad.str  <- paste0("MAD",  suffix)

    rowData(sce)[, mean.str] <- rowMeans(values, na.rm = TRUE)
    rowData(sce)[, var.str]  <- matrixStats::rowVars(values, na.rm = TRUE)
    rowData(sce)[, cv.str]   <- sqrt(rowData(sce)[, var.str]) /
        rowData(sce)[, mean.str]
    rowData(sce)[, med.str]  <- matrixStats::rowMedians(values, na.rm = TRUE)
    rowData(sce)[, mad.str]  <- matrixStats::rowMads(values, na.rm = TRUE)

    return(sce)
}

#' Add gene lengths
#'
#' Add gene lengths to an SingleCellExperiment object
#'
#' @param sce SingleCellExperiment to add gene lengths to.
#' @param method Method to use for creating lengths.
#' @param loc Location parameter for the generate method.
#' @param scale Scale parameter for the generate method.
#' @param lengths Vector of lengths for the sample method.
#'
#' @details
#' This function adds simulated gene lengths to the
#' \code{\link{rowData}} slot of a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object that can be
#' used for calculating length normalised expression values such as TPM or FPKM.
#' The \code{generate} method simulates lengths using a (rounded) log-normal
#' distribution, with the default \code{loc} and \code{scale} parameters based
#' on human protein-coding genes. Alternatively the \code{sample} method can be
#' used which randomly samples lengths (with replacement) from a supplied
#' vector.
#'
#' @return SingleCellExperiment with added gene lengths
#' @examples
#' # Default generate method
#' sce <- simpleSimulate()
#' sce <- addGeneLengths(sce)
#' head(rowData(sce))
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
    checkmate::assertClass(sce, "SingleCellExperiment")
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

    rowData(sce)$Length <- sim.lengths

    return(sce)
}

#' Get counts
#'
#' Get counts matrix from a SingleCellExperiment object. If counts is missing
#' a warning is issued and the first assay is returned.
#'
#' @param sce SingleCellExperiment object
#'
#' @return Counts matrix
getCounts <- function(sce) {

    checkmate::assertClass(sce, "SingleCellExperiment")

    if ("counts" %in% SummarizedExperiment::assayNames(sce)) {
        counts <- SingleCellExperiment::counts(sce)
    } else {
        warning("counts assay is missing, using the first assay instead")
        counts <- SummarizedExperiment::assay(sce)
    }

    return(counts)
}
