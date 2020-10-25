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

    if (is(values, "dgCMatrix")) {
        values <- as.matrix(values)
    }

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

#' Minimise SCE
#'
#' Reduce the size of a SingleCellExperiment object by unneeded information.
#'
#' @param sce SingleCellExperiment object
#' @param rowData.keep Either TRUE (keep all rowData columns), FALSE (remove all
#' rowData columns) or a character vector with the names of the rowData columns
#' to keep
#' @param colData.keep Either TRUE (keep all colData columns), FALSE (remove all
#' colData columns) or a character vector with the names of the colData columns
#' to keep
#' @param metadata.keep Either TRUE (keep all metadata), FALSE (remove all
#' metadata) or a character vector with the names of the metadata items to keep
#' @param assays.keep Either TRUE (keep all assays), FALSE (remove all
#' assays) or a character vector with the names of the assays to keep
#' @param sparsify Whether to convert assay matrices to sparse format. Either
#' "all", "none" or "auto" (default) to only convert those matrices that will
#' result in a size reduction
#' @param verbose Whether to print status messages
#'
#' @return SingleCellExperiment object
#'
#' @examples
#'
#' sce <- splatSimulate(verbose = FALSE)
#' sce.min <- minimiseSCE(sce, verbose = FALSE)
#' object.size(sce)
#' object.size(sce.min)
#'
#' @export
#' @importFrom utils object.size
minimiseSCE <- function(sce, rowData.keep = FALSE, colData.keep = FALSE,
                        metadata.keep = FALSE, assays.keep = "counts",
                        sparsify = c("auto", "all", "none"), verbose = TRUE) {

    sparsify <- match.arg(sparsify)

    if (verbose) {
        start.size <- object.size(sce)
        message("Minimising SingleCellExperiment...")
        message("Original size: ", format(start.size, unit = "auto"))
    }

    if (isFALSE(rowData.keep)) {
        if (verbose) {message("Removing all rowData columns")}
        SummarizedExperiment::rowData(sce) <- NULL
    } else if (is.character(rowData.keep)) {
        rowData <- SummarizedExperiment::rowData(sce)
        keep <- colnames(rowData) %in% rowData.keep
        if (verbose) {
            message("Keeping ", sum(keep), " rowData columns: ",
                    paste(colnames(rowData)[keep], collapse = ", "))
            message("Removing ", sum(!keep), " rowData columns: ",
                    paste(colnames(rowData)[!keep], collapse = ", "))
        }
        SummarizedExperiment::rowData(sce) <- rowData[, keep, drop = FALSE]
    }

    if (isFALSE(colData.keep)) {
        if (verbose) {message("Removing all colData columns")}
        SummarizedExperiment::colData(sce) <- NULL
    } else if (is.character(colData.keep)) {
        colData <- SummarizedExperiment::colData(sce)
        keep <- colnames(colData) %in% colData.keep
        if (verbose) {
            message("Keeping ", sum(keep), " colData columns: ",
                    paste(colnames(colData)[keep], collapse = ", "))
            message("Removing ", sum(!keep), " colData columns: ",
                    paste(colnames(colData)[!keep], collapse = ", "))
        }
        SummarizedExperiment::colData(sce) <- colData[, keep, drop = FALSE]
    }

    if (isFALSE(metadata.keep)) {
        if (verbose) {message("Removing all metadata items")}
        S4Vectors::metadata(sce) <- list()
    } else if (is.character(metadata.keep)) {
        metadata <- S4Vectors::metadata(sce)
        keep <- names(metadata) %in% metadata.keep
        if (verbose) {
            message("Keeping ", sum(keep), " metadata items: ",
                    paste(names(metadata)[keep], collapse = ", "))
            message("Removing ", sum(!keep), " metadata items: ",
                    paste(names(metadata)[!keep], collapse = ", "))
        }
        S4Vectors::metadata(sce) <- metadata[keep]
    }

    if (isFALSE(assays.keep)) {
        if (verbose) {message("Removing all assays")}
        SummarizedExperiment::assays(sce) <- list()
    } else if (is.character(assays.keep)) {
        assays <- SummarizedExperiment::assays(sce)
        keep <- names(assays) %in% assays.keep
        if (verbose) {
            message("Keeping ", sum(keep), " assays: ",
                    paste(names(assays)[keep], collapse = ", "))
            message("Removing ", sum(!keep), " assays: ",
                    paste(names(assays)[!keep], collapse = ", "))
        }
        SummarizedExperiment::assays(sce) <- assays[keep]
    }

    if (sparsify != "none") {
        if (verbose) {message("Sparsifying assays...")}
        SummarizedExperiment::assays(sce) <- sparsifyMatrices(
            SummarizedExperiment::assays(sce), auto = sparsify == "auto",
            verbose = verbose)
    }

    if (verbose) {
        final.size <- object.size(sce)
        message("Minimised size: ", format(final.size, unit = "auto"),
                " (", round(final.size / start.size * 100), "% of original)")
    }

    return(sce)
}

#' Sparsify matrices
#'
#' Convert a list of matrices to sparse matrices.
#'
#' @param matrix.list List of matrices
#' @param auto Whether to automatically choose which matrices to convert based
#' on how big the size reduction will be
#' @param threshold Threshold for automatically selecting matrices to convert,
#' any matrix with an estimated sparse size less than this proportion of the
#' original size will be converted
#' @param verbose Whether to print status messages
#'
#' @return List of converted matrices
sparsifyMatrices <- function(matrix.list, auto = TRUE, threshold = 0.95,
                             verbose = TRUE) {

    # If not auto, just do it and return
    if (!auto) {

        if (verbose) {message("Converting all matrices to sparse format")}
        matrix.list <- lapply(matrix.list, as, Class = "dgCMatrix")

        return(matrix.list)
    }

    if (verbose) {
        message("Automatically converting to sparse matrices, ",
                "threshold = ", threshold)
    }
    for (mat.name in names(matrix.list)) {
        mat <- matrix.list[[mat.name]]
        prop.zero <- sum(mat == 0) / length(mat)
        if (is.integer(mat)) {
            size.factor <- 3 - 3 * prop.zero
        } else if (is.numeric(mat)) {
            size.factor <- 1.5 - 1.5 * prop.zero
        } else if (is(mat, "dgCMatrix")) {
            if (verbose) {
                message("Skipping '", mat.name,
                        "' as it is already a dgCMatrix")
            }
            next
        } else {
            warning("matrix '", mat.name, "' is class '", class(mat),
                    "', unable to estimate size reduction factor")
            size.factor <- NA
        }
        if (is.na(size.factor) | size.factor < threshold) {
            if (verbose) {
                message("Converting '", mat.name, "' to sparse matrix: ",
                        "estimated sparse size ", round(size.factor, 2),
                        " * dense matrix")
            }
            matrix.list[[mat.name]] <- as(mat, "dgCMatrix")
        } else {
            if (verbose) {
                message("Skipping '", mat.name, "': estimated sparse size ",
                        round(size.factor, 2), " * dense matrix")
            }
        }
    }

    return(matrix.list)
}
