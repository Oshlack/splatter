estimateSimpleParams <- function(data, params = newSimpleParams()) {
    UseMethod("estimateSimpleParams")
}

#' @rdname simpleEstimate
#' @export
estimateSimpleParams.SCESet <- function(data, params = newSimpleParams()) {
    counts <- scater::counts(x)
    estimateSimpleParams(counts, params)
}

#' @rdname simpleEstimate
#' @importFrom stats median
#' @export
estimateSimpleParams.matrix <- function(data, params = newSimpleParams()) {

    checkmate::assertClass(params, "SimpleParams")

    # Normalise for library size and remove all zero genes
    lib.sizes <- colSums(data)
    lib.med <- median(lib.sizes)
    norm.counts <- t(t(data) / lib.sizes * lib.med)
    norm.counts <- norm.counts[rowSums(norm.counts > 0) > 1, ]

    means <- rowMeans(norm.counts)

    means.fit <- fitdistrplus::fitdist(means, "gamma", method = "mme")

    params <- setParams(params, nGenes = nrow(data), nCells = nrow(data),
                        mean.shape = unname(means.fit$estimate["shape"]),
                        mean.rate = unname(means.fit$estimate["rate"]))

    return(params)
}