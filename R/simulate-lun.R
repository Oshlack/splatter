simLun <- function(params = defaultParams(), verbose = TRUE, ...) {

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- mergeParams(params, defaultParams())
    params <- expandParams(params)

    # Set random seed
    seed <- getParams(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nGenes <- getParams(params, "nGenes")
    nCells <- getParams(params, "nCells")
    mean.shape <- getParams(params, "mean.shape")
    mean.rate <- getParams(params, "mean.rate")
    dispersion <- getParams(params, "bcv.common")

    if (verbose) {message("Simulating means...")}
    gene.means <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)

    if (verbose) {message("Simulating cell means...")}
    cell.facs <- 2 ^ rnorm(nCells, sd = 0.5)
    cell.means <- outer(gene.means, cell.facs, "*")

    if (verbose) {message("Simulating counts...")}
    counts <- matrix(rnbinom(nGenes * nCells, mu = cell.means,
                             size = 1 / dispersion),
                     nrow = nGenes, ncol = nCells)

    if (verbose) {message("Creating final SCESet...")}
    cell.names <- paste0("Cell", 1:nCells)
    gene.names <- paste0("Gene", 1:nGenes)
    rownames(counts) <- gene.names
    colnames(counts) <- cell.names

    phenos <- new("AnnotatedDataFrame",
                  data = data.frame(Cell = cell.names, CellFac = cell.facs))
    rownames(phenos) <- cell.names
    features <- new("AnnotatedDataFrame",
                    data = data.frame(Gene = gene.names, GeneMean = gene.means))
    rownames(features) <- gene.names
    sim <- newSCESet(countData = counts, phenoData = phenos,
                     featureData = features)

    rownames(cell.means) <- gene.names
    colnames(cell.means) <- cell.names
    assayData(sim)$CellMeans <- cell.means

    return(sim)
}

# mean.rate = 2
# mean.shape = 2