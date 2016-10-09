# setup metadata
# Means
# Groups
# Paths
# Lib size
# Base Means
# BCV
# Means
# Counts
# Dropout
# Add metadata
#
# Add length
# Median outliers
# set.seed
# group DE



#' FUNCTION TITLE
#'
#' FUNCTION DESCRIPTION
#'
#' @param params DESCRIPTION.
#' @param method DESCRIPTION.
#' @param add.assay DESCRIPTION.
#' @param verbose DESCRIPTION.
#' @param ... DESCRIPTION.
#'
#' @return RETURN DESCRIPTION
#' @examples
#' # ADD EXAMPLES HERE
#' @importFrom BioBase pData fData
#' @importFrom scater newSCESet counts
splat <- function(params = defaultParams(), method = c("groups", "paths"),
                  add.assay = TRUE, verbose = TRUE, ...) {

    method <- match.arg(method)

    params <- setParams(params, ...)
    params <- mergeParams(params, defaultParams())

    # Get the parameters we are going to use
    n.cells <- getParams(params, "nCells")
    n.genes <- getParams(params, "nGenes")
    n.groups <- getParams(params, "nGroups")
    group.cells <- getParams(params, "groupCells")

    # Set up name vectors
    cell.names <- paste0("Cell", 1:n.cells)
    gene.names <- paste0("Gene", 1:n.genes)
    group.names <- paste0("Group", 1:n.groups)

    # Create SCESet with dummy counts to store simulation
    dummy.counts <- matrix(1, ncol = n.cells, nrow = n.genes)
    rownames(dummy.counts) <- gene.names
    colnames(dummy.counts) <- cell.names
    phenos <- new("AnnotatedDataFrame", data = data.frame(Cell = cell.names))
    rownames(phenos) <- cell.names
    features <- new("AnnotatedDataFrame", data = data.frame(Gene = gene.names))
    rownames(features) <- gene.names
    sim <- newSCESet(countData = dummy.counts, phenoData = phenos,
                     featureData = features)

    # Make groups vector which is the index of param$groupCells repeated
    # params$groupCells[index] times
    groups <- lapply(1:n.groups, function(i, g) {rep(i, g[i])},
                     g = group.cells)
    groups <- unlist(groups)
    pData(sim)$Group <- group.names[groups]

    # Create new SCESet to make sure values are calculated correctly
    sce <- newSCESet(countData = counts(sim),
                     phenoData = new("AnnotatedDataFrame", data = pData(sim)),
                     featureData = new("AnnotatedDataFrame", data = fData(sim)))

    return(sce)
}