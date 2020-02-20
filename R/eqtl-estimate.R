#' Estimate eQTL simulation parameters
#'
#' Estimate simulation parameters for the eQTL population simulation from
#' real data. See the individual estimation functions for more details on
#' how this is done. 
#'
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' @param all.pairs Txt file with all eQTL pairs from a bulk eQTL analysis. Must
#'        include 3 columns: 'gene_id', 'pval_nominal', and 'slope'.
#' @param gene.means Dataframe of real gene means across a population, where 
#'        each row is a gene and each column is an individual in the population.
#'
#' @seealso
#' \code{\link{eQTLEstES}},  \code{\link{eQTLEstMeanCV}}
#'
#' @return eQTLParams object containing the estimated parameters.
#'
#' @examples
#' # Load example data
#' data(ex_means)
#' data(ex_pairs)
#' eQTLparams <- eQTLEstimate()
#' eQTLparams
#' 
#' @export
eQTLEstimate <- function(eQTLparams = neweQTLParams(),
                         all.pairs = ex_pairs,
                         gene.means = ex_means){
    
    checkmate::assertClass(eQTLparams, "eQTLParams")
    
    # Get parameters for eQTL Effect Size distribution
    eQTLparams <- eQTLEstES(all.pairs, eQTLparams)
        
    # Get parameters for population wide gene mean and variance distributions
    eQTLparams <- eQTLEstMeanCV(gene.means, eQTLparams)

    return(eQTLparams)
}

#' Estimate eQTL Effect Size parameters
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate eQTL (eSNP-eGene) effect sizes.
#'
#' @param all.pairs Txt file with all eQTL pairs from a bulk eQTL analysis. Must
#'        include 3 columns: 'gene_id', 'pval_nominal', and 'slope'. 
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#'
#' @details
#' Parameters for the gamma distribution are estimated by fitting the top eSNP-
#' eGene pair effect sizes using \code{\link[fitdistrplus]{fitdist}}. The 
#' 'maximum goodness-of-fit estimation' method is used to minimise the 
#' Cramer-von Mises distance. This can fail in some situations, in which case 
#' the 'method of moments estimation' method is used instead. 
#'
#' @return eQTLparams object with estimated values.
#' @importFrom data.table data.table .I
eQTLEstES <- function(all.pairs, eQTLparams) {

    # Select top eSNP for each gene (i.e. lowest p.value)
    all.pairs <- data.table::data.table(all.pairs)
    pairs_top <- all.pairs[all.pairs[, .I[which.min(pval_nominal)], 
                                     by='gene_id']$V1]

    # Fit absolute value of effect sizes to gamma distribution
    e.sizes <- abs(pairs_top$slope)
    fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mge", gof = "CvM")
    
    if (fit$convergence > 0) {
        warning("Fitting effect sizes using the Goodness of Fit method failed,",
                " using the Method of Moments instead")
        fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mme")
    }
    
    eQTLparams <- setParams(eQTLparams, 
                            eqtlES.shape = unname(fit$estimate["shape"]),
                            eqtlES.rate = unname(fit$estimate["rate"]))

    return(eQTLparams)
}

#' Estimate gene mean and gene mean variance parameters
#'
#' The Shapiro-Wilks test is used to determine if the library sizes are
#' normally distributed. If so a normal distribution is fitted to the library
#' sizes, if not (most cases) a log-normal distribution is fitted and the
#' estimated parameters are added to the params object. See
#' \code{\link[fitdistrplus]{fitdist}} for details on the fitting.
#'
#' @param eQTLparams eQTLParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{eQTLParams}} for details.
#' @param gene.means Dataframe of real gene means across a population, where 
#'        each row is a gene and each column is an individual in the population.
#'
#' @details
#' Parameters for the mean gamma distribution are estimated by fitting the mean
#' (across the population) expression of genes that meet the criteria (<50% of
#' samples have exp <0.1) and parameters for the cv gamma distribution are 
#' estimated for each bin of mean expression using the cv of expression across
#' the population for genes in that bin. Both are fit using 
#' \code{\link[fitdistrplus]{fitdist}}. The "Nelder-Mead" method is used to fit
#' the mean gamma distribution and the 'maximum goodness-of-fit estimation' 
#' method is used to minimise the Cramer-von Mises distance for the CV 
#' distribution.
#' 
#' @return eQTLParams object with estimated values.
#' @importFrom stats quantile
eQTLEstMeanCV <- function(gene.means, eQTLparams) {
    
    # Remove genes with low variance/low means
    abv_thr <- data.frame(perc = (rowSums(gene.means >= 0.1)/ncol(gene.means)))
    genes <- row.names(subset(abv_thr, perc > 0.5))
    gene.means <- gene.means[genes, ]
    
    # Calculate mean expression parameters
    means <- rowMeans(gene.means)
    mfit <- fitdistrplus::fitdist(means, "gamma", optim.method="Nelder-Mead")
    
    # Calculate CV parameters for genes based on 10 expresion mean bins
    print(eQTLparams)
    nbins <- getParam(eQTLparams, "bulkcv.bins")
    bins <- split(means, cut(means, quantile(means,(0:nbins)/nbins), include.lowest=T))
    cvparams <- data.frame(start = character(), end = character(),
                           shape = character(), rate = character(), 
                           stringsAsFactors=FALSE)
    for(b in names(bins)){
        stst <- unlist(strsplit(gsub("\\)|\\(|\\[|\\]", "", b), ','))
        b_genes <- names(unlist(bins[b], use.names = T))
        b_genes <- gsub(paste0(b, '.'), '', b_genes, fixed=T)
        b_gene.means <- gene.means[b_genes, ]
        cv <- apply(b_gene.means, 1, co.var)
        cv[is.na(cv)] <- 0
        cvfit <- fitdistrplus::fitdist(cv, "gamma", method = "mge", gof = "CvM")
        cvparams <- rbind(cvparams, 
                          list(start= as.numeric(stst[1]),
                               end= as.numeric(stst[2]), 
                               shape=cvfit$estimate['shape'],
                               rate=cvfit$estimate['rate']),
                          stringsAsFactors=FALSE)
    }
    cvparams[1, 'start'] <- 0
    cvparams[nrow(cvparams), 'end'] <- 1e100
    eQTLparams <- setParams(eQTLparams, 
                             bulkmean.shape = unname(mfit$estimate["shape"]),
                             bulkmean.rate = unname(mfit$estimate["rate"]),
                             bulkcv.param = cvparams)

    return(eQTLparams)
}

