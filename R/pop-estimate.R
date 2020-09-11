#' Estimate population/eQTL simulation parameters
#'
#' Estimate simulation parameters for the eQTL population simulation from
#' real data. See the individual estimation functions for more details on
#' how this is done. 
#'
#' @param gene.means Dataframe of real gene means across a population, where 
#'        each row is a gene and each column is an individual in the population.
#' @param eqtl Txt file with all or top eQTL pairs from a real eQTL analysis.
#'         Must include columns: 'gene_id', 'pval_nominal', and 'slope'.
#' @param popParams popParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{popParams}} for details.
#' 
#'
#' @seealso
#' \code{\link{popEstimate.ES}},  \code{\link{popEstimate.MeanCV}}
#'
#' @return popParams object containing the estimated parameters.
#' 
#' @export
#' 
popEstimate <- function(gene.means, eqtl, popParams = newPopParams()){
    
    checkmate::assertClass(popParams, "popParams")
    
    # Get parameters for eQTL Effect Size distribution
    popParams <- popEstimate.ES(eqtl, popParams)
        
    # Get parameters for population wide gene mean and variance distributions
    popParams <- popEstimate.MeanCV(gene.means, popParams)

    return(popParams)
}

#' Estimate eQTL Effect Size parameters
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate eQTL (eSNP-eGene) effect sizes.
#'
#' @param eqtl Txt file with all or top eQTL pairs from a real eQTL analysis. 
#'         Must include columns: gene_id, pval_nominal, and slope. 
#' @param popParams popParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{popParams}} for details.
#'
#' @details
#' Parameters for the gamma distribution are estimated by fitting the top eSNP-
#' eGene pair effect sizes using \code{\link[fitdistrplus]{fitdist}}. The 
#' maximum goodness-of-fit estimation method is used to minimise the 
#' Cramer-von Mises distance. This can fail in some situations, in which case 
#' the method of moments estimation method is used instead. 
#'
#' @return popParams object with estimated values.
#' @importFrom dplyr group_by filter "%>%"
#' 
popEstimate.ES <- function(eqtl, popParams) {

    # Test input eSNP-eGene pairs
    if (!("gene_id" %in% names(eqtl) &
          "pval_nominal" %in% names(eqtl))) {
        stop("Incorrect format for eqtl data.")
    }
    
    # Select top eSNP for each gene (i.e. lowest p.value)
    pairs_top <- eqtl %>% group_by(gene_id) %>%
        filter(pval_nominal == min(pval_nominal))

    # Fit absolute value of effect sizes to gamma distribution
    e.sizes <- abs(pairs_top$slope)
    fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mge", gof = "CvM")
    
    if (fit$convergence > 0) {
        warning("Fitting effect sizes using the Goodness of Fit method failed,",
                " using the Method of Moments instead")
        fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mme")
    }
    
    popParams <- setParams(popParams, 
                            eqtl.ES.shape = unname(fit$estimate["shape"]),
                            eqtl.ES.rate = unname(fit$estimate["rate"]))

    return(popParams)
}

#' Estimate gene mean and gene mean variance parameters
#'
#' The Shapiro-Wilks test is used to determine if the library sizes are
#' normally distributed. If so a normal distribution is fitted to the library
#' sizes, if not (most cases) a log-normal distribution is fitted and the
#' estimated parameters are added to the params object. See
#' \code{\link[fitdistrplus]{fitdist}} for details on the fitting.
#'
#' @param popParams popParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{popParams}} for details.
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
#' the mean gamma distribution and the maximum goodness-of-fit estimation 
#' method is used to minimise the Cramer-von Mises distance for the CV 
#' distribution.
#' 
#' @return popParams object with estimated values.
#' @importFrom stats quantile
#' @importFrom grDevices boxplot.stats
#' @importFrom matrixStats rowMedians
#' 
popEstimate.MeanCV <- function(gene.means, popParams) {
    
    # Test input gene means
    if ((anyNA(gene.means) | !(validObject(rowSums(gene.means))))) {
        stop("Incorrect format or NAs present in gene.means. See example data.")
    }
    
    # Remove genes with low variance/low means
    abv_thr <- data.frame(perc = (rowSums(gene.means >= 0.1)/ncol(gene.means)),
                          gene_id = row.names(gene.means))
    
    genes <- row.names(abv_thr[abv_thr$perc > 0.5, ])
    gene.means <- gene.means[genes, ]
    
    # Calculate mean expression parameters
    means <- rowMedians(as.matrix(gene.means))
    names(means) <- row.names(gene.means)
    mfit <- fitdistrplus::fitdist(means, "gamma", optim.method="Nelder-Mead")
    
    # Calculate CV parameters for genes based on 10 expresion mean bins
    nbins <- getParam(popParams, "pop.cv.bins")
    bins <- split(means, cut(means, quantile(means,(0:nbins)/nbins), 
                             include.lowest=TRUE))
    cvparams <- data.frame(start = character(), end = character(),
                           shape = character(), rate = character(), 
                           stringsAsFactors=FALSE)
    for(b in names(bins)){
        stst <- unlist(strsplit(gsub("\\)|\\(|\\[|\\]", "", b), ","))
        b_genes <- names(unlist(bins[b], use.names = T))
        b_genes <- gsub(paste0(b, "."), "", b_genes, fixed=T)
        b_gene.means <- gene.means[b_genes, ]
        
        cv <- apply(b_gene.means, 1, co.var)
        cv[is.na(cv)] <- 0
        # outliers <- boxplot(cv)$out
        # cv <- cv[-outliers]
        cv <- cv[!cv %in% boxplot.stats(cv)$out]
        cvfit <- fitdistrplus::fitdist(cv, "gamma", method = "mge", gof = "CvM")
        cvparams <- rbind(cvparams, 
                          list(start= as.numeric(stst[1]),
                               end= as.numeric(stst[2]), 
                               shape=cvfit$estimate["shape"],
                               rate=cvfit$estimate["rate"]),
                          stringsAsFactors=FALSE)
    }
    cvparams[1, "start"] <- 0
    cvparams[nrow(cvparams), "end"] <- 1e100
    popParams <- setParams(popParams, 
                            pop.mean.shape = unname(mfit$estimate["shape"]),
                            pop.mean.rate = unname(mfit$estimate["rate"]),
                            pop.cv.param = cvparams)
    
    return(popParams)
}
