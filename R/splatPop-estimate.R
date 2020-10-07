#' Estimate population/eQTL simulation parameters
#'
#' Estimate simulation parameters for the eQTL population simulation from
#' real data. See the individual estimation functions for more details on
#' how this is done. 
#'
#' @param params splatPopParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details.
#' @param counts either a counts matrix or a SingleCellExperiment object
#'        containing count data to estimate parameters from.
#' @param means Dataframe of real gene means across a population, where 
#'        each row is a gene and each column is an individual in the population.
#' @param eqtl Txt file with all or top eQTL pairs from a real eQTL analysis.
#'         Must include columns: 'gene_id', 'pval_nominal', and 'slope'.
#' 
#' @seealso
#' \code{\link{splatPopEstimateEffectSize}},
#' \code{\link{splatPopEstimateMeanCV}}
#'
#' @return splatPopParams object containing the estimated parameters.
#' 
#' @export
#' 
splatPopEstimate <- function(params = newSplatPopParams(), 
                             counts = "skip",
                             means = "skip",
                             eqtl = "skip"){
    
    checkmate::assertClass(params, "splatPopParams")

    # Estimate single-cell parameters using base splatEstimate function
    # THIS DOESN'T WORK BECAUSE IT RETURNS A SPLATPARAM NOT A SPLATPOPPARAM...
    #if(!is.character(counts)){params <- splatEstimate(counts)}
    
    # Get parameters for eQTL Effect Size distribution
    if(!is.character(eqtl)){params <- splatPopEstimateEffectSize(params, eqtl)}

    # Get parameters for population wide gene mean and variance distributions
    if(!is.character(means)){params <- splatPopEstimateMeanCV(params, means)}
    
    return(params)
}

#' Estimate eQTL Effect Size parameters
#'
#' Estimate rate and shape parameters for the gamma distribution used to
#' simulate eQTL (eSNP-eGene) effect sizes.
#'
#' @param eqtl Txt file with all or top eQTL pairs from a real eQTL analysis. 
#'         Must include columns: gene_id, pval_nominal, and slope. 
#' @param params splatPopParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details.
#'
#' @details
#' Parameters for the gamma distribution are estimated by fitting the top eSNP-
#' eGene pair effect sizes using \code{\link[fitdistrplus]{fitdist}}. The 
#' maximum goodness-of-fit estimation method is used to minimise the 
#' Cramer-von Mises distance. This can fail in some situations, in which case 
#' the method of moments estimation method is used instead. 
#'
#' @return params object with estimated values.
#' 
splatPopEstimateEffectSize <- function(params, eqtl) {

    # Test input eSNP-eGene pairs
    if (!("gene_id" %in% names(eqtl) & 
          "pval_nominal" %in% names(eqtl) & 
          "slope" %in% names(eqtl))) {
        stop("Incorrect format for eqtl data.")}
    
    # Select top eSNP for each gene (i.e. lowest p.value)
    eqtl_top <- eqtl[order(eqtl$gene_id, eqtl$pval_nominal), ]
    eqtl_top <- eqtl_top[!duplicated(eqtl_top$gene_id), ]

    # Fit absolute value of effect sizes to gamma distribution
    e.sizes <- abs(eqtl_top$slope)
    fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mge", gof = "CvM")
    if (fit$convergence > 0) {
        warning("Fitting effect sizes using the Goodness of Fit method failed,",
                " using the Method of Moments instead")
        fit <- fitdistrplus::fitdist(e.sizes, "gamma", method = "mme")
    }

    params <- setParams(params, eqtl.ES.shape = unname(fit$estimate["shape"]),
                        eqtl.ES.rate = unname(fit$estimate["rate"]))

    return(params)
}

#' Estimate gene mean and gene mean variance parameters
#'
#' The Shapiro-Wilks test is used to determine if the library sizes are
#' normally distributed. If so a normal distribution is fitted to the library
#' sizes, if not (most cases) a log-normal distribution is fitted and the
#' estimated parameters are added to the params object. See
#' \code{\link[fitdistrplus]{fitdist}} for details on the fitting.
#'
#' @param params splatPopParams object containing parameters for the 
#'         simulation of the mean expression levels for the population.
#'        See \code{\link{splatPopParams}} for details.
#' @param means Dataframe of real gene means across a population, where 
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
#' @return params object with estimated values.
#' @importFrom stats quantile
#' @importFrom grDevices boxplot.stats
#' @importFrom matrixStats rowMedians
#' 
splatPopEstimateMeanCV <- function(params, means) {
    
    # Test input gene means
    if ((anyNA(means) | !(validObject(rowSums(means))))) {
        stop("Incorrect format or NAs present in gene.means. See example data.")
    }
    
    # Remove genes with low variance/low means
    abv_thr <- data.frame(perc = (rowSums(means >= 0.1)/ncol(means)),
                          gene_id = row.names(means))
    
    genes <- row.names(abv_thr[abv_thr$perc > 0.5, ])
    means <- means[genes, ]

    # Calculate mean expression parameters
    row_means <- rowMedians(as.matrix(means))
    names(row_means) <- row.names(means)
    mfit <- fitdistrplus::fitdist(row_means, "gamma",
                                  optim.method="Nelder-Mead")
    
    # Calculate CV parameters for genes based on 10 expresion mean bins
    nbins <- getParam(params, "pop.cv.bins")
    bins <- split(row_means, cut(row_means, quantile(row_means,(0:nbins)/nbins), 
                             include.lowest=TRUE))
    cvparams <- data.frame(start = character(), end = character(),
                           shape = character(), rate = character(), 
                           stringsAsFactors=FALSE)

    for(b in names(bins)){
        stst <- unlist(strsplit(gsub("\\)|\\(|\\[|\\]", "", b), ","))
        b_genes <- names(unlist(bins[b], use.names = T))
        b_genes <- gsub(paste0(b, "."), "", b_genes, fixed=T)
        b_gene.means <- means[b_genes, ]
        
        cv <- apply(b_gene.means, 1, co.var)
        cv[is.na(cv)] <- 0
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
    params <- setParams(params, 
                        pop.mean.shape = unname(mfit$estimate["shape"]),
                        pop.mean.rate = unname(mfit$estimate["rate"]),
                        pop.cv.param = cvparams)
    
    return(params)
}
