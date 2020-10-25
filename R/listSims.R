#' List simulations
#'
#' List all the simulations that are currently available in Splatter with a
#' brief description.
#'
#' @param print logical. Whether to print to the console.
#'
#' @return Invisibly returns a data.frame containing the information that is
#' displayed.
#' @examples
#' listSims()
#' @export
listSims <- function(print = TRUE) {

    sims <- list(c("Splat", "splat", "10.1186/s13059-017-1305-0",
                   "Oshlack/splatter",
                   "The Splat simulation generates means from a gamma
                   distribution, adjusts them for BCV and generates counts from
                   a gamma-poisson. Dropout and batch effects can be optionally
                   added.",
                   ""),
                 c("Splat Single", "splatSingle", "10.1186/s13059-017-1305-0",
                   "Oshlack/splatter",
                   "The Splat simulation with a single population.",
                   ""),
                 c("Splat Groups", "splatGroups", "10.1186/s13059-017-1305-0",
                   "Oshlack/splatter",
                   "The Splat simulation with multiple groups. Each group can
                   have it's own differential expression probability and
                   fold change distribution.",
                   ""),
                 c("Splat Paths", "splatPaths", "10.1186/s13059-017-1305-0",
                   "Oshlack/splatter",
                   "The Splat simulation with differentiation paths. Each
                   path can have it's own length, skew and probability. Genes
                   can change in non-linear ways.",
                   ""),
                 c("Kersplat", "kersplat", "",
                   "Oshlack/splatter",
                   "The Kersplat simulation extends the Splat model by adding a
                   gene network, more complex cell structure, doublets and
                   empty cells (Experimental).",
                   "scuttle, igraph"),
                 c("splatPop", "splatPop", "",
                   "Oshlack/splatter",
                   "The splatPop simulation enables splat simulations to be
                   generated for multiple individuals in a population,
                   accounting for correlation structure by simulating
                   expression quantitative trait loci (eQTL).",
                   "VariantAnnotation, preprocessCore"),
                 c("Simple", "simple", "10.1186/s13059-017-1305-0",
                   "Oshlack/splatter",
                   "A simple simulation with gamma means and negative binomial
                   counts.",
                   ""),
                 c("Lun", "lun", "10.1186/s13059-016-0947-7",
                   "MarioniLab/Deconvolution2016",
                   "Gamma distributed means and negative binomial counts. Cells
                   are given a size factor and differential expression can be
                   simulated with fixed fold changes.",
                   ""),
                 c("Lun 2", "lun2", "10.1093/biostatistics/kxw055",
                   "MarioniLab/PlateEffects2016",
                   "Negative binomial counts where the means and dispersions
                   have been sampled from a real dataset. The core feature of
                   the Lun 2 simulation is the addition of plate effects.
                   Differential expression can be added between two groups of
                   plates and optionally a zero-inflated negative-binomial can
                   be used.",
                   "scran, lme4, pscl, limSolve"),
                 c("scDD", "scDD", "10.1186/s13059-016-1077-y",
                   "kdkorthauer/scDD",
                   "The scDD simulation samples a given dataset and can
                   simulate differentially expressed and differentially
                   distributed genes between two conditions.",
                   "scDD"),
                 c("BASiCS", "BASiCS", "10.1371/journal.pcbi.1004333",
                   "catavallejos/BASiCS",
                   "The BASiCS simulation is based on a bayesian model used to
                   deconvolve biological and technical variation and
                   includes spike-ins and batch effects.",
                   "BASiCS"),
                 c("mfa", "mfa", "10.12688/wellcomeopenres.11087.1",
                   "kieranrcampbell/mfa",
                   "The mfa simulation produces a bifurcating pseudotime
                   trajectory. This can optionally include genes with transient
                   changes in expression and added dropout.",
                   "mfa"),
                 c("PhenoPath", "pheno", "10.1101/159913",
                   "kieranrcampbell/phenopath",
                   "The PhenoPath simulation produces a pseudotime trajectory
                   with different types of genes.",
                   "phenopath"),
                 c("ZINB-WaVE", "zinb", "10.1101/125112",
                   "drisso/zinbwave",
                   "The ZINB-WaVE simulation simulates counts from a
                   sophisticated zero-inflated negative-binomial distribution
                   including cell and gene-level covariates.",
                   "zinbwave"),
                 c("SparseDC", "sparseDC", "10.1093/nar/gkx1113",
                   "cran/SparseDC",
                   "The SparseDC simulation simulates a set of clusters
                   across two conditions, where some clusters may be present in
                   only one condition.",
                   "SparseDC")
                 )

    sims.table <- data.frame(Name         = rep(NA, length(sims)),
                             Prefix       = rep(NA, length(sims)),
                             DOI          = rep(NA, length(sims)),
                             GitHub       = rep(NA, length(sims)),
                             Description  = rep(NA, length(sims)),
                             Dependencies = rep(NA, length(sims)))

    for (idx in seq_along(sims)) {
        entry <- sims[[idx]]
        entry[5] <- gsub("[[:space:]]+", " ", entry[5])
        sims.table[idx, ] <- entry
    }

    if (print) {
        cat("Splatter currently contains", length(sims), "simulations", "\n\n")
        for (idx in seq_len(nrow(sims.table))) {
            sim <- as.character(sims.table[idx, ])
            cat(crayon::bold(sim[1]), crayon::yellow(paste0("(", sim[2], ")")),
                "\n")
            cat(crayon::bold("DOI:"), crayon::cyan(sim[3]), "\t",
                crayon::bold("GitHub:"), crayon::cyan(sim[4]), "\t",
                crayon::bold("Dependencies:"), crayon::cyan(sim[6]), "\n")
            cat(sim[5], "\n\n")
        }
    }

    invisible(sims.table)
}
