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

    sims <- list(c("Splat", "splat", "", "",
                   "The Splat simulation generates means from a gamma
                   distribution, adjusts them for BCV and generates counts from
                   a gamma-poisson. Dropout can be optionally added."),
                 c("Splat Single", "splatSingle", "", "",
                   "The Splat simulation with a single population."),
                 c("Splat Groups", "splatGroups", "", "",
                   "The Splat simulation with multiple groups. Each group can
                   have it's own differential expression probability and
                   fold change distribution."),
                 c("Splat Paths", "splatPaths", "", "",
                   "The Splat simulation with differentiation paths. Each
                   path can have it's own length, skew and probability. Genes
                   can change in non-linear ways."),
                 c("Simple", "simple", "", "",
                   "A simple simulation with gamma means and negative binomial
                   counts."),
                 c("Lun", "lun", "10.1186/s13059-016-0947-7",
                   "MarioniLab/Deconvolution2016",
                   "Gamma distributed means and negative binomial counts. Cells
                   are given a size factor and differential expression can be
                   simulated with fixed fold changes."),
                 c("Lun 2", "lun2", "10.1101/073973",
                   "MarioniLab/PlateEffects2016",
                   "Negative binomial counts where the means and dispersions
                   have been sampled from a real dataset. The core feature of
                   the Lun 2 simulation is the addition of plate effects.
                   Differential expression can be added between two groups of
                   plates and optionally a zero-inflated negative-binomial can
                   be used."),
                 c("scDD", "scDD", "10.1186/s13059-016-1077-y",
                   "kdkorthauer/scDD",
                   "The scDD simulation samples a given dataset and can
                   simulate differentially expressed and differentially
                   distributed genes between two conditions."),
                 c("BASiCS", "BASiCS", "10.1371/journal.pcbi.1004333",
                   "catavallejos/BASiCS",
                   "The BASiCS simulation is based on a bayesian model used to
                   deconvolve biological and technical variation and
                   includes spike-ins and batch effects."))

    sims.table <- data.frame(Name        = rep(NA, length(sims)),
                             Prefix      = rep(NA, length(sims)),
                             DOI         = rep(NA, length(sims)),
                             Github      = rep(NA, length(sims)),
                             Description = rep(NA, length(sims)))

    for (idx in seq_along(sims)) {
        entry <- sims[[idx]]
        entry[5] <- gsub("[[:space:]]+", " ", entry[5])
        sims.table[idx, ] <- entry
    }

    if (print) {
        cat("Splatter currently contains", length(sims), "simulations", "\n\n")
        for (idx in seq_len(nrow(sims.table))) {
            sim <- as.character(sims.table[idx, ])
            cat(sim[1], paste0("(", sim[2], ")"), "\n")
            cat("DOI:", sim[3], "\t", "Github:", sim[4], "\n")
            cat(sim[5], "\n\n")
        }
    }

    invisible(sims.table)
}
