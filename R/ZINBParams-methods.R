#' @rdname newParams
#' @importFrom methods new
#' @export
newZINBParams <- function(...) {

    if (!requireNamespace("zinbwave", quietly = TRUE)) {
        stop("The ZINB-WaVE simulation requires the 'zinbwave' package.")
    }

    params <- new("ZINBParams")

    params <- setParams(params, ...)

    return(params)
}

setValidity("ZINBParams", function(object) {

    v <- getParams(object, slotNames(object))

    checks <- c(nGenes = checkmate::checkInt(v$nGenes, lower = 1),
                nCells = checkmate::checkInt(v$nCells, lower = 1),
                seed = checkmate::checkInt(v$seed, lower = 0),
                model = checkmate::checkClass(v$model, "ZinbModel"),
                model_valid = validObject(v$model, test = TRUE))

    if (all(checks == TRUE)) {
        valid <- TRUE
    } else {
        valid <- checks[checks != TRUE]
        valid <- paste(names(valid), valid, sep = ": ")
    }

    return(valid)
})

#' @rdname setParam
setMethod("setParam", "ZINBParams", function(object, name, value) {
    checkmate::assertString(name)

    if (name %in% names(getSlots("ZinbModel"))) {
        model <- getParam(object, "model")
        slot(model, name) <- value
        object <- setParam(object, "model", model)
    } else {
        object <- callNextMethod()
    }

    return(object)
})

setMethod("show", "ZINBParams", function(object) {

    pp <- list("Design:"         = c("(Samples)"       = "X",
                                     "(Genes)"         = "V"),
               "Offsets:"        = c("(Mu)"            = "O_mu",
                                     "(Pi)"            = "O_pi"),
               "Indices:"        = c("(Sample Mu)"     = "which_X_mu",
                                     "(Gene Mu)"       = "which_V_mu",
                                     "(Sample Pi)"     = "which_X_pi",
                                     "(Gene Pi)"       = "which_V_pi"),
               "Intercepts:"     = c("(Sample Mu)"     = "X_mu_intercept",
                                     "(Gene Mu)"       = "V_mu_intercept",
                                     "(Sample Pi)"     = "X_pi_intercept",
                                     "(Gene Pi)"       = "V_pi_intercept"),
               "Latent factors:" = c("(W)"             = "W"),
               "Coefficients:"   = c("(Sample Mu)"     = "beta_mu",
                                     "(Gene Mu)"       = "gamma_mu",
                                     "(Latent Mu)"     = "alpha_mu",
                                     "(Sample Pi)"     = "beta_pi",
                                     "(Gene Pi)"       = "gamma_pi",
                                     "(Latent Pi)"     = "alpha_pi"),
               "Regularisation:" = c("(Sample Mu)"     = "epsilon_beta_mu",
                                     "(Gene Mu)"       = "epsilon_gamma_mu",
                                     "(Sample Pi)"     = "epsilon_beta_pi",
                                     "(Gene Pi)"       = "epsilon_gamma_pi",
                                     "(Latent)"        = "epsilon_W",
                                     "(Latent coeffs)" = "epsilon_alpha",
                                     "(Zeta)"          = "epsilon_zeta",
                                     "(Logit)"         = "epsilon_min_logit"))

    callNextMethod()

    model <- getParam(object, "model")
    cat("Model:", "\n")
    cat("ZinbModel with", zinbwave::nFeatures(model), "features,",
        zinbwave::nSamples(model), "samples,", zinbwave::nFactors(model),
        "factors and", zinbwave::nParams(model), "parameters", "\n\n")

    default <- zinbwave::zinbModel()
    for (category in names(pp)) {
        parameters <- pp[[category]]
        values <- lapply(parameters, function(x) {slot(model, x)})
        short.values <- sapply(values, function(x) {
            if ("matrix" %in% class(x)) {
                if (nrow(x) == 1) {
                    paste0(paste(head(x[1, ], n = 4), collapse = ", "), ",...")
                } else if (ncol(x) == 1) {
                    paste0(paste(head(x[, 1], n = 4), collapse = ", "), ",...")
                } else {
                    paste(nrow(x), "x", ncol(x), "matrix")
                }
            } else if (length(x) > 4) {
                paste0(paste(head(x, n = 4), collapse = ", "), ",...")
            } else {
                paste(x, collapse = ", ")
            }
        })
        values <- sapply(values, paste, collapse = ", ")
        default.values <- lapply(parameters, function(x) {slot(default, x)})
        default.values <- sapply(default.values, paste, collapse = ", ")
        not.default <- values != default.values
        names(short.values)[not.default] <- toupper(names(values[not.default]))
        cat("Model", category, "\n")
        print(noquote(short.values), print.gap = 2)
        cat("\n")
    }

})
