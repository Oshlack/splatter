#' @rdname newParams
#' @importFrom methods new
#' @export
newZINBParams <- function(...) {

    if (!requireNamespace("zinbwave", quietly = TRUE)) {
        stop("The ZINB-WaVE simulation requires the 'zinbwave' package.")
    }

    params <- new("ZINBParams")
    default.model <- zinbwave::zinbModel()

    params <- setParams(params, model = default.model)
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

    if (name %in% c("nGenes", "nCells")) {
        stop(name, " cannot be set directly, set model instead")
    }

    if (name == "model") {
        checkmate::assertClass(value, "ZinbModel")
        object <- setParamUnchecked(object, "nGenes",
                                    as.numeric(zinbwave::nFeatures(value)))
        object <- setParamUnchecked(object, "nCells",
                                    as.numeric(zinbwave::nSamples(value)))
    }

    object <- callNextMethod()

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
    default <- zinbwave::zinbModel()
    not.default <- !identical(model, default)
    cat(crayon::bold("Model:"), "\n")
    msg <- paste("ZinbModel with", zinbwave::nFeatures(model), "features,",
                 zinbwave::nSamples(model), "samples,",
                 zinbwave::nFactors(model), "latent factors and",
                 zinbwave::nParams(model), "parameters")
    if (not.default) {
        msg <- crayon::bold(crayon::green(msg))
    }
    cat(msg, "\n\n")


    cat(crayon::bold("Parameters of the ZinbModel"), "\n\n")
    for (category in names(pp)) {
        parameters <- pp[[category]]
        values <- lapply(parameters, function(x) {slot(model, x)})
        short.values <- vapply(values, function(x) {
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
        }, c(Value = "None"))
        values <- vapply(values, paste, c(Value = "None"), collapse = ", ")
        default.values <- lapply(parameters, function(x) {slot(default, x)})
        default.values <- vapply(default.values, paste, c(Value = "None"),
                                 collapse = ", ")
        not.default <- values != default.values
        cat(crayon::bold(c("Model", category)), "\n")
        showValues(short.values, not.default)
        cat("\n")
    }

})
