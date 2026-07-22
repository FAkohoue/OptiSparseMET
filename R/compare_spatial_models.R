#' Score one design under several spatial models
#'
#' @description
#' A design that is best under one spatial assumption is not guaranteed to be
#' best under another, and at design time the spatial structure is rarely known.
#' `compare_spatial_models()` evaluates a single field layout under a set of
#' spatial models and returns the criteria side by side, so a design can be
#' checked for robustness rather than tuned to one assumption. It is the spatial
#' analogue of [sensitivity_varcomp()].
#'
#' @details
#' Available models:
#' \describe{
#'   \item{`"IID"`}{Independent residuals -- no spatial structure.}
#'   \item{`"AR1"`, `"AR1xAR1"`}{Separable autoregressive residual correlation
#'     along rows, or along rows and columns.}
#'   \item{`"AR1xAR1_nugget"`}{Separable AR1xAR1 plus an independent
#'     measurement-error (nugget) component -- the standard applied model, since
#'     pure AR1xAR1 often fits field data poorly.}
#'   \item{`"exponential"`, `"gaussian"`, `"matern"`}{Isotropic correlation as a
#'     function of Euclidean plot distance, useful for irregular layouts.}
#'   \item{`"pspline"`}{A SpATS-style two-dimensional penalised-spline surface
#'     fitted as a random effect with IID residuals (Rodriguez-Alvarez et al.
#'     2018), rather than as a residual covariance.}
#' }
#'
#' Because the criteria are conditional on the assumed parameters, compare
#' designs on the same settings. A design that wins across every model is a safe
#' choice; one that wins under a single model is a warning.
#'
#' @param field_book Field book to evaluate, as passed to the evaluator.
#' @param n_rows,n_cols Field dimensions.
#' @param check_treatments Character vector of check identifiers.
#' @param models Character vector of models to evaluate (see Details).
#' @param evaluator Which engine to use: `"alpha"` for
#'   [met_evaluate_alpha_efficiency()] or `"famoptg"` for
#'   [met_evaluate_famoptg_efficiency()].
#' @param rho_row,rho_col AR1 parameters used by the AR1-type models.
#' @param nugget Nugget proportion used by `"AR1xAR1_nugget"` and the kernels.
#' @param kernel_range,matern_nu Parameters for the distance kernels.
#' @param spline_lambda_row,spline_lambda_col Smoothing parameters for
#'   `"pspline"`.
#' @param ... Further arguments passed to the evaluator (for example
#'   `varcomp`, `treatment_effect`, `prediction_type`, `K`).
#'
#' @return A data frame with one row per model and columns `model`,
#'   `A_criterion`, `D_criterion`, `CDmean` and `mean_PEV` (`NA` where a
#'   criterion does not apply), plus `error` for any model that failed.
#'
#' @references
#' Rodriguez-Alvarez, M. X., Boer, M. P., van Eeuwijk, F. A., & Eilers, P. H. C.
#' (2018). Correcting for spatial heterogeneity in plant breeding experiments
#' with P-splines. *Spatial Statistics*, 23, 52-71.
#'
#' @seealso [sensitivity_varcomp()], [met_evaluate_alpha_efficiency()].
#' @export
compare_spatial_models <- function(field_book, n_rows, n_cols,
                                   check_treatments,
                                   models = c("IID", "AR1", "AR1xAR1",
                                              "AR1xAR1_nugget", "exponential",
                                              "pspline"),
                                   evaluator = c("alpha", "famoptg"),
                                   rho_row = 0.3, rho_col = 0.3,
                                   nugget = 0.2,
                                   kernel_range = 3, matern_nu = 1.5,
                                   spline_lambda_row = 1,
                                   spline_lambda_col = 1,
                                   ...) {
  evaluator <- match.arg(evaluator)
  fun <- if (evaluator == "alpha") met_evaluate_alpha_efficiency else
    met_evaluate_famoptg_efficiency

  pick <- function(x, nm) if (!is.null(x[[nm]])) x[[nm]] else NA_real_

  rows <- lapply(models, function(m) {
    args <- list(field_book = field_book, n_rows = n_rows, n_cols = n_cols,
                 check_treatments = check_treatments, ...)
    if (m == "pspline") {
      args$residual_structure   <- "IID"
      args$spatial_random       <- "pspline"
      args$spline_lambda_row    <- spline_lambda_row
      args$spline_lambda_col    <- spline_lambda_col
    } else {
      args$residual_structure <- m
      args$rho_row <- rho_row; args$rho_col <- rho_col
      if (m == "AR1xAR1_nugget") args$nugget <- nugget
      if (m %in% c("exponential", "gaussian", "matern")) {
        args$kernel_range <- kernel_range
        args$matern_nu    <- matern_nu
        args$nugget       <- nugget
      }
    }
    eff <- tryCatch(do.call(fun, args), error = function(e) e)
    if (inherits(eff, "error")) {
      data.frame(model = m, A_criterion = NA_real_, D_criterion = NA_real_,
                 CDmean = NA_real_, mean_PEV = NA_real_,
                 error = conditionMessage(eff), stringsAsFactors = FALSE)
    } else {
      data.frame(model = m,
                 A_criterion = pick(eff, "A_criterion"),
                 D_criterion = pick(eff, "D_criterion"),
                 CDmean      = pick(eff, "CDmean"),
                 mean_PEV    = pick(eff, "mean_PEV"),
                 error       = NA_character_, stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}
