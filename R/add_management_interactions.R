#' Build management interaction covariates
#'
#' @description
#' Some management effects only make sense in combination -- fertiliser *dose*
#' matters differently under each fertiliser *type*, or planting *mode* interacts
#' with *irrigation*. `add_management_interactions()` adds those interaction
#' columns to a management data frame before it is passed to
#' [build_enviromic_covariates()], so the environmental covariate matrix can
#' capture them.
#'
#' @section When to use (and when it is redundant):
#' Explicit interactions are **not always needed**. If the environmental kinship
#' uses the Gaussian (RBF) kernel in [build_environment_relationship()], that
#' kernel already *represents* covariate interactions implicitly, so adding them
#' here is usually redundant and can even degrade the kernel (more dimensions,
#' distance concentration). Reach for this function when you (i) use the linear
#' `"correlation"` kernel, which captures no implicit interactions; (ii) want to
#' emphasise or select a specific agronomic combination as a limiting factor
#' (via the `variables`/`weights` arguments of
#' [build_environment_relationship()]); or (iii) need interpretable interaction
#' features for reporting or a reaction-norm model. Keep the list short and
#' agronomy- or [screen_interactions()]-driven. Management **doses** are a
#' separate matter -- always encode them with [nest_dose_within_type()], which
#' corrects a mis-encoded variable (100 kg NPK is not 100 kg urea) rather than
#' adding a modelling interaction.
#'
#' @details
#' Each interaction is the element-wise product of its variables after expansion:
#' a numeric variable is used as-is; a categorical variable is expanded to its
#' one-hot indicators. So `list(c("fert_type", "fert_dose"))` produces one column
#' per fertiliser type holding the dose within that type (0 otherwise), e.g.
#' `fert_type_NPK:fert_dose` and `fert_type_urea:fert_dose`. Two categorical
#' variables give the combined-level indicators; two numeric variables give their
#' product. Higher-order interactions (three or more variables) are supported by
#' listing them together. The original main-effect columns are kept.
#'
#' @param management A data frame with an `environment` column plus management
#'   covariates.
#' @param interactions A list, each element a character vector of two or more
#'   column names in `management` to interact.
#' @param drop_na Replace non-finite products (from missing inputs) with 0.
#'   Default `TRUE`.
#' @return The `management` data frame with the interaction columns appended
#'   (names joined by `":"`).
#' @seealso [build_enviromic_covariates()].
#' @examples
#' mg <- data.frame(environment = c("E1", "E2", "E3"),
#'                  fert_type = c("NPK", "urea", "NPK"),
#'                  fert_dose = c(120, 80, 100),
#'                  planting_mode = c("transplanting", "direct", "direct"))
#' add_management_interactions(mg, list(c("fert_type", "fert_dose")))
#' @export
add_management_interactions <- function(management, interactions,
                                        drop_na = TRUE) {
  df <- as.data.frame(management)
  if (!is.list(interactions)) interactions <- list(interactions)

  expand_var <- function(v) {
    if (!v %in% names(df)) stop("Interaction variable not found: ", v, ".")
    col <- df[[v]]
    if (is.numeric(col)) {
      stats::setNames(list(as.numeric(col)), v)
    } else {
      f <- factor(col)
      stats::setNames(lapply(levels(f), function(lv) as.integer(f == lv)),
                      paste0(v, "_", levels(f)))
    }
  }

  for (term in interactions) {
    if (length(term) < 2L)
      stop("Each interaction needs at least two variables.")
    parts <- lapply(term, expand_var)
    combos <- Reduce(function(a, b) {
      out <- list()
      for (na in names(a)) for (nb in names(b))
        out[[paste0(na, ":", nb)]] <- a[[na]] * b[[nb]]
      out
    }, parts)
    for (nm in names(combos)) {
      v <- combos[[nm]]
      if (drop_na) v[!is.finite(v)] <- 0
      df[[nm]] <- v
    }
  }
  df
}
