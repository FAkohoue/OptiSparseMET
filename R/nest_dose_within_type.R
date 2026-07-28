#' Nest management doses within their type (fertiliser, pesticide, irrigation)
#'
#' @description
#' A raw dose column must never be used as a standalone environmental covariate:
#' 100 kg of NPK is not the same input as 100 kg of urea, and 10 mm of drip
#' irrigation is not 10 mm of manual spraying. Describing two sites by an
#' ungrouped dose would falsely treat different inputs as equal.
#' `nest_dose_within_type()` enforces the correct encoding -- it replaces each
#' dose with **type-nested dose** columns (the dose carried only within its own
#' fertiliser/pesticide/irrigation type, zero otherwise) and drops the raw dose,
#' so a dose is always interpreted relative to the product/method that delivered
#' it.
#'
#' @details
#' For a dose `dose` nested in a type `type` with levels \eqn{\ell}, the output
#' has one column per level, `type_<level>:dose`, equal to the dose where
#' `type == level` and 0 elsewhere (including where `type` is missing -- treated
#' as "no application"). This applies identically to fertiliser type, pesticide
#' type, and irrigation method (drip, sprinkler, manual spray, supplemental,
#' ...). The raw dose column is removed by default (`drop_raw_dose = TRUE`); the
#' type indicator is kept (`keep_type = TRUE`) so the presence of a type is still
#' captured for [build_enviromic_covariates()] to one-hot encode.
#'
#' @param management A data frame of per-site management (with an `environment`
#'   column).
#' @param nesting A named list/vector mapping each dose column (names) to its
#'   type column (values), e.g.
#'   `list(fert_dose = "fert_type", water_dose = "irrigation_method")`.
#' @param drop_raw_dose Remove the ungrouped dose column(s). Default `TRUE`
#'   (recommended -- a raw dose is not comparable across types).
#' @param keep_type Keep the type column(s) for one-hot encoding downstream.
#'   Default `TRUE`.
#' @return The `management` data frame with type-nested dose columns added (and
#'   the raw dose removed unless `drop_raw_dose = FALSE`).
#' @seealso [add_management_interactions()], [build_enviromic_covariates()].
#' @examples
#' mg <- data.frame(environment = c("E1", "E2", "E3"),
#'                  fert_type = c("NPK", "urea", "NPK"),
#'                  fert_dose = c(120, 80, 100),
#'                  irrigation_method = c("drip", "spray", "drip"),
#'                  water_dose = c(10, 20, 10))
#' nest_dose_within_type(mg, list(fert_dose = "fert_type",
#'                                water_dose = "irrigation_method"))
#' @export
nest_dose_within_type <- function(management, nesting,
                                  drop_raw_dose = TRUE, keep_type = TRUE) {
  df <- as.data.frame(management)
  if (is.null(names(nesting)) || any(!nzchar(names(nesting))))
    stop("`nesting` must be a named list/vector mapping dose columns to type columns.")

  for (dose in names(nesting)) {
    type <- nesting[[dose]]
    if (!dose %in% names(df)) stop("Dose column not found: ", dose, ".")
    if (!type %in% names(df)) stop("Type column not found: ", type, ".")
    d <- as.numeric(df[[dose]])
    f <- factor(df[[type]])
    for (lv in levels(f)) {
      ind <- !is.na(f) & f == lv                      # matching, non-missing type
      df[[paste0(type, "_", lv, ":", dose)]] <- ifelse(ind, d, 0)
    }
    if (drop_raw_dose) df[[dose]] <- NULL             # never keep ungrouped dose
    if (!keep_type)    df[[type]] <- NULL
  }
  df
}
