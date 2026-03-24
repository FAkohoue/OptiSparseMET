# Internal utilities for OptiSparseMET
# nocov start
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # field book columns (within-environment design)
    "Plot",
    "Row",
    "Col",
    "Block",
    "Replicate",
    "Entry",
    "Treatment",
    "Family",
    "Check",
    "IsCheck",
    "StreamPos",
    # MET-level allocation columns
    "Environment",
    "Allocation",
    "CommonTreatment",
    "TargetRep",
    "AssignedRep",
    "CoOccurrence",
    # seed-aware replication columns
    "SeedAvailable",
    "SeedRequired",
    "ReplicationRole",
    "SeedFeasible",
    # field book assembly
    "FieldBook",
    "EnvID",
    # relationship matrix helpers
    "..keep_cols",
    # dispersion optimisation
    "SpatialGroup",
    "RelatednessScore"
  ))
}
# nocov end

#' Internal debug message helper
#'
#' @param debug Logical.
#' @param ... Arguments passed to `sprintf()`.
#'
#' @keywords internal
.optisparsemet_dbg <- function(debug, ...) {
  if (isTRUE(debug)) message(sprintf(...))
}

#' Internal z-score helper
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector.
#' @keywords internal
.optisparsemet_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' Validate a square matrix
#'
#' @param M Matrix.
#' @param p Expected dimension.
#' @param nm Object name for messages.
#'
#' @keywords internal
.validate_square_matrix <- function(M, p, nm = "matrix") {
  if (!is.matrix(M)) stop(sprintf("%s must be a matrix.", nm), call. = FALSE)
  if (nrow(M) != p || ncol(M) != p) {
    stop(sprintf("%s must be a %d x %d matrix.", nm, p, p), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate that a named vector covers all required names
#'
#' @param x Named vector or list.
#' @param required Character vector of required names.
#' @param nm Object name for messages.
#'
#' @keywords internal
.validate_named_coverage <- function(x, required, nm = "object") {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(sprintf(
      "%s is missing required entries: %s",
      nm, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate that a data frame or data.table contains required columns
#'
#' @param dt A data.frame or data.table.
#' @param cols Character vector of required column names.
#' @param nm Object name for messages.
#'
#' @keywords internal
.validate_cols <- function(dt, cols, nm = "data") {
  missing <- setdiff(cols, names(dt))
  if (length(missing) > 0L) {
    stop(sprintf(
      "%s is missing required columns: %s",
      nm, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' OPTISPARSEMET: Sparse Multi-Environment Trial Design
#'
#' Package providing tools for:
#' \itemize{
#'   \item sparse treatment allocation across environments (M3 and M4 strategies),
#'   \item feasibility checking and capacity helpers,
#'   \item seed-aware replication planning,
#'   \item within-environment field layout construction,
#'   \item MET-level field book assembly,
#'   \item end-to-end pipeline execution.
#' }
#'
#' @keywords internal
"_PACKAGE"