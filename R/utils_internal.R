# ==============================================================================
# utils_internal.R
# Internal utilities for OptiSparseMET
#
# Contents:
#   - Global variable declarations (suppresses R CMD CHECK NOTEs)
#   - .optisparsemet_dbg()       -- conditional debug message helper
#   - .optisparsemet_z()         -- z-score normalisation helper
#   - .validate_square_matrix()  -- square matrix validator
#   - .validate_named_coverage() -- named vector coverage validator
#   - .validate_cols()           -- data frame column validator
# ==============================================================================

# nocov start
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(

    # -- Field book columns (within-environment design) -----------------------
    # Columns present in field books returned by met_prep_famoptg() and
    # met_alpha_rc_stream(). Declared to suppress R CMD CHECK NOTEs when these
    # names appear unquoted in NSE contexts.
    "Plot",
    "Row",
    "Col",
    "Column",
    "Block",
    "IBlock",
    "BlockInRep",
    "Rep",
    "Replicate",
    "Entry",
    "Treatment",
    "Family",
    "Gcluster",
    "Check",
    "IsCheck",
    "StreamPos",

    # -- MET-level allocation columns ----------------------------------------
    # Columns added by allocate_sparse_met() and combine_met_fieldbooks().
    "Environment",
    "Allocation",
    "CommonTreatment",
    "IsCommonTreatment",
    "TargetRep",
    "AssignedRep",
    "CoOccurrence",
    "AllocationGroup",
    "LocalDesign",
    "ReplicationMode",
    "SparseMethod",

    # -- Seed-aware replication columns --------------------------------------
    # Columns produced by assign_replication_by_seed().
    "SeedAvailable",
    "SeedRequired",
    "ReplicationRole",
    "SeedFeasible",

    # -- Field book assembly -------------------------------------------------
    "FieldBook",
    "EnvID",

    # -- Relationship matrix helpers -----------------------------------------
    "..keep_cols",

    # -- Dispersion optimisation ---------------------------------------------
    "SpatialGroup",
    "RelatednessScore",

    # -- pracma import -------------------------------------------------------
    # mod() is used for serpentine traversal parity in met_prep_famoptg()
    # and met_alpha_rc_stream().
    "mod"
  ))
}
# nocov end


# ------------------------------------------------------------------------------
# .optisparsemet_dbg
# ------------------------------------------------------------------------------
#' Internal conditional debug message helper
#'
#' Emits a formatted message via [base::message()] when `debug = TRUE`.
#' Used internally to provide optional verbose output during design
#' construction and optimisation.
#'
#' @param debug Logical. If `TRUE`, the message is emitted.
#' @param ... Arguments passed to [base::sprintf()].
#'
#' @return Invisibly returns `NULL`.
#' @keywords internal
.optisparsemet_dbg <- function(debug, ...) {
  if (isTRUE(debug)) message(sprintf(...))
}


# ------------------------------------------------------------------------------
# .optisparsemet_z
# ------------------------------------------------------------------------------
#' Internal z-score normalisation helper
#'
#' Computes the z-score of a numeric vector. If the standard deviation is
#' zero or non-finite, returns a zero vector of the same length.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector of z-scores, or zeros if `sd(x) == 0`.
#' @keywords internal
.optisparsemet_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}


# ------------------------------------------------------------------------------
# .validate_square_matrix
# ------------------------------------------------------------------------------
#' Internal square matrix validator
#'
#' Checks that an object is a base R matrix with the expected square dimensions.
#' Used internally to validate relationship matrices (GRM, A, K).
#'
#' @param M Object to validate.
#' @param p Expected dimension. Must have `nrow(M) == p` and `ncol(M) == p`.
#' @param nm Character scalar. Object name used in error messages.
#'
#' @return Invisibly returns `TRUE` if validation passes.
#' @keywords internal
.validate_square_matrix <- function(M, p, nm = "matrix") {
  if (!is.matrix(M))
    stop(sprintf("%s must be a matrix.", nm), call. = FALSE)
  if (nrow(M) != p || ncol(M) != p)
    stop(sprintf("%s must be a %d x %d matrix.", nm, p, p), call. = FALSE)
  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# .validate_named_coverage
# ------------------------------------------------------------------------------
#' Internal named vector coverage validator
#'
#' Checks that a named vector or list contains all required names.
#' Used internally to validate that `env_design_specs` covers all environments.
#'
#' @param x Named vector or list.
#' @param required Character vector of required names.
#' @param nm Character scalar. Object name used in error messages.
#'
#' @return Invisibly returns `TRUE` if validation passes.
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


# ------------------------------------------------------------------------------
# .validate_cols
# ------------------------------------------------------------------------------
#' Internal data frame column validator
#'
#' Checks that a data frame contains all required column names.
#' Used internally across allocation and seed-planning functions.
#'
#' @param dt A data frame.
#' @param cols Character vector of required column names.
#' @param nm Character scalar. Object name used in error messages.
#'
#' @return Invisibly returns `TRUE` if validation passes.
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


# ------------------------------------------------------------------------------
# .normalise_tpe_weights
# ------------------------------------------------------------------------------
# Keep the definition of the target population of environments identical in
# the information, simulation, and optimisation paths.
.normalise_tpe_weights <- function(tpe_weights, envs) {
  E <- length(envs)
  if (is.null(tpe_weights)) {
    return(stats::setNames(rep(1 / E, E), envs))
  }
  if (!is.numeric(tpe_weights) || length(tpe_weights) != E ||
      any(!is.finite(tpe_weights)) || any(tpe_weights < 0) ||
      sum(tpe_weights) <= 0)
    stop("`tpe_weights` must contain one finite non-negative value per ",
         "environment and must not sum to zero.")
  if (!is.null(names(tpe_weights))) {
    if (any(names(tpe_weights) == "") || anyDuplicated(names(tpe_weights)) ||
        !setequal(names(tpe_weights), envs))
      stop("Named `tpe_weights` must uniquely cover every environment.")
    tpe_weights <- tpe_weights[envs]
  }
  stats::setNames(as.numeric(tpe_weights) / sum(tpe_weights), envs)
}


# ------------------------------------------------------------------------------
# .normalise_environment_variance
# ------------------------------------------------------------------------------
# Accept a scalar or an environment-specific residual-variance vector.
.normalise_environment_variance <- function(x, envs, arg = "sigma_e2") {
  E <- length(envs)
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x <= 0))
    stop("`", arg, "` must contain finite positive values.")
  if (!is.null(names(x))) {
    if (any(names(x) == "") || anyDuplicated(names(x)) ||
        !all(envs %in% names(x)))
      stop("Named `", arg, "` must uniquely cover every environment.")
    x <- x[envs]
  } else if (length(x) == 1L) {
    x <- rep(x, E)
  } else if (length(x) != E) {
    stop("`", arg, "` must be scalar or contain one value per environment.")
  }
  stats::setNames(as.numeric(x), envs)
}
