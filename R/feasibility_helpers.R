#' Feasibility helpers for sparse MET allocation
#'
#' Utility functions for assessing whether a proposed sparse MET configuration
#' has sufficient per-environment capacity to assign every non-common treatment
#' at least once, and for deriving safe values for
#' `n_test_entries_per_environment` prior to calling [allocate_sparse_met()].
#'
#' These helpers address a specific failure mode: passing `n_test_entries_per_environment`
#' that is too small to cover all treatments, which causes [allocate_sparse_met()]
#' to either error, produce incomplete allocation, or require `allow_approximate = TRUE`
#' in ways that leave some treatments unassigned. Calling these helpers before
#' allocation lets the user verify feasibility and choose a safe capacity without
#' trial and error.
#'
#' The three exported helpers serve distinct purposes. [min_k_for_full_coverage()]
#' computes the arithmetic minimum. [suggest_safe_k()] adds a user-defined
#' buffer on top of that minimum and returns a single ready-to-use integer.
#' [warn_if_k_too_small()] checks an existing configuration and emits a warning
#' with a corrective suggestion if the capacity is insufficient. The internal
#' function [.check_full_coverage_feasibility()] is called by [allocate_sparse_met()]
#' to enforce the same condition with a hard stop.
#'
#' @name feasibility_helpers
NULL


#' Compute the minimum per-environment capacity for full treatment coverage
#'
#' Returns the minimum value of `n_test_entries_per_environment` under a
#' uniform capacity assumption that guarantees every non-common treatment can
#' be assigned to at least one environment. This is the floor below which any
#' sparse allocation must leave some treatments unassigned.
#'
#' @description
#' Let \eqn{J} be the total number of treatments, \eqn{C} the number of
#' common treatments assigned to every environment, and \eqn{I} the number of
#' environments. The number of sparse (non-common) treatments is:
#'
#' \deqn{J^* = J - C}
#'
#' Common treatments consume \eqn{C} slots in every environment before sparse
#' allocation begins. The number of sparse-allocatable slots per environment is:
#'
#' \deqn{k_e^* = k_e - C}
#'
#' For every sparse treatment to appear in at least one environment, the total
#' number of sparse slots across all environments must satisfy:
#'
#' \deqn{\sum_{e=1}^{I} k_e^* \geq J^*}
#'
#' Under equal per-environment capacity, the minimum number of sparse slots
#' per environment is:
#'
#' \deqn{k^*_{\min} = \left\lceil J^* / I \right\rceil}
#'
#' and the corresponding minimum total entries per environment is:
#'
#' \deqn{k_{\min} = k^*_{\min} + C}
#'
#' This function returns both quantities. The value `min_total_entries_per_environment`
#' is what should be passed to `n_test_entries_per_environment` to guarantee
#' full coverage at minimum capacity.
#'
#' @param n_treatments_total Positive integer. Total number of treatments
#'   including common treatments. Corresponds to `length(treatments)` in
#'   [allocate_sparse_met()].
#'
#' @param n_environments Positive integer. Number of environments in the MET.
#'   Must be at least 1.
#'
#' @param n_common_treatments Non-negative integer, default `0`. Number of
#'   treatments assigned to every environment before sparse allocation.
#'   Must not exceed `n_treatments_total`.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{`n_sparse_treatments`}{Integer. Number of non-common treatments
#'     \eqn{J^* = J - C}. These are the treatments subject to sparse
#'     allocation.}
#'   \item{`min_sparse_slots_per_environment`}{Integer. Minimum number of
#'     sparse-allocatable slots required per environment under equal capacities:
#'     \eqn{\lceil J^* / I \rceil}. This is the minimum value of
#'     \eqn{k_e^* = k_e - C}.}
#'   \item{`min_total_entries_per_environment`}{Integer. Minimum total entries
#'     per environment including common treatments:
#'     \eqn{\lceil J^* / I \rceil + C}. Pass this as
#'     `n_test_entries_per_environment` to guarantee that every non-common
#'     treatment can be assigned at least once.}
#' }
#'
#' @examples
#' ## No common treatments: minimum k is ceil(80 / 3) = 27
#' min_k_for_full_coverage(
#'   n_treatments_total  = 80,
#'   n_environments      = 3,
#'   n_common_treatments = 0
#' )
#'
#' ## With 5 common treatments: J* = 75, minimum sparse k = ceil(75 / 3) = 25,
#' ## minimum total k = 25 + 5 = 30
#' min_k_for_full_coverage(
#'   n_treatments_total  = 80,
#'   n_environments      = 3,
#'   n_common_treatments = 5
#' )
#'
#' @export
#' 
min_k_for_full_coverage <- function(
    n_treatments_total,
    n_environments,
    n_common_treatments = 0
) {
  if (!is.numeric(n_treatments_total) || length(n_treatments_total) != 1 ||
      is.na(n_treatments_total) || n_treatments_total < 1) {
    stop("`n_treatments_total` must be a single positive integer.")
  }
  
  if (!is.numeric(n_environments) || length(n_environments) != 1 ||
      is.na(n_environments) || n_environments < 1) {
    stop("`n_environments` must be a single positive integer.")
  }
  
  if (!is.numeric(n_common_treatments) || length(n_common_treatments) != 1 ||
      is.na(n_common_treatments) || n_common_treatments < 0) {
    stop("`n_common_treatments` must be a single non-negative integer.")
  }
  
  n_treatments_total  <- as.integer(n_treatments_total)
  n_environments      <- as.integer(n_environments)
  n_common_treatments <- as.integer(n_common_treatments)
  
  if (n_common_treatments > n_treatments_total) {
    stop("`n_common_treatments` cannot exceed `n_treatments_total`.")
  }
  
  n_sparse <- n_treatments_total - n_common_treatments
  k_sparse_min <- if (n_sparse > 0L) ceiling(n_sparse / n_environments) else 0L
  k_total_min <- k_sparse_min + n_common_treatments
  
  list(
    n_sparse_treatments = n_sparse,
    min_sparse_slots_per_environment = as.integer(k_sparse_min),
    min_total_entries_per_environment = as.integer(k_total_min)
  )
}


#' Suggest a safe uniform per-environment capacity for sparse MET allocation
#'
#' Returns a single integer suitable for passing to `n_test_entries_per_environment`
#' in [allocate_sparse_met()]. The value is the minimum capacity that guarantees
#' full treatment coverage plus a user-defined buffer, making it appropriate for
#' examples, tests, and demonstration workflows where a comfortably feasible
#' setting is preferred over the exact minimum.
#'
#' @description
#' The suggested value is:
#'
#' \deqn{k_{\text{safe}} = \left\lceil J^* / I \right\rceil + C + b}
#'
#' where \eqn{J^* = J - C} is the number of non-common treatments, \eqn{I} is
#' the number of environments, \eqn{C} is the number of common treatments, and
#' \eqn{b} is `buffer`. Setting `buffer = 0` returns the strict minimum from
#' [min_k_for_full_coverage()]. A buffer of 3 to 5 is typically sufficient to
#' absorb rounding and ensure [allocate_sparse_met()] runs without needing
#' `allow_approximate = TRUE`.
#'
#' @param treatments Character vector of treatment IDs. Duplicates are silently
#'   removed.
#'
#' @param environments Character vector of environment names. Duplicates are
#'   silently removed.
#'
#' @param common_treatments Optional character vector of common treatment IDs.
#'   Values not present in `treatments` are silently dropped. If `NULL`, no
#'   common treatments are assumed.
#'
#' @param buffer Non-negative integer, default `3`. Extra slots added on top of
#'   the strict minimum capacity. A small buffer reduces the risk of
#'   infeasibility from rounding and makes the configuration more robust to
#'   minor changes in trial dimensions.
#'
#' @return Integer scalar. A safe uniform value for
#'   `n_test_entries_per_environment`.
#'
#' @examples
#' trt <- paste0("L", sprintf("%03d", 1:80))
#' env <- c("Env1", "Env2", "Env3")
#'
#' ## No common treatments, default buffer of 3:
#' ## ceil(80 / 3) + 0 + 3 = 27 + 3 = 30
#' suggest_safe_k(
#'   treatments        = trt,
#'   environments      = env,
#'   common_treatments = NULL,
#'   buffer            = 3
#' )
#'
#' ## With 5 common treatments, buffer of 2:
#' ## ceil(75 / 3) + 5 + 2 = 25 + 5 + 2 = 32
#' suggest_safe_k(
#'   treatments        = trt,
#'   environments      = env,
#'   common_treatments = trt[1:5],
#'   buffer            = 2
#' )
#'
#' @export
#' 
suggest_safe_k <- function(
    treatments,
    environments,
    common_treatments = NULL,
    buffer = 3L
) {
  treatments   <- unique(as.character(treatments))
  environments <- unique(as.character(environments))
  
  if (length(treatments) < 1L) {
    stop("`treatments` must contain at least one treatment ID.")
  }
  if (length(environments) < 1L) {
    stop("`environments` must contain at least one environment.")
  }
  
  if (is.null(common_treatments)) {
    common_treatments <- character(0)
  } else {
    common_treatments <- intersect(unique(as.character(common_treatments)), treatments)
  }
  
  if (!is.numeric(buffer) || length(buffer) != 1 || is.na(buffer) || buffer < 0) {
    stop("`buffer` must be a single non-negative integer.")
  }
  buffer <- as.integer(buffer)
  
  min_info <- min_k_for_full_coverage(
    n_treatments_total = length(treatments),
    n_environments = length(environments),
    n_common_treatments = length(common_treatments)
  )
  
  as.integer(min_info$min_total_entries_per_environment + buffer)
}


#' Warn when per-environment capacity is insufficient for full treatment coverage
#'
#' Checks whether the supplied `n_test_entries_per_environment` provides enough
#' sparse slots to assign every non-common treatment at least once. If not,
#' emits a warning that states the deficit and the minimum capacity needed to
#' resolve it. Execution continues regardless — this function is a diagnostic
#' aid, not a hard stop.
#'
#' @description
#' The feasibility condition checked is:
#'
#' \deqn{\sum_{e=1}^{I} (k_e - C) \geq J^*}
#'
#' where \eqn{k_e} is the capacity of environment \eqn{e}, \eqn{C} is the
#' number of common treatments, and \eqn{J^*} is the number of non-common
#' treatments. When the condition fails, the warning message includes the
#' current total sparse slots, the number of sparse treatments, and the minimum
#' uniform capacity that would satisfy the condition.
#'
#' Use this function interactively before calling [allocate_sparse_met()] when
#' you want a non-fatal check, or inside wrapper functions that should guide the
#' user toward a valid configuration without stopping.
#'
#' @param treatments Character vector of treatment IDs. Duplicates are silently
#'   removed.
#'
#' @param environments Character vector of environment names. Duplicates are
#'   silently removed.
#'
#' @param n_test_entries_per_environment Integer scalar or integer vector.
#'   Total number of entries per environment including common treatments. If a
#'   scalar, the same capacity is applied to all environments. If a vector, its
#'   length must equal the number of environments. All values must be at least
#'   `length(common_treatments)`.
#'
#' @param common_treatments Optional character vector of common treatment IDs.
#'   Values not present in `treatments` are silently dropped. If `NULL`, no
#'   common treatments are assumed.
#'
#' @return Invisibly returns `NULL`. Called for its side effect (warning) when
#'   the capacity is insufficient.
#'
#' @examples
#' trt <- paste0("L", sprintf("%03d", 1:80))
#' env <- c("Env1", "Env2", "Env3")
#'
#' ## k = 20: total sparse slots = 3 x 20 = 60 < 80 — warns and suggests k = 27
#' warn_if_k_too_small(
#'   treatments                     = trt,
#'   environments                   = env,
#'   n_test_entries_per_environment = 20
#' )
#'
#' ## k = 27: total sparse slots = 3 x 27 = 81 >= 80 — no warning
#' warn_if_k_too_small(
#'   treatments                     = trt,
#'   environments                   = env,
#'   n_test_entries_per_environment = 27
#' )
#'
#' ## Heterogeneous capacities: total = 20 + 25 + 20 = 65 < 80 — warns
#' warn_if_k_too_small(
#'   treatments                     = trt,
#'   environments                   = env,
#'   n_test_entries_per_environment = c(20, 25, 20)
#' )
#'
#' @export
#' 
warn_if_k_too_small <- function(
    treatments,
    environments,
    n_test_entries_per_environment,
    common_treatments = NULL
) {
  treatments   <- unique(as.character(treatments))
  environments <- unique(as.character(environments))
  
  if (length(treatments) < 1L) {
    stop("`treatments` must contain at least one treatment ID.")
  }
  if (length(environments) < 1L) {
    stop("`environments` must contain at least one environment.")
  }
  
  if (is.null(common_treatments)) {
    common_treatments <- character(0)
  } else {
    common_treatments <- intersect(unique(as.character(common_treatments)), treatments)
  }
  
  n_env <- length(environments)
  n_common <- length(common_treatments)
  n_sparse <- length(treatments) - n_common
  
  if (length(n_test_entries_per_environment) == 1L) {
    k_vec <- rep(as.integer(n_test_entries_per_environment), n_env)
  } else {
    k_vec <- as.integer(n_test_entries_per_environment)
    if (length(k_vec) != n_env) {
      stop("`n_test_entries_per_environment` must have length 1 or length(environments).")
    }
  }
  
  if (any(is.na(k_vec)) || any(k_vec < n_common)) {
    stop("All environments must have capacity at least equal to the number of common treatments.")
  }
  
  total_sparse_slots <- sum(k_vec - n_common)
  
  if (total_sparse_slots < n_sparse) {
    min_info <- min_k_for_full_coverage(
      n_treatments_total = length(treatments),
      n_environments = n_env,
      n_common_treatments = n_common
    )
    
    warning(
      paste0(
        "Current design is infeasible for full treatment coverage: ",
        total_sparse_slots, " sparse slots available for ", n_sparse,
        " non-common treatments. ",
        "For a uniform design, use at least ",
        min_info$min_total_entries_per_environment,
        " entries per environment."
      ),
      call. = FALSE
    )
  }
  
  invisible(NULL)
}


#' Enforce full-coverage feasibility for sparse MET allocation
#'
#' Internal function called by [allocate_sparse_met()] to verify that the
#' supplied `n_test_entries_per_environment` provides enough sparse slots to
#' assign every non-common treatment at least once. Stops with an informative
#' error when the condition is not met, naming the deficit and the minimum
#' uniform capacity that would resolve it.
#'
#' Unlike [warn_if_k_too_small()], this function stops rather than warns and is
#' not intended for direct use. It is documented here for package maintainers.
#'
#' @param treatments Character vector of all treatment IDs.
#' @param environments Character vector of environment names.
#' @param n_test_entries_per_environment Integer scalar or integer vector.
#' @param common_treatments Character vector of common treatment IDs. Values
#'   not in `treatments` are silently dropped.
#'
#' @return Invisibly returns `TRUE` when the condition is satisfied. Stops with
#'   an error otherwise.
#'
#' @keywords internal
#' 
.check_full_coverage_feasibility <- function(
    treatments,
    environments,
    n_test_entries_per_environment,
    common_treatments = character(0)
) {
  treatments   <- unique(as.character(treatments))
  environments <- unique(as.character(environments))
  common_treatments <- intersect(unique(as.character(common_treatments)), treatments)
  
  n_env <- length(environments)
  n_common <- length(common_treatments)
  n_sparse <- length(treatments) - n_common
  
  if (length(n_test_entries_per_environment) == 1L) {
    k_vec <- rep(as.integer(n_test_entries_per_environment), n_env)
  } else {
    k_vec <- as.integer(n_test_entries_per_environment)
    if (length(k_vec) != n_env) {
      stop("`n_test_entries_per_environment` must have length 1 or length(environments).")
    }
  }
  
  if (any(is.na(k_vec)) || any(k_vec < n_common)) {
    stop("Each environment must have capacity at least equal to the number of common treatments.")
  }
  
  total_sparse_slots <- sum(k_vec - n_common)
  
  if (total_sparse_slots < n_sparse) {
    min_info <- min_k_for_full_coverage(
      n_treatments_total = length(treatments),
      n_environments = n_env,
      n_common_treatments = n_common
    )
    
    stop(
      paste0(
        "Not enough slots to place every non-common treatment at least once. ",
        "Available sparse slots = ", total_sparse_slots,
        ", sparse treatments = ", n_sparse, ". ",
        "For a uniform design, the minimum `n_test_entries_per_environment` is ",
        min_info$min_total_entries_per_environment, "."
      )
    )
  }
  
  invisible(TRUE)
}