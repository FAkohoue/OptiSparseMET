#' Evaluate feasibility of an exact balanced incomplete sparse MET allocation
#'
#' `check_balanced_incomplete_feasibility()` determines whether the parameters
#' of a proposed sparse MET design admit an exact balanced incomplete
#' allocation — that is, whether the total number of sparse-allocatable
#' treatment slots across environments exactly equals the number of non-common
#' treatments multiplied by the target replication. The function is a
#' diagnostic companion to `allocate_sparse_met()` with
#' `allocation_method = "balanced_incomplete"`, and should be called before
#' allocation when the user wants to verify feasibility or understand the
#' magnitude of any imbalance before deciding whether to set
#' `allow_approximate = TRUE`.
#'
#' @description
#' For a balanced incomplete allocation to be exact, the slot identity:
#'
#' \deqn{J^* \times r = \sum_{e=1}^{I} k_e^*}
#'
#' must hold, where \eqn{J^*} is the number of non-common treatments,
#' \eqn{r} is the target number of environments per treatment, \eqn{I} is
#' the number of environments, and \eqn{k_e^*} is the number of
#' sparse-allocatable slots in environment \eqn{e} after common treatments
#' have been subtracted from its total capacity. When environments have
#' equal capacity \eqn{k^*}, this reduces to the standard BIBD identity
#' \eqn{J^* \times r = I \times k^*}.
#'
#' The function evaluates this condition and returns both a logical indicator
#' and the signed difference \eqn{\sum k_e^* - J^* \times r}, which quantifies
#' the degree of imbalance when exact feasibility fails. A positive difference
#' means there are more available slots than required — some treatments will
#' receive an extra replication in an approximate allocation. A negative
#' difference means slots are insufficient — some treatments will be
#' under-replicated. Both cases can be handled by `allocate_sparse_met()` with
#' `allow_approximate = TRUE`, but the magnitude of the difference informs
#' whether the resulting imbalance is practically acceptable.
#'
#' @details
#' Common treatments are assigned to every environment before sparse
#' allocation begins. Their count is subtracted from each environment's total
#' capacity to obtain the sparse-allocatable slots:
#'
#' \deqn{k_e^* = k_e - C}
#'
#' where \eqn{k_e} is the total number of test entries assigned to environment
#' \eqn{e} (`n_test_entries_per_environment`) and \eqn{C} is the number of
#' common treatments (`n_common_treatments`). The number of non-common
#' treatments is:
#'
#' \deqn{J^* = J - C}
#'
#' where \eqn{J} is `n_treatments_total`. If any environment has
#' \eqn{k_e < C}, that environment cannot accommodate all common treatments
#' and the function stops with an error before evaluating the balance
#' condition.
#'
#' This function performs no allocation. It only evaluates the arithmetic
#' condition and returns the components needed to diagnose feasibility. For
#' the actual construction of the incidence matrix, see
#' `allocate_sparse_met()`.
#'
#' @param n_treatments_total Positive integer. Total number of test treatments,
#'   including any common treatments. This is the full candidate pool before
#'   any subdivision into common and sparse subsets.
#'
#' @param n_environments Integer greater than or equal to 2. Number of
#'   environments in the trial.
#'
#' @param n_test_entries_per_environment Integer scalar or integer vector.
#'   Total number of test treatments assigned to each environment, including
#'   common treatments. If a scalar, the same capacity is applied to all
#'   environments. If a vector, its length must equal `n_environments`. All
#'   values must be positive and must not be less than `n_common_treatments`.
#'
#' @param target_replications Positive integer. Desired number of environments
#'   in which each non-common treatment should appear. This is the \eqn{r}
#'   term in the slot identity. For exact feasibility, the product
#'   \eqn{J^* \times r} must equal the total number of sparse slots
#'   \eqn{\sum k_e^*}.
#'
#' @param n_common_treatments Non-negative integer, default `0`. Number of
#'   treatments forced into all environments before sparse allocation. These
#'   treatments consume capacity in every environment and are excluded from
#'   the sparse allocation pool. Setting this to `0` corresponds to a design
#'   with no forced common entries.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`feasible`}{Logical. `TRUE` if and only if the slot identity
#'     holds exactly, i.e. `difference == 0`.}
#'   \item{`n_sparse_treatments`}{Integer. Number of non-common treatments
#'     \eqn{J^* = J - C}. These are the treatments subject to sparse
#'     allocation.}
#'   \item{`k_sparse`}{Integer vector of length `n_environments`. Sparse-
#'     allocatable slots per environment \eqn{k_e^* = k_e - C}.}
#'   \item{`total_sparse_slots`}{Integer. Total available slots for sparse
#'     allocation: \eqn{\sum_{e=1}^{I} k_e^*}.}
#'   \item{`required_sparse_slots`}{Integer. Slots required for exact balance:
#'     \eqn{J^* \times r}.}
#'   \item{`difference`}{Integer. Signed difference
#'     \eqn{\sum k_e^* - J^* \times r}. Zero when allocation is exactly
#'     feasible. Positive when slots exceed requirement; negative when slots
#'     are insufficient.}
#'   \item{`message`}{Character scalar. Human-readable summary of the
#'     feasibility check, including the values of the key quantities and
#'     the direction of any imbalance.}
#' }
#'
#' @references
#' Montesinos-López, O. A., Mosqueda-González, B. A., Salinas-Ruiz, J.,
#' Montesinos-López, A., & Crossa, J. (2023). Sparse multi-trait genomic
#' prediction under balanced incomplete block design. *The Plant Genome*,
#' 16, e20305.
#'
#' @examples
#' ## Example 1: exact feasibility — 110 sparse treatments x 2 reps = 220,
#' ## 4 environments x 55 sparse slots each = 220.
#' check_balanced_incomplete_feasibility(
#'   n_treatments_total             = 120,
#'   n_environments                 = 4,
#'   n_test_entries_per_environment = 65,
#'   target_replications            = 2,
#'   n_common_treatments            = 10
#' )
#' # feasible = TRUE, difference = 0
#'
#' ## Example 2: slot deficit — 110 sparse treatments x 2 reps = 220,
#' ## 4 environments x 50 sparse slots each = 200. The deficit of 20 means
#' ## approximately 20 treatments will receive only 1 replication under an
#' ## approximate allocation.
#' check_balanced_incomplete_feasibility(
#'   n_treatments_total             = 120,
#'   n_environments                 = 4,
#'   n_test_entries_per_environment = 60,
#'   target_replications            = 2,
#'   n_common_treatments            = 10
#' )
#' # feasible = FALSE, difference = -20
#'
#' ## Example 3: slot surplus — 110 sparse treatments x 2 reps = 220,
#' ## 4 environments x 60 sparse slots each = 240. The surplus of 20 means
#' ## approximately 20 treatments will receive a third replication under an
#' ## approximate allocation.
#' check_balanced_incomplete_feasibility(
#'   n_treatments_total             = 120,
#'   n_environments                 = 4,
#'   n_test_entries_per_environment = 70,
#'   target_replications            = 2,
#'   n_common_treatments            = 10
#' )
#' # feasible = FALSE, difference = 20
#'
#' ## Example 4: heterogeneous environment capacities — useful when
#' ## environments differ in the number of plots available.
#' check_balanced_incomplete_feasibility(
#'   n_treatments_total             = 100,
#'   n_environments                 = 4,
#'   n_test_entries_per_environment = c(40, 45, 40, 45),
#'   target_replications            = 3,
#'   n_common_treatments            = 5
#' )
#'
#' @export
#' 
check_balanced_incomplete_feasibility <- function(
    n_treatments_total,
    n_environments,
    n_test_entries_per_environment,
    target_replications,
    n_common_treatments = 0L
) {
  
  if (!is.numeric(n_treatments_total) || length(n_treatments_total) != 1 ||
      is.na(n_treatments_total) || n_treatments_total < 1) {
    stop("`n_treatments_total` must be a single positive integer.")
  }
  if (!is.numeric(n_environments) || length(n_environments) != 1 ||
      is.na(n_environments) || n_environments < 2) {
    stop("`n_environments` must be a single integer >= 2.")
  }
  if (!is.numeric(target_replications) || length(target_replications) != 1 ||
      is.na(target_replications) || target_replications < 1) {
    stop("`target_replications` must be a single positive integer.")
  }
  if (!is.numeric(n_common_treatments) || length(n_common_treatments) != 1 ||
      is.na(n_common_treatments) || n_common_treatments < 0) {
    stop("`n_common_treatments` must be a single integer >= 0.")
  }
  
  n_treatments_total <- as.integer(n_treatments_total)
  n_environments <- as.integer(n_environments)
  target_replications <- as.integer(target_replications)
  n_common_treatments <- as.integer(n_common_treatments)
  
  if (length(n_test_entries_per_environment) == 1) {
    k_vec <- rep(as.integer(n_test_entries_per_environment), n_environments)
  } else {
    k_vec <- as.integer(n_test_entries_per_environment)
    if (length(k_vec) != n_environments) {
      stop("If `n_test_entries_per_environment` is a vector, its length must match `n_environments`.")
    }
  }
  
  if (any(is.na(k_vec)) || any(k_vec < 1)) {
    stop("All values in `n_test_entries_per_environment` must be positive integers.")
  }
  
  if (n_common_treatments > n_treatments_total) {
    stop("`n_common_treatments` cannot exceed `n_treatments_total`.")
  }
  if (any(k_vec < n_common_treatments)) {
    stop("At least one environment has fewer slots than the number of common treatments.")
  }
  
  n_sparse_treatments <- n_treatments_total - n_common_treatments
  k_sparse <- k_vec - n_common_treatments
  total_sparse_slots <- sum(k_sparse)
  required_sparse_slots <- n_sparse_treatments * target_replications
  difference <- total_sparse_slots - required_sparse_slots
  feasible <- (difference == 0L)
  
  msg <- if (feasible) {
    paste0(
      "Exact balanced incomplete allocation is feasible: ",
      n_sparse_treatments, " sparse treatments x ",
      target_replications, " replications = ",
      required_sparse_slots, " sparse slots, matching the available total."
    )
  } else {
    paste0(
      "Exact balanced incomplete allocation is not feasible: available sparse slots = ",
      total_sparse_slots, ", required sparse slots = ",
      required_sparse_slots, ", difference = ", difference, "."
    )
  }
  
  list(
    feasible = feasible,
    n_sparse_treatments = n_sparse_treatments,
    k_sparse = k_sparse,
    total_sparse_slots = total_sparse_slots,
    required_sparse_slots = required_sparse_slots,
    difference = difference,
    message = msg
  )
}