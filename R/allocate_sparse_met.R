#' Allocate test treatments across environments for sparse MET designs
#'
#' `allocate_sparse_met()` constructs the treatment-by-environment incidence
#' structure for a sparse multi-environment trial. Given a set of candidate
#' lines and a set of environments, the function determines which lines enter
#' which environments subject to replication targets, capacity constraints,
#' optional common treatment requirements, and optional genetic structure
#' constraints. The function guarantees that every non-common treatment is
#' assigned to at least one environment before attempting to fill additional
#' replication slots.
#'
#' @description
#' Two allocation strategies are available. `"random_balanced"` implements an
#' M3-inspired stochastic allocation that approximates balance without requiring
#' exact BIBD parameters to be satisfiable -- appropriate when environment
#' capacities differ or when the trial dimensions do not admit an exact balanced
#' solution. Unlike the original M3 of Montesinos-Lopez et al. (2023), which
#' allocates location by location and may silently leave some lines with fewer
#' than \eqn{r} replications, this implementation uses a two-phase
#' coverage-first strategy that additionally guarantees every non-common
#' treatment appears in at least one environment before replication filling
#' begins.
#'
#' The \code{"balanced_incomplete"} method implements the M4 allocation of
#' Montesinos-Lopez et al. (2023): every non-common treatment appears in
#' exactly \eqn{r} environments (equal replication) and every environment
#' receives exactly \eqn{k^*} sparse treatments (equal environment sizes),
#' so the resource identity \eqn{J^* \times r = I \times k^*} holds exactly.
#' This equal-replication, equal-environment-size guarantee is what
#' distinguishes M4 from M3 in the paper, and is always the goal in plant
#' breeding programs where thousands of lines are tested across a few
#' environments.
#' 
#' #' `allow_approximate = FALSE` (the default) is the standard M4 path: the
#' slot identity must hold exactly, and the function stops with an informative
#' error if it cannot be met, so the caller always knows whether equal
#' replication was achieved. Construction uses a greedy load-balanced
#' constructor that assigns each sparse treatment to environments in
#' decreasing order of remaining capacity, guaranteeing equal replication
#' for every non-common treatment. `allow_approximate = TRUE` relaxes the
#' slot identity and allows minor replication imbalances across lines; it is
#' a fallback for exploratory use, not the intended primary path.
#' 
#' Before allocation begins, the function calls
#' `.check_full_coverage_feasibility()` to verify that the total number of
#' sparse slots across all environments is sufficient to assign every non-common
#' treatment at least once. If it is not, the function stops with an informative
#' error naming the deficit and the minimum capacity that would resolve it. Use
#' [suggest_safe_k()] or [min_k_for_full_coverage()] to choose a feasible
#' `n_test_entries_per_environment` before calling this function.
#'
#' @details
#' ## Relation to sparse testing theory
#'
#' The allocation strategies implemented here correspond to M3 and M4 in
#' Montesinos-Lopez et al. (2023). The underlying resource identity is:
#'
#' \deqn{N = J \times r = I \times k}
#'
#' where \eqn{J} is the number of lines, \eqn{I} is the number of environments,
#' \eqn{k} is the number of lines per environment, \eqn{r} is the number of
#' environments each line enters, and \eqn{N} is the total number of
#' line-by-environment observations. This identity makes the tradeoff between
#' coverage (\eqn{k}) and replication depth (\eqn{r}) explicit given fixed
#' \eqn{N}.
#'
#' ## Common treatments
#'
#' Treatments in `common_treatments` are assigned to every environment before
#' sparse allocation begins and are excluded from the sparse allocation pool.
#' The per-environment sparse capacity is:
#'
#' \deqn{k_e^* = k_e - C}
#'
#' where \eqn{C} is the number of common treatments and \eqn{k_e} is the total
#' capacity of environment \eqn{e}. Common treatments establish direct
#' cross-environment connectivity that does not depend on the relationship
#' structure in the covariance model, and provide stable references for
#' estimating environment-level effects. They are most important when
#' environments are weakly correlated.
#'
#' ## Two-phase allocation
#'
#' For `"random_balanced"` and approximate `"balanced_incomplete"`, phase one
#' iterates over sparse treatments in random order and assigns each to one
#' environment, preferring environments where the treatment's genetic group is
#' not yet represented when `allocation_group_source` is active. Phase two
#' iterates over environments in decreasing order of remaining capacity and
#' fills each to its target, guided by line-level replication deficit scores and
#' group-balance penalties.
#'
#' ## Within-environment layout
#'
#' Once allocation is complete, the set of treatments assigned to each
#' environment is passed to the within-environment design function.
#' [met_prep_famoptg()] constructs repeated-check augmented, partially
#' replicated (p-rep), and RCBD-type block designs. [met_alpha_rc_stream()]
#' constructs stream-based alpha row-column designs for fixed-grid field
#' geometry. Both accept output from [assign_replication_by_seed()] directly.
#' The end-to-end pipeline that coordinates allocation and local design
#' construction is [plan_sparse_met_design()].
#'
#' ## Allocation groups
#'
#' When `allocation_group_source` is not `"none"`, both phases use genetic
#' group membership to guide assignments. Family-based grouping reads labels
#' from `treatment_info$Family`. GRM-based and A-based grouping derive clusters
#' from the leading principal components of the respective relationship matrix.
#' In phase two, the allocator penalizes assignments that concentrate a group
#' in too few environments (`balance_groups_across_env`) and gives preference to
#' groups that have not yet reached `min_env_per_group`
#' (`force_group_connectivity`). Line-level replication targets are preserved
#' within these group-level constraints.
#'
#' ## M4 allocation and the role of allow_approximate
#'
#' The M4 method (Montesinos-Lopez et al., 2023) enforces two conditions:
#'
#' \enumerate{
#'   \item Equal replication: every non-common treatment appears in exactly
#'         \eqn{r} environments.
#'   \item Equal environment sizes: every environment receives exactly
#'         \eqn{k^*} sparse treatments, so \eqn{J^* \times r = I \times k^*}.
#' }
#' 
#' `allow_approximate = FALSE` (the default) enforces both conditions strictly.
#' If the slot identity \eqn{J^* \times r = I \times k^*} does not hold for
#' the chosen `n_test_entries_per_environment` and `target_replications`, the
#' function stops with an informative error. Use
#' [check_balanced_incomplete_feasibility()] to verify the slot identity before
#' calling, or adjust \eqn{k} and \eqn{r} so that \eqn{J^* \times r =
#' I \times k^*}. Construction uses a greedy load-balanced constructor that
#' assigns each sparse treatment to the least-loaded eligible environments,
#' guaranteeing equal replication for every non-common treatment.
#'
#' `allow_approximate = TRUE` relaxes the slot identity: the function
#' constructs the most balanced allocation it can without stopping on
#' infeasibility, accepting that some lines may receive more or fewer
#' replications than \eqn{r}. This mode is useful for exploratory analysis
#' but does not provide the equal-replication guarantee that is the defining
#' property of M4.
#'
#' @param treatments Character vector of test treatment IDs to allocate across
#'   environments. Check treatments should not be included here; they are
#'   handled separately within the within-environment design functions. Duplicate
#'   values are silently removed.
#'
#' @param environments Character vector of environment names. Must contain at
#'   least two elements. Duplicate values are silently removed.
#'
#' @param allocation_method Character scalar. Sparse allocation strategy.
#'   Accepted values are `"random_balanced"` (M3-inspired stochastic
#'   allocation) and `"balanced_incomplete"` (M4-type BIBD allocation). The
#'   aliases `"M3"` and `"M4"` are also accepted and translated internally to
#'   their canonical names before any further processing.
#'
#' @param n_test_entries_per_environment Integer scalar or integer vector. Total
#'   number of test treatments per environment including common treatments. If a
#'   scalar, it applies uniformly to all environments. If a vector, its length
#'   must equal the number of environments. All values must be positive and must
#'   be at least `length(common_treatments)`. Use [suggest_safe_k()] to choose
#'   a value that guarantees full treatment coverage.
#'
#' @param target_replications Optional positive integer. Target number of
#'   environments in which each non-common treatment should appear. For
#'   `"random_balanced"`, this is a soft target that guides the phase-two
#'   filling. For `"balanced_incomplete"`, this is the strict replication level
#'   required for an exact balanced solution. If `NULL`, the function infers a
#'   value as `floor(total_sparse_slots / n_sparse_treatments)`, with a minimum
#'   of 1.
#'
#' @param common_treatments Optional character vector of treatment IDs to assign
#'   to all environments. Placed before sparse allocation begins and excluded
#'   from the sparse pool. Values not present in `treatments` are silently
#'   dropped.
#'
#' @param allocation_group_source Character scalar. Controls whether and how
#'   genetic group structure guides allocation. `"none"` disables group-guided
#'   allocation. `"Family"` uses `treatment_info$Family`. `"GRM"` derives
#'   clusters from `GRM`. `"A"` derives clusters from `A`. Active in both
#'   allocation phases when set to anything other than `"none"`.
#'
#' @param treatment_info Optional data frame. Required when
#'   `allocation_group_source = "Family"`. Must contain columns `Treatment` and
#'   `Family`. When `allocation_group_source %in% c("GRM", "A")`, this argument
#'   is optional but used to anchor the number of clusters when supplied.
#'
#' @param GRM Optional numeric matrix. Genomic relationship matrix. Required
#'   when `allocation_group_source = "GRM"`. Row and column names must match
#'   treatment IDs or be reachable through `id_map`.
#'
#' @param A Optional numeric matrix. Pedigree-based numerator relationship
#'   matrix. Required when `allocation_group_source = "A"`. Same naming
#'   requirements as `GRM`.
#'
#' @param id_map Optional data frame with columns `Treatment` and `LineID`.
#'   Required only when `allocation_group_source %in% c("GRM", "A")` and
#'   treatment IDs do not match the row/column names of the relationship matrix.
#'
#' @param group_method Character scalar. Clustering algorithm applied to the
#'   principal components of `GRM` or `A`. `"kmeans"` or `"hclust"`. Ignored
#'   when `allocation_group_source %in% c("none", "Family")`.
#'
#' @param group_seed Integer. Seed for k-means initialization. Active only when
#'   `allocation_group_source %in% c("GRM", "A")` and
#'   `group_method = "kmeans"`.
#'
#' @param group_attempts Integer. Number of random restarts for k-means.
#'   Active only when `allocation_group_source %in% c("GRM", "A")` and
#'   `group_method = "kmeans"`.
#'
#' @param n_pcs_use Integer or `Inf`. Number of leading principal components
#'   retained for matrix-based clustering. `Inf` retains all components
#'   corresponding to positive eigenvalues. Ignored when
#'   `allocation_group_source %in% c("none", "Family")`.
#'
#' @param min_groups_per_environment Optional positive integer. Minimum number
#'   of allocation groups that should be represented in each environment, where
#'   feasible given available treatments. Active in phase two when
#'   `allocation_group_source` is not `"none"`.
#'
#' @param min_env_per_group Optional positive integer. Minimum number of
#'   environments in which each allocation group should appear, where feasible.
#'   Used in phase two when `force_group_connectivity = TRUE`.
#'
#' @param balance_groups_across_env Logical. If `TRUE`, phase two preferentially
#'   assigns treatments from groups whose current allocation count is below
#'   their proportional target. Active when `allocation_group_source` is not
#'   `"none"`.
#'
#' @param force_group_connectivity Logical. If `TRUE`, phase two preferentially
#'   assigns treatments from groups that have not yet appeared in
#'   `min_env_per_group` environments. Active when `allocation_group_source` is
#'   not `"none"` and `min_env_per_group` is not `NULL`.
#'
#' @param allow_approximate Logical, default `FALSE`. When `FALSE` and
#'   `allocation_method = "balanced_incomplete"`, the slot identity
#'   \eqn{J^* \times r = I \times k^*} must hold exactly; the function stops
#'   with an error if it does not. This is the standard M4 path and guarantees
#'   equal replication for every non-common treatment. When `TRUE`, the slot
#'   identity is not enforced and minor replication imbalances are accepted; this
#'   is a relaxed fallback for exploratory use, not the primary mode. For
#'   `"random_balanced"`, this argument has no effect.
#'
#' @param seed Optional integer. Random seed for reproducibility. Controls
#'   the random order in which sparse treatments are processed in phase one
#'   and the stochastic tie-breaking in phase two under `"random_balanced"`,
#'   as well as tie-breaking in the strict exact constructor. If `NULL`, no
#'   seed is set and results may differ across runs; the seed used internally
#'   is returned as `seed_used`.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`allocation_matrix`}{Binary integer matrix of dimension
#'     `n_treatments x n_environments` with `dimnames` set to `treatments` and
#'     `environments`. Entry `[i, e]` is `1L` if treatment `i` is assigned to
#'     environment `e`, and `0L` otherwise. Every non-common treatment is
#'     guaranteed to have row sum at least 1.}
#'   \item{`allocation_long`}{Long-format data frame with one row per
#'     treatment-by-environment combination. Columns: `Treatment`,
#'     `Environment`, `Assigned` (integer 0/1), `IsCommonTreatment` (logical),
#'     and `AllocationGroup` (character, present when
#'     `allocation_group_source` is not `"none"`).}
#'   \item{`overlap_matrix`}{Square integer matrix of dimension
#'     `n_environments x n_environments` giving the number of treatments shared
#'     between each pair of environments. Diagonal entries give the total
#'     treatment count per environment.}
#'   \item{`line_replications`}{Named integer vector of length `n_treatments`.
#'     Row sums of `allocation_matrix`: the number of environments each
#'     treatment enters.}
#'   \item{`environment_sizes`}{Named integer vector of length `n_environments`.
#'     Column sums of `allocation_matrix`: the number of treatments assigned to
#'     each environment.}
#'   \item{`group_assignment`}{Data frame with columns `Treatment` and
#'     `AllocationGroup`, one row per treatment. `NULL` when
#'     `allocation_group_source = "none"`.}
#'   \item{`group_by_environment`}{Data frame summarizing the count of assigned
#'     treatments from each allocation group in each environment. Columns:
#'     `Environment`, `AllocationGroup`, `n_treatments`. `NULL` when
#'     `allocation_group_source = "none"`.}
#'   \item{`group_overlap_matrix`}{Square integer matrix of dimension
#'     `n_environments x n_environments` giving the number of shared allocation
#'     groups between each pair of environments. `NULL` when
#'     `allocation_group_source = "none"`.}
#'   \item{`summary`}{Named list with allocation metadata: `allocation_method`,
#'     `allocation_group_source`, `target_replications`, `n_treatments_total`,
#'     `n_sparse_treatments`, `n_common_treatments`, `total_sparse_slots`,
#'     `environment_sizes`, `min_replication`, `max_replication`,
#'     `mean_replication`, `min_sparse_replication`, `max_sparse_replication`,
#'     `mean_sparse_replication`, `min_common_replication`,
#'     `max_common_replication`, and `mean_common_replication`.}
#'   \item{`seed_used`}{The integer seed passed to `set.seed()` internally, or
#'     `NULL` if no seed was supplied.}
#' }
#'
#' @seealso
#' [suggest_safe_k()] and [min_k_for_full_coverage()] for choosing a feasible
#' `n_test_entries_per_environment` before calling this function.
#' [check_balanced_incomplete_feasibility()] for verifying the slot identity
#' before attempting a `"balanced_incomplete"` allocation.
#' [derive_allocation_groups()] for inspecting the group structure that guides
#' allocation when `allocation_group_source` is not `"none"`.
#' [met_prep_famoptg()] and [met_alpha_rc_stream()] for the within-environment
#' design functions that consume the allocation output.
#' [plan_sparse_met_design()] for the end-to-end two-stage MET pipeline.
#'
#' @references
#' Montesinos-Lopez, O. A., Mosqueda-Gonzalez, B. A., Salinas-Ruiz, J.,
#' Montesinos-Lopez, A., & Crossa, J. (2023). Sparse multi-trait genomic
#' prediction under balanced incomplete block design. \emph{The Plant Genome},
#' 16, e20305. \doi{10.1002/tpg2.20305}
#'
#' @examples
#' treatments <- paste0("L", sprintf("%03d", 1:120))
#' envs       <- c("Env1", "Env2", "Env3", "Env4")
#' fam        <- rep(paste0("F", 1:6), each = 20)
#'
#' treatment_info <- data.frame(
#'   Treatment = treatments,
#'   Family    = fam,
#'   stringsAsFactors = FALSE
#' )
#'
#' ## Check a safe per-environment capacity before running allocation
#' suggest_safe_k(treatments, envs, buffer = 3)  # 33
#'
#' ## Example 1: random balanced allocation with family-guided group spreading
#' out1 <- allocate_sparse_met(
#'   treatments                     = treatments,
#'   environments                   = envs,
#'   allocation_method              = "random_balanced",
#'   n_test_entries_per_environment = 33,
#'   target_replications            = 1,
#'   allocation_group_source        = "Family",
#'   treatment_info                 = treatment_info,
#'   min_groups_per_environment     = 4,
#'   min_env_per_group              = 2,
#'   seed                           = 123
#' )
#'
#' out1$summary
#' head(out1$allocation_long)
#' head(out1$group_by_environment)
#'
#' ## Example 2: M4 balanced incomplete allocation (paper method).
#' ## Equal replication (r=2) and equal environment sizes (k*=55 per env).
#' ## Slot identity: J* x r = I x k* => 110 x 2 = 4 x 55 = 220. Valid.
#' ## allow_approximate = FALSE enforces the slot identity strictly, stopping
#' ## with an error if it is not met (the default safe behaviour).
#' out2 <- allocate_sparse_met(
#'   treatments                     = treatments,
#'   environments                   = envs,
#'   allocation_method              = "balanced_incomplete",
#'   n_test_entries_per_environment = 65,
#'   target_replications            = 2,
#'   common_treatments              = treatments[1:10],
#'   allow_approximate              = FALSE,
#'   seed                           = 123
#' )
#'
#' out2$summary
#' # Every sparse line appears in exactly 2 environments
#' range(out2$line_replications[!(names(out2$line_replications) %in%
#'                                 treatments[1:10])])
#'
#' @export
allocate_sparse_met <- function(
    treatments,
    environments,
    allocation_method = c("random_balanced", "balanced_incomplete", "M3", "M4"),
    n_test_entries_per_environment,
    target_replications = NULL,
    common_treatments = NULL,
    allocation_group_source = c("none", "Family", "GRM", "A"),
    treatment_info = NULL,
    GRM = NULL,
    A = NULL,
    id_map = NULL,
    group_method = c("kmeans", "hclust"),
    group_seed = 1,
    group_attempts = 25,
    n_pcs_use = Inf,
    min_groups_per_environment = NULL,
    min_env_per_group = NULL,
    balance_groups_across_env = TRUE,
    force_group_connectivity = TRUE,
    allow_approximate = FALSE,
    seed = NULL
) {
  
  # ============================================================
  # 0. RNG
  # ============================================================
  seed_used <- seed
  if (!is.null(seed_used)) set.seed(seed_used)
  
  allocation_method <- match.arg(allocation_method)
  if (allocation_method == "M3") allocation_method <- "random_balanced"
  if (allocation_method == "M4") allocation_method <- "balanced_incomplete"
  
  allocation_group_source <- match.arg(allocation_group_source)
  group_method            <- match.arg(group_method)
  
  # ============================================================
  # 1. Basic validation and normalisation
  # ============================================================
  treatments   <- unique(as.character(treatments))
  environments <- unique(as.character(environments))
  
  if (length(treatments) < 1L)   stop("`treatments` must contain at least one treatment ID.")
  if (length(environments) < 2L) stop("`environments` must contain at least two environment names.")
  
  n_treat <- length(treatments)
  n_env   <- length(environments)
  
  if (length(n_test_entries_per_environment) == 1L) {
    k_vec <- rep(as.integer(n_test_entries_per_environment), n_env)
  } else {
    k_vec <- as.integer(n_test_entries_per_environment)
  }
  
  if (length(k_vec) != n_env)
    stop("`n_test_entries_per_environment` must have length 1 or length(environments).")
  if (any(is.na(k_vec)) || any(k_vec < 1L))
    stop("All values of `n_test_entries_per_environment` must be positive integers.")
  
  if (!is.null(min_groups_per_environment)) {
    if (!is.numeric(min_groups_per_environment) || length(min_groups_per_environment) != 1L ||
        is.na(min_groups_per_environment) || min_groups_per_environment < 1L)
      stop("`min_groups_per_environment` must be NULL or a single positive integer.")
    min_groups_per_environment <- as.integer(min_groups_per_environment)
  }
  
  if (!is.null(min_env_per_group)) {
    if (!is.numeric(min_env_per_group) || length(min_env_per_group) != 1L ||
        is.na(min_env_per_group) || min_env_per_group < 1L)
      stop("`min_env_per_group` must be NULL or a single positive integer.")
    min_env_per_group <- as.integer(min_env_per_group)
  }
  
  if (!(is.numeric(n_pcs_use) && length(n_pcs_use) == 1L &&
        (is.finite(n_pcs_use) || is.infinite(n_pcs_use)) && n_pcs_use > 0))
    stop("`n_pcs_use` must be a single positive number or Inf.")
  
  # ============================================================
  # 2. Common treatments
  # ============================================================
  if (is.null(common_treatments)) {
    common_treatments <- character(0)
  } else {
    common_treatments <- intersect(treatments, unique(as.character(common_treatments)))
  }
  
  n_common          <- length(common_treatments)
  sparse_treatments <- setdiff(treatments, common_treatments)
  n_sparse          <- length(sparse_treatments)
  
  .check_full_coverage_feasibility(
    treatments                     = treatments,
    environments                   = environments,
    n_test_entries_per_environment = k_vec,
    common_treatments              = common_treatments
  )
  
  k_sparse           <- k_vec - n_common
  total_sparse_slots <- sum(k_sparse)
  
  # ============================================================
  # 3. Replication target
  # ============================================================
  if (is.null(target_replications)) {
    target_replications <- if (n_sparse == 0L) 0L else max(1L, floor(total_sparse_slots / n_sparse))
  } else {
    if (!is.numeric(target_replications) || length(target_replications) != 1L ||
        is.na(target_replications) || target_replications < 1L)
      stop("`target_replications` must be NULL or a single positive integer.")
    target_replications <- as.integer(target_replications)
  }
  
  # ============================================================
  # 4. M4 feasibility checks
  # ============================================================
  if (allocation_method == "balanced_incomplete") {
    
    required_slots <- n_sparse * target_replications
    
    if (!allow_approximate && required_slots != total_sparse_slots) {
      stop(paste0(
        "Exact balanced incomplete allocation infeasible. ",
        "Required sparse slots = ", required_slots,
        ", available sparse slots = ", total_sparse_slots, ". ",
        "Adjust n_test_entries_per_environment, target_replications, or set ",
        "allow_approximate = TRUE."
      ))
    }
  }
  
  # ============================================================
  # 5. Optional group derivation
  # ============================================================
  group_assignment     <- NULL
  sparse_groups        <- NULL
  unique_sparse_groups <- character(0)
  n_groups             <- 0L
  
  if (allocation_group_source != "none") {
    
    group_assignment_sparse <- derive_allocation_groups(
      treatments              = sparse_treatments,
      allocation_group_source = allocation_group_source,
      treatment_info          = treatment_info,
      GRM                     = GRM,
      A                       = A,
      id_map                  = id_map,
      group_method            = group_method,
      group_seed              = group_seed,
      group_attempts          = group_attempts,
      n_pcs_use               = n_pcs_use
    )
    
    if (!is.data.frame(group_assignment_sparse) ||
        !all(c("Treatment", "AllocationGroup") %in% names(group_assignment_sparse)))
      stop("`derive_allocation_groups()` must return a data frame with `Treatment` and `AllocationGroup`.")
    if (nrow(group_assignment_sparse) != n_sparse)
      stop("`derive_allocation_groups()` did not return one row per sparse treatment.")
    if (any(is.na(group_assignment_sparse$AllocationGroup)))
      stop("Missing allocation groups detected for sparse treatments.")
    
    group_assignment <- data.frame(
      Treatment       = treatments,
      AllocationGroup = NA_character_,
      stringsAsFactors = FALSE
    )
    group_assignment$AllocationGroup[
      match(group_assignment_sparse$Treatment, group_assignment$Treatment)
    ] <- as.character(group_assignment_sparse$AllocationGroup)
    
    if (n_common > 0L) {
      common_grp <- if (!is.null(treatment_info) && is.data.frame(treatment_info) &&
                        allocation_group_source == "Family" &&
                        all(c("Treatment", "Family") %in% names(treatment_info))) {
        grp <- treatment_info$Family[match(common_treatments, treatment_info$Treatment)]
        grp[is.na(grp)] <- "COMMON"
        grp
      } else {
        rep("COMMON", n_common)
      }
      group_assignment$AllocationGroup[
        match(common_treatments, group_assignment$Treatment)
      ] <- common_grp
    }
    
    sparse_groups        <- stats::setNames(
      group_assignment$AllocationGroup[match(sparse_treatments, group_assignment$Treatment)],
      sparse_treatments
    )
    unique_sparse_groups <- unique(unname(sparse_groups))
    n_groups             <- length(unique_sparse_groups)
  }
  
  # ============================================================
  # 6. Initialise allocation matrix
  # ============================================================
  alloc <- matrix(0L, nrow = n_treat, ncol = n_env,
                  dimnames = list(treatments, environments))
  if (n_common > 0L) alloc[common_treatments, ] <- 1L
  
  # ============================================================
  # 7. Allocation
  # ============================================================
  # ------------------------------------------------------------------
  # 7A. Strict exact M4 constructor (allow_approximate = FALSE)
  # ------------------------------------------------------------------
  if (allocation_method == "balanced_incomplete" && !allow_approximate) {
    
    sparse_env_load <- stats::setNames(integer(n_env), environments)
    target_env_load <- stats::setNames(as.integer(k_sparse), environments)
    line_rep        <- stats::setNames(integer(n_sparse), sparse_treatments)
    
    group_env_current <- if (n_groups > 0L) {
      matrix(
        0L,
        nrow = n_groups,
        ncol = n_env,
        dimnames = list(unique_sparse_groups, environments)
      )
    } else {
      NULL
    }
    
    sparse_order <- sample(sparse_treatments)
    
    for (trt in sparse_order) {
      chosen_envs <- character(0)
      
      for (rr in seq_len(target_replications)) {
        candidate_envs <- environments[
          !(environments %in% chosen_envs) &
            sparse_env_load < target_env_load
        ]
        
        if (length(candidate_envs) == 0L) {
          stop(
            paste0(
              "Cannot meet strict balanced incomplete replication target. ",
              "The slot totals are feasible, but the exact constructor ran out of ",
              "candidate environments. This indicates an implementation bug ",
              "for this slot-feasible case."
            )
          )
        }
        
        if (n_groups > 0L) {
          grp         <- sparse_groups[[trt]]
          grp_presence <- group_env_current[grp, candidate_envs]
          pref        <- ifelse(grp_presence == 0L, 1L, 0L)
          ord         <- order(-pref, sparse_env_load[candidate_envs], candidate_envs)
          chosen_env  <- candidate_envs[ord][1L]
          group_env_current[grp, chosen_env] <- 1L
        } else {
          min_load   <- min(sparse_env_load[candidate_envs])
          best_envs  <- candidate_envs[sparse_env_load[candidate_envs] == min_load]
          chosen_env <- sample(best_envs, 1L)
        }
        
        chosen_envs                     <- c(chosen_envs, chosen_env)
        alloc[trt, chosen_env]          <- 1L
        sparse_env_load[chosen_env]     <- sparse_env_load[chosen_env] + 1L
        line_rep[trt]                   <- line_rep[trt] + 1L
      }
    }
    
    if (!all(line_rep == target_replications)) {
      stop(
        "Strict balanced incomplete allocation failed: sparse treatments do not all have equal replication."
      )
    }
    
    if (!all(sparse_env_load == target_env_load)) {
      stop(
        "Strict balanced incomplete allocation failed: sparse environment loads are not exact."
      )
    }
    
    env_load <- colSums(alloc)
    
  } else {
    
    # ------------------------------------------------------------------
    # 7B. Coverage-first + greedy fill (M3 and approximate M4)
    # ------------------------------------------------------------------
    env_load <- colSums(alloc)
    
    if (n_sparse > 0L) {
      
      group_env_current <- if (n_groups > 0L)
        matrix(0L, nrow = n_groups, ncol = n_env,
               dimnames = list(unique_sparse_groups, environments))
      else NULL
      
      line_rep <- stats::setNames(integer(n_sparse), sparse_treatments)
      
      # Phase 1: force minimum coverage
      for (trt in sample(sparse_treatments)) {
        spare          <- k_vec - env_load
        candidate_envs <- environments[spare > 0L]
        
        if (length(candidate_envs) == 0L)
          stop("Internal allocation error: insufficient capacity during forced minimum coverage.")
        
        if (n_groups > 0L) {
          grp          <- sparse_groups[[trt]]
          grp_presence <- group_env_current[grp, candidate_envs]
          pref         <- ifelse(grp_presence == 0L, 1L, 0L)
          ord          <- order(-pref, env_load[candidate_envs], candidate_envs)
          chosen_env   <- candidate_envs[ord][1L]
          group_env_current[grp, chosen_env] <- 1L
        } else {
          chosen_env <- candidate_envs[which.min(env_load[candidate_envs])]
        }
        
        alloc[trt, chosen_env] <- 1L
        env_load[chosen_env]   <- env_load[chosen_env] + 1L
        line_rep[trt]          <- 1L
      }
      
    } else {
      line_rep          <- stats::setNames(integer(0), character(0))
      group_env_current <- NULL
    }
    
    remaining_slots <- k_vec - env_load
    env_order       <- order(remaining_slots, decreasing = TRUE)
    
    group_current <- if (n_groups > 0L) {
      gc             <- stats::setNames(integer(n_groups), unique_sparse_groups)
      sparse_assigned <- sparse_treatments[line_rep[sparse_treatments] > 0L]
      if (length(sparse_assigned) > 0L) {
        tt <- table(sparse_groups[sparse_assigned])
        gc[names(tt)] <- as.integer(tt)
      }
      gc
    } else integer(0)
    
    group_sizes  <- if (n_groups > 0L) as.integer(table(sparse_groups)) else integer(0)
    group_target <- if (n_groups > 0L)
      stats::setNames(group_sizes * target_replications, names(table(sparse_groups)))
    else integer(0)
    
    # Phase 2: fill remaining capacity
    for (e in env_order) {
      env_name <- environments[e]
      
      while (env_load[env_name] < k_vec[e]) {
        candidates <- sparse_treatments[alloc[sparse_treatments, env_name] == 0L]
        if (length(candidates) == 0L) break
        
        deficit <- target_replications - line_rep[candidates]
        score   <- as.numeric(deficit)
        names(score) <- candidates
        
        if (n_groups > 0L && balance_groups_across_env) {
          cand_groups <- sparse_groups[candidates]
          size_lookup <- stats::setNames(group_sizes, names(group_target))
          grp_term    <- pmax(0, as.numeric(
            group_target[cand_groups] - group_current[cand_groups]
          )) / pmax(1, as.numeric(size_lookup[cand_groups]))
          grp_term[is.na(grp_term)] <- 0
          score <- score + grp_term
        }
        
        if (n_groups > 0L && !is.null(min_groups_per_environment)) {
          env_groups_now   <- unique(sparse_groups[sparse_treatments[
            alloc[sparse_treatments, env_name] == 1L
          ]])
          n_env_groups_now <- sum(!is.na(env_groups_now) & nzchar(env_groups_now))
          cand_groups      <- sparse_groups[candidates]
          lacking_bonus    <- ifelse(!(cand_groups %in% env_groups_now), 1, 0)
          bonus_weight     <- if (n_env_groups_now < min(min_groups_per_environment, n_groups))
            100 else 5
          score <- score + bonus_weight * lacking_bonus
        }
        
        if (n_groups > 0L && force_group_connectivity && !is.null(min_env_per_group)) {
          env_count_group <- rowSums(group_env_current > 0L)
          cand_groups     <- sparse_groups[candidates]
          conn_term       <- 50 * pmax(0, as.numeric(
            min_env_per_group - env_count_group[cand_groups]
          ))
          conn_term[is.na(conn_term)] <- 0
          score <- score + conn_term
        }
        
        if (allocation_method == "balanced_incomplete" && allow_approximate)
          score[deficit[names(score)] < 0] <- score[deficit[names(score)] < 0] - 1000
        
        max_score <- max(score)
        best      <- names(score)[score >= (max_score - 1e-8)]
        
        chosen <- if (allocation_method == "random_balanced") {
          sample(best, 1L)
        } else {
          ord <- order(-score[best], line_rep[best], best)
          best[ord][1L]
        }
        
        alloc[chosen, env_name] <- 1L
        line_rep[chosen]        <- line_rep[chosen] + 1L
        env_load[env_name]      <- env_load[env_name] + 1L
        
        if (n_groups > 0L) {
          gp                             <- sparse_groups[[chosen]]
          group_current[gp]              <- group_current[gp] + 1L
          group_env_current[gp, env_name] <- 1L
        }
      }
    }
  }
  
  # ============================================================
  # 8. Outputs
  # ============================================================
  allocation_long <- expand.grid(Treatment   = treatments,
                                 Environment = environments,
                                 stringsAsFactors = FALSE)
  allocation_long$Assigned <- as.integer(
    alloc[cbind(allocation_long$Treatment, allocation_long$Environment)]
  )
  allocation_long$IsCommonTreatment <- allocation_long$Treatment %in% common_treatments
  
  if (!is.null(group_assignment)) {
    allocation_long$AllocationGroup <- group_assignment$AllocationGroup[
      match(allocation_long$Treatment, group_assignment$Treatment)
    ]
  }
  
  line_replications <- rowSums(alloc)
  environment_sizes <- colSums(alloc)
  overlap_matrix    <- t(alloc) %*% alloc
  
  group_by_environment <- NULL
  group_overlap_matrix <- NULL
  
  if (!is.null(group_assignment)) {
    gbe <- allocation_long[
      allocation_long$Assigned == 1L,
      c("Environment", "AllocationGroup", "Treatment"),
      drop = FALSE
    ]
    
    if (nrow(gbe) > 0L) {
      group_by_environment <- stats::aggregate(
        Treatment ~ Environment + AllocationGroup, data = gbe, FUN = length
      )
      names(group_by_environment)[
        names(group_by_environment) == "Treatment"
      ] <- "n_treatments"
      
      all_groups       <- unique(group_assignment$AllocationGroup[
        !is.na(group_assignment$AllocationGroup)
      ])
      grp_env_incidence <- matrix(
        0L, nrow = n_env, ncol = length(all_groups),
        dimnames = list(environments, all_groups)
      )
      for (i in seq_len(nrow(group_by_environment)))
        grp_env_incidence[
          group_by_environment$Environment[i],
          group_by_environment$AllocationGroup[i]
        ] <- 1L
      
      group_overlap_matrix <- grp_env_incidence %*% t(grp_env_incidence)
      
    } else {
      group_by_environment <- data.frame(
        Environment = character(0), AllocationGroup = character(0),
        n_treatments = integer(0), stringsAsFactors = FALSE
      )
      group_overlap_matrix <- matrix(
        0L, nrow = n_env, ncol = n_env,
        dimnames = list(environments, environments)
      )
    }
  }
  
  sparse_replications <- if (n_sparse > 0L) line_replications[sparse_treatments] else integer(0)
  common_replications <- if (n_common > 0L) line_replications[common_treatments] else integer(0)
  
  
  summary_out <- list(
    allocation_method        = allocation_method,
    allocation_group_source  = allocation_group_source,
    target_replications      = target_replications,
    n_treatments_total       = n_treat,
    n_sparse_treatments      = n_sparse,
    n_common_treatments      = n_common,
    total_sparse_slots       = total_sparse_slots,
    environment_sizes        = environment_sizes,
    min_replication          = if (length(line_replications))  min(line_replications)  else NA_integer_,
    max_replication          = if (length(line_replications))  max(line_replications)  else NA_integer_,
    mean_replication         = if (length(line_replications))  mean(line_replications) else NA_real_,
    min_sparse_replication   = if (length(sparse_replications)) min(sparse_replications)  else NA_integer_,
    max_sparse_replication   = if (length(sparse_replications)) max(sparse_replications)  else NA_integer_,
    mean_sparse_replication  = if (length(sparse_replications)) mean(sparse_replications) else NA_real_,
    min_common_replication   = if (length(common_replications)) min(common_replications)  else NA_integer_,
    max_common_replication   = if (length(common_replications)) max(common_replications)  else NA_integer_,
    mean_common_replication  = if (length(common_replications)) mean(common_replications) else NA_real_,
    n_groups                 = n_groups
  )
  
  list(
    allocation_matrix    = alloc,
    allocation_long      = allocation_long,
    overlap_matrix       = overlap_matrix,
    line_replications    = line_replications,
    environment_sizes    = environment_sizes,
    group_assignment     = group_assignment,
    group_by_environment = group_by_environment,
    group_overlap_matrix = group_overlap_matrix,
    summary              = summary_out,
    seed_used            = seed_used
  )
}
