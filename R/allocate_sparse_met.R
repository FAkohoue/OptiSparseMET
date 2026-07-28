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
#' The \code{"equireplicate"} method implements the M4 allocation of
#' Montesinos-Lopez et al. (2023): every non-common treatment appears in
#' exactly \eqn{r} environments (equal across-environment incidence). Here,
#' \eqn{r} counts different environments, not repeated plots within one
#' environment. Repeating a plot locally never increases connectivity. Equal
#' across-environment incidence distinguishes M4 from the random M3 method and
#' controls unequal information among treatments. Pairwise connectivity is a
#' separate property determined by how many distinct treatments each
#' environment pair shares. Equal incidence is preserved even
#' when environments have **unequal sizes** (partners offering different
#' amounts of space): the sparse loads \eqn{k^*_e} may differ across
#' environments, subject only to the resource identity
#' \eqn{\sum_e k^*_e = J^* \times r} and to the equal-replication degree sequence
#' being realisable (a Gale-Ryser condition checked internally). When all
#' environments are equal-sized this reduces to the classic
#' \eqn{J^* \times r = I \times k^*} identity. So a breeder can keep M4's
#' equal-replication strength while occupying exactly the space each location
#' offers.
#' 
#' `allow_approximate = FALSE` (the default) is the standard M4 path: the
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
#' \deqn{N = J \times r = \sum_{e=1}^{I} k_e}
#'
#' where \eqn{J} is the number of lines, \eqn{I} is the number of environments,
#' \eqn{k_e} is the number of lines in environment \eqn{e}, \eqn{r} is the
#' number of environments each line enters, and \eqn{N} is the total number of
#' line-by-environment observations. For equal environment sizes,
#' \eqn{\sum_e k_e = I k}. This identity makes the tradeoff between coverage
#' and replication depth explicit given fixed \eqn{N}.
#'
#' ## Common treatments
#'
#' Treatments in `common_treatments` are assigned to every environment before
#' sparse allocation begins and are excluded from the sparse allocation pool.
#' The supplied identities are fixed for this allocation: the function does not
#' reselect or replace them. “Common” fixes presence in every environment, not
#' equal local replication; site-specific replication is handled by the
#' downstream field-design plan.
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
#' For `"random_balanced"` and approximate `"equireplicate"`, phase one
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
#'   \item Exact requested environment sizes: each environment receives its
#'         requested sparse load \eqn{k_e^*}; loads may differ provided
#'         \eqn{\sum_e k_e^* = J^*r} and the degree sequence is realisable.
#' }
#' 
#' `allow_approximate = FALSE` (the default) enforces both conditions strictly.
#' If the total slot identity or the Gale-Ryser degree condition does not hold,
#' the function stops with an informative error. Use
#' [check_equireplicate_feasibility()] before calling. Construction assigns each
#' sparse treatment to the eligible environments with the largest remaining
#' capacity, guaranteeing exact loads and equal replication.
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
#'   handled separately within the within-environment design functions. IDs
#'   must be unique, non-missing, and non-empty.
#'
#' @param environments Character vector of environment names. Must contain at
#'   least two unique, non-missing, non-empty elements.
#'
#' @param allocation_method Character scalar. Sparse allocation strategy.
#'   Accepted values are `"random_balanced"` (M3-inspired stochastic
#'   allocation) and `"equireplicate"` (M4-type equal-replication allocation). The
#'   aliases `"M3"` and `"M4"` are also accepted and translated internally to
#'   their canonical names before any further processing.
#'
#' @param n_test_entries_per_environment Integer scalar or vector giving the
#'   field capacity (total test treatments, including common treatments) of each
#'   environment. A scalar applies uniformly. A vector may be **positional**
#'   (length = number of environments, in the same order) or **named** (matched
#'   to environment names, any order), which is the natural way to encode
#'   *site-specific* capacity when partners offer different amounts of space:
#'   the breeder simply occupies all the space each location offers, so
#'   environments need not be equal-sized. All values must be positive and at
#'   least `length(common_treatments)`. Variable capacities work with **both**
#'   methods: `"random_balanced"` fills each site to its capacity, and
#'   `"equireplicate"` (M4) still guarantees **equal replication** (each
#'   treatment in exactly `target_replications` environments) across
#'   unequal-sized environments whenever the resulting degree sequence is
#'   realisable (checked internally). For a location with *no* capacity limit,
#'   use [suggest_site_capacity()] to choose a plot number first, then pass it
#'   here. Use [suggest_safe_k()] to choose a value that guarantees full
#'   coverage.
#'
#' @param target_replications Optional positive integer. Target number of
#'   environments in which each non-common treatment should appear. For
#'   `"random_balanced"`, this is a soft target that guides the phase-two
#'   filling. For `"equireplicate"`, this is the strict replication level
#'   required for an exact balanced solution. If `NULL`, the function infers a
#'   value as `floor(total_sparse_slots / n_sparse_treatments)`, with a minimum
#'   of 1. This argument controls distinct-environment coverage only; it does
#'   not add plots within an environment.
#'
#' @param common_treatments Optional character vector of preselected treatment
#'   IDs. Supplied identities are fixed as present in all environments, placed
#'   before sparse allocation begins, and excluded from the sparse pool. This
#'   fixes presence, not equal site-specific replication. Values not present in
#'   `treatments` are silently dropped.
#'
#' @param seed_available Optional network-wide seed inventory. Supply either a
#'   named numeric vector or a data frame with columns `Treatment` and
#'   `SeedAvailable`. The inventory is shared across all environments; it is not
#'   reset independently at each site.
#'
#' @param seed_required_per_environment Seed consumed when one treatment is
#'   assigned to each environment: a positive scalar, a named numeric vector,
#'   or a data frame with columns `Environment` and `SeedRequiredPerPlot`.
#'   Values may include mandatory local replication (for example, twice the
#'   per-plot requirement for a two-replicate alpha design). Required with
#'   `seed_available`.
#'
#' @param minimum_seed_buffer Non-negative reserve retained for every treatment
#'   after mandatory allocation. Default 0.
#'
#' @param balance Character scalar. Post-construction refinement toward a
#'   near-balanced (regular-graph) design, preserving equal replication and
#'   every requested environment size. `"none"` (default) skips it; `"env_pair"`
#'   equalises the number of lines shared by each pair of environments
#'   (cross-environment connectivity); `"line_pair"` equalises pairwise
#'   treatment concurrence; `"both"` applies `"env_pair"` then `"line_pair"`.
#'   With genetic grouping, only within-group swaps are accepted, so
#'   group-by-environment counts remain unchanged. With seed constraints, swaps
#'   that would exceed either treatment's network inventory are rejected. A
#'   strict BIBD is generally unattainable when treatments greatly outnumber
#'   environments, so this targets achievable balance rather than constant
#'   \eqn{\lambda}.
#'
#' @param balance_iter Integer, default `2000`. Number of swap iterations per
#'   balance pass.
#'
#' @param balance_seed Optional integer seed for the balance swap search
#'   (independent of `seed`).
#'
#' @param Sigma_E Optional named environmental covariance matrix, or named list
#'   of covariance candidates. When supplied, a margin-preserving refinement
#'   allocates distinct shared treatments towards environment pairs according
#'   to their correlation-dependent targets. Local plot replication is not
#'   used in this connectivity calculation.
#' @param pair_target_se Positive target standard error used to convert each
#'   environment-pair correlation into a target number of distinct shared
#'   treatments.
#' @param pair_aggregate Pairwise protection rule for correlation-adaptive
#'   refinement: `"maximin"` (default), `"cvar"`, or `"mean"`.
#' @param pair_cvar_alpha Lower-tail fraction used when
#'   `pair_aggregate = "cvar"`.
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
#'   `allocation_method = "equireplicate"`, the total slot identity and
#'   equal-replication degree sequence must be feasible; the function stops
#'   otherwise. This is the standard M4 path and guarantees equal replication
#'   for every non-common treatment. When `TRUE`, those conditions are relaxed
#'   and minor replication imbalances are accepted; this is an exploratory
#'   fallback. For `"random_balanced"`, this argument has no effect.
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
#'   \item{`pairwise_connectivity`}{When `Sigma_E` is supplied, a data frame
#'     reporting the target and achieved number of distinct treatments shared
#'     by every environment pair. Repeated plots never increase this count.}
#'   \item{`pair_refinement_report`}{Before-and-after maximin, CVaR, or mean
#'     target-attainment scores and preserved-constraint indicators for the
#'     correlation-adaptive refinement.}
#'   \item{`seed_summary`}{When seed constraints are active, treatment-level
#'     seed available, allocated, buffered, remaining, and feasibility values;
#'     otherwise `NULL`.}
#'   \item{`summary`}{Named list with allocation metadata: `allocation_method`,
#'     `allocation_group_source`, `target_replications`, `n_treatments_total`,
#'     `n_sparse_treatments`, `n_common_treatments`, `total_sparse_slots`,
#'     `environment_sizes`, `min_replication`, `max_replication`,
#'     `mean_replication`, `min_sparse_replication`, `max_sparse_replication`,
#'     `mean_sparse_replication`, `min_common_replication`,
#'     `max_common_replication`, `mean_common_replication`, `pair_aggregate`,
#'     and `pair_target_se`.}
#'   \item{`seed_used`}{The integer seed passed to `set.seed()` internally, or
#'     `NULL` if no seed was supplied.}
#' }
#'
#' @seealso
#' [suggest_safe_k()] and [min_k_for_full_coverage()] for choosing a feasible
#' `n_test_entries_per_environment` before calling this function.
#' [check_equireplicate_feasibility()] for verifying the slot identity
#' before attempting a `"equireplicate"` allocation.
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
#'   allocation_method              = "equireplicate",
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
    allocation_method = c("random_balanced", "equireplicate", "M3", "M4"),
    n_test_entries_per_environment,
    target_replications = NULL,
    common_treatments = NULL,
    seed_available = NULL,
    seed_required_per_environment = NULL,
    minimum_seed_buffer = 0,
    balance = c("none", "env_pair", "line_pair", "both"),
    balance_iter = 2000L,
    balance_seed = NULL,
    Sigma_E = NULL,
    pair_target_se = 0.15,
    pair_aggregate = c("maximin", "cvar", "mean"),
    pair_cvar_alpha = 0.25,
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
  if (!is.null(seed_used)) {
    if (!is.numeric(seed_used) || length(seed_used) != 1L ||
        !is.finite(seed_used) || abs(seed_used - round(seed_used)) > 1e-8)
      stop("`seed` must be one finite integer or NULL.")
    seed_used <- as.integer(round(seed_used))
    set.seed(seed_used)
  }
  
  allocation_method <- match.arg(allocation_method)
  # "M3"/"M4" are convenience aliases for the two Montesinos-Lopez (2023)
  # methods. "equireplicate" is the canonical name for the M4-type constructor:
  # it guarantees equal treatment replication and exact requested environment
  # sizes, which may be unequal (NOT a strict BIBD, which generally cannot exist
  # when treatments greatly outnumber environments). Use `balance` to drive
  # concurrence toward a near-balanced design.
  if (allocation_method == "M3") allocation_method <- "random_balanced"
  if (allocation_method == "M4") allocation_method <- "equireplicate"

  balance <- match.arg(balance)
  pair_aggregate <- match.arg(pair_aggregate)

  allocation_group_source <- match.arg(allocation_group_source)
  group_method            <- match.arg(group_method)
  
  # ============================================================
  # 1. Basic validation and normalisation
  # ============================================================
  treatments <- as.character(treatments)
  environments <- as.character(environments)
  if (anyNA(treatments) || any(!nzchar(treatments)) ||
      anyNA(environments) || any(!nzchar(environments)))
    stop("Treatment and environment IDs must be non-missing and non-empty.")
  if (anyDuplicated(treatments))
    stop("`treatments` must contain unique IDs.")
  if (anyDuplicated(environments))
    stop("`environments` must contain unique names.")
  
  if (length(treatments) < 1L)   stop("`treatments` must contain at least one treatment ID.")
  if (length(environments) < 2L) stop("`environments` must contain at least two environment names.")
  
  n_treat <- length(treatments)
  n_env   <- length(environments)
  if (!is.numeric(pair_target_se) || length(pair_target_se) != 1L ||
      !is.finite(pair_target_se) || pair_target_se <= 0)
    stop("`pair_target_se` must be one finite positive value.")
  if (!is.numeric(pair_cvar_alpha) || length(pair_cvar_alpha) != 1L ||
      !is.finite(pair_cvar_alpha) || pair_cvar_alpha <= 0 ||
      pair_cvar_alpha > 1)
    stop("`pair_cvar_alpha` must be one value in (0, 1].")

  pair_targets <- NULL
  if (!is.null(Sigma_E)) {
    covariance_candidates <- if (
      is.list(Sigma_E) && !is.matrix(Sigma_E)
    ) Sigma_E else list(central = Sigma_E)
    if (!length(covariance_candidates))
      stop("`Sigma_E` cannot be an empty list.")
    target_by_candidate <- vapply(
      covariance_candidates,
      function(M) {
        M <- as.matrix(M)
        if (!is.numeric(M) || nrow(M) != ncol(M) ||
            is.null(rownames(M)) || is.null(colnames(M)) ||
            !all(environments %in% rownames(M)) ||
            !all(environments %in% colnames(M)))
          stop("Every `Sigma_E` candidate must be a named square matrix ",
               "covering all environments.")
        M <- M[environments, environments, drop = FALSE]
        if (any(!is.finite(M)) ||
            !isTRUE(all.equal(M, t(M), tolerance = 1e-8)) ||
            any(diag(M) <= 0))
          stop("Every `Sigma_E` candidate must be finite, symmetric, and ",
               "have a positive diagonal.")
        eigenvalues <- eigen(
          (M + t(M)) / 2, symmetric = TRUE, only.values = TRUE
        )$values
        if (min(eigenvalues) < -1e-8 *
            max(1, max(abs(eigenvalues))))
          stop("Every `Sigma_E` candidate must be positive semidefinite.")
        rho <- stats::cov2cor(M)[upper.tri(M)]
        pmax(2, ((1 - rho^2) / pair_target_se)^2)
      },
      numeric(n_env * (n_env - 1L) / 2L)
    )
    if (is.null(dim(target_by_candidate)))
      target_by_candidate <- matrix(
        target_by_candidate,
        nrow = n_env * (n_env - 1L) / 2L
      )
    pair_targets <- apply(target_by_candidate, 1L, max)
  }
  
  if (!is.numeric(n_test_entries_per_environment) ||
      any(!is.finite(n_test_entries_per_environment)) ||
      any(abs(n_test_entries_per_environment -
                round(n_test_entries_per_environment)) > 1e-8))
    stop("`n_test_entries_per_environment` must contain finite integers.")
  if (length(n_test_entries_per_environment) == 1L) {
    k_vec <- rep(as.integer(round(n_test_entries_per_environment)), n_env)
  } else if (!is.null(names(n_test_entries_per_environment))) {
    # Named per-site capacities: align to the environment order so a breeder can
    # supply "as much space as each partner offers" per location, in any order.
    miss <- setdiff(environments, names(n_test_entries_per_environment))
    if (length(miss))
      stop("`n_test_entries_per_environment` is missing capacities for ",
           "environment(s): ", paste(miss, collapse = ", "), ".")
    k_vec <- as.integer(round(n_test_entries_per_environment[environments]))
  } else {
    k_vec <- as.integer(round(n_test_entries_per_environment))
  }

  if (length(k_vec) != n_env)
    stop("`n_test_entries_per_environment` must have length 1, length(environments), ",
         "or be a named vector covering every environment.")
  if (any(is.na(k_vec)) || any(k_vec < 1L))
    stop("All values of `n_test_entries_per_environment` must be positive integers. ",
         "For a site with no capacity limit, use suggest_site_capacity() to pick a ",
         "plot number first.")
  if (any(k_vec > n_treat))
    stop("No environment capacity may exceed the number of unique treatments.")
  if (!is.logical(allow_approximate) || length(allow_approximate) != 1L ||
      is.na(allow_approximate) ||
      !is.logical(balance_groups_across_env) ||
      length(balance_groups_across_env) != 1L ||
      is.na(balance_groups_across_env) ||
      !is.logical(force_group_connectivity) ||
      length(force_group_connectivity) != 1L ||
      is.na(force_group_connectivity))
    stop("Logical control arguments must be single non-missing values.")
  if (!is.numeric(balance_iter) || length(balance_iter) != 1L ||
      !is.finite(balance_iter) || balance_iter < 0 ||
      abs(balance_iter - round(balance_iter)) > 1e-8)
    stop("`balance_iter` must be one non-negative integer.")
  balance_iter <- as.integer(round(balance_iter))
  if (!is.null(balance_seed) &&
      (!is.numeric(balance_seed) || length(balance_seed) != 1L ||
       !is.finite(balance_seed) ||
       abs(balance_seed - round(balance_seed)) > 1e-8))
    stop("`balance_seed` must be one finite integer or NULL.")
  
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

  # Optional network-wide seed constraint. `seed_cost[e]` is the mandatory
  # amount consumed when a treatment is assigned to environment e. Tracking one
  # shared balance prevents the same inventory from being spent independently
  # at several sites.
  seed_constrained <- !is.null(seed_available) ||
    !is.null(seed_required_per_environment)
  seed_budget <- seed_cost <- seed_consumed <- NULL
  if (!is.numeric(minimum_seed_buffer) || length(minimum_seed_buffer) != 1L ||
      !is.finite(minimum_seed_buffer) || minimum_seed_buffer < 0)
    stop("`minimum_seed_buffer` must be a finite non-negative scalar.")
  if (!seed_constrained && minimum_seed_buffer > 0)
    stop("Seed inputs are required when `minimum_seed_buffer` is positive.")
  if (seed_constrained) {
    if (is.null(seed_available) || is.null(seed_required_per_environment))
      stop("Supply both `seed_available` and `seed_required_per_environment`.")
    if (allocation_method == "equireplicate")
      stop("Seed-constrained allocation currently supports `random_balanced`; ",
           "use that method when environment seed costs differ.")
    if (is.data.frame(seed_available)) {
      if (!all(c("Treatment", "SeedAvailable") %in% names(seed_available)))
        stop("Seed inventory data must contain `Treatment` and `SeedAvailable`.")
      if (anyDuplicated(as.character(seed_available$Treatment)))
        stop("Seed inventory contains duplicate treatment IDs.")
      seed_budget <- stats::setNames(
        as.numeric(seed_available$SeedAvailable),
        as.character(seed_available$Treatment))
    } else {
      seed_budget <- as.numeric(seed_available)
      names(seed_budget) <- names(seed_available)
    }
    if (is.null(names(seed_budget)) || any(names(seed_budget) == "") ||
        anyDuplicated(names(seed_budget)) ||
        !all(treatments %in% names(seed_budget)))
      stop("`seed_available` must be named and cover every treatment.")
    seed_budget <- seed_budget[treatments] - minimum_seed_buffer
    if (any(!is.finite(seed_budget)) || any(seed_budget < 0))
      stop("Seed availability must be finite and at least the requested buffer.")

    if (is.data.frame(seed_required_per_environment)) {
      if (!all(c("Environment", "SeedRequiredPerPlot") %in%
               names(seed_required_per_environment)))
        stop("Seed-cost data must contain `Environment` and ",
             "`SeedRequiredPerPlot`.")
      if (anyDuplicated(
            as.character(seed_required_per_environment$Environment)))
        stop("Seed-cost data contain duplicate environments.")
      seed_cost <- stats::setNames(
        as.numeric(seed_required_per_environment$SeedRequiredPerPlot),
        as.character(seed_required_per_environment$Environment))
    } else if (length(seed_required_per_environment) == 1L) {
      seed_cost <- stats::setNames(
        rep(as.numeric(seed_required_per_environment), n_env), environments)
    } else {
      seed_cost <- as.numeric(seed_required_per_environment)
      names(seed_cost) <- names(seed_required_per_environment)
    }
    if (is.null(names(seed_cost)) || any(names(seed_cost) == "") ||
        anyDuplicated(names(seed_cost)) ||
        !all(environments %in% names(seed_cost)))
      stop("`seed_required_per_environment` must cover every environment.")
    seed_cost <- seed_cost[environments]
    if (any(!is.finite(seed_cost)) || any(seed_cost <= 0))
      stop("Environment seed costs must be finite and positive.")

    seed_consumed <- stats::setNames(numeric(n_treat), treatments)
    common_cost <- sum(seed_cost)
    if (n_common > 0L) {
      bad <- common_treatments[
        seed_budget[common_treatments] + 1e-8 < common_cost]
      if (length(bad))
        stop("Common treatments lack seed for mandatory placement at all sites: ",
             paste(bad, collapse = ", "), ". Required per common treatment: ",
             signif(common_cost, 6), ".")
      seed_consumed[common_treatments] <- common_cost
    }
  }
  
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
        !is.finite(target_replications) || target_replications < 1L ||
        abs(target_replications - round(target_replications)) > 1e-8)
      stop("`target_replications` must be NULL or a single positive integer.")
    target_replications <- as.integer(round(target_replications))
  }
  
  # ============================================================
  # 4. Equireplicate (M4) feasibility checks
  # ============================================================
  if (allocation_method == "equireplicate") {

    required_slots <- n_sparse * target_replications

    if (!allow_approximate && required_slots != total_sparse_slots) {
      stop(paste0(
        "Exact equireplicate allocation infeasible. ",
        "Required sparse slots = ", required_slots,
        ", available sparse slots = ", total_sparse_slots, ". ",
        "Adjust n_test_entries_per_environment, target_replications, or set ",
        "allow_approximate = TRUE."
      ))
    }

    # Equireplicate = equal REPLICATION (each sparse treatment in exactly r
    # environments), which is M4's strength. That does NOT require equal
    # environment sizes: unequal per-site capacities are allowed as long as the
    # equal-replication degree sequence is realizable (Gale-Ryser). The strict
    # constructor below assigns each treatment to its r environments by largest
    # remaining capacity, which realizes any feasible sequence.
    if (!allow_approximate && n_sparse > 0L &&
        !.equireplicate_degree_feasible(k_sparse, target_replications, n_sparse)) {
      stop(paste0(
        "Equireplicate with these per-environment capacities is infeasible: ",
        "the equal-replication degree sequence (each of ", n_sparse,
        " sparse treatments in exactly ", target_replications,
        " environments, environment loads = ", paste(k_sparse, collapse = ", "),
        ") is not realizable. Rebalance the per-site capacities, change ",
        "target_replications, or use allocation_method = \"random_balanced\"."
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
  if (allocation_method == "equireplicate" && !allow_approximate) {
    
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
        
        # Gale-Ryser greedy: prefer the environments with the largest remaining
        # capacity so unequal environment sizes are filled without stranding a
        # large environment at the end.
        remaining <- target_env_load[candidate_envs] - sparse_env_load[candidate_envs]
        if (n_groups > 0L) {
          grp         <- sparse_groups[[trt]]
          grp_presence <- group_env_current[grp, candidate_envs]
          pref        <- ifelse(grp_presence == 0L, 1L, 0L)
          # Havel-Hakimi validity requires using a largest remaining column
          # degree. Genetic-group dispersion may break ties, but must not take
          # precedence over degree realization in strict mode.
          ord         <- order(-remaining, -pref, candidate_envs)
          chosen_env  <- candidate_envs[ord][1L]
          group_env_current[grp, chosen_env] <- 1L
        } else {
          max_rem    <- max(remaining)
          best_envs  <- candidate_envs[remaining == max_rem]
          chosen_env <- if (length(best_envs) == 1L) best_envs
                        else sample(best_envs, 1L)
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
      
      # Phase 1: force minimum coverage. Under seed constraints, allocate lines
      # with the fewest affordable environments first; otherwise a flexible,
      # well-supplied line can occupy the only slot available to a constrained
      # line and make a feasible network appear infeasible.
      coverage_order <- if (seed_constrained) {
        feasible_count <- vapply(
          sparse_treatments,
          function(trt) sum(seed_cost <= seed_budget[trt] + 1e-8),
          integer(1))
        if (any(feasible_count == 0L))
          stop("At least one treatment cannot afford a plot in any environment.")
        tie_break <- stats::runif(n_sparse)
        sparse_treatments[
          order(feasible_count, seed_budget[sparse_treatments], tie_break)]
      } else {
        sample(sparse_treatments)
      }
      for (trt in coverage_order) {
        spare          <- k_vec - env_load
        candidate_envs <- environments[spare > 0L]
        if (seed_constrained)
          candidate_envs <- candidate_envs[
            seed_consumed[trt] + seed_cost[candidate_envs] <=
              seed_budget[trt] + 1e-8
          ]
        
        if (length(candidate_envs) == 0L)
          stop("Cannot give treatment `", trt,
               "` minimum coverage under the joint field-capacity and seed ",
               "constraints. Increase seed, reduce site capacities/mandatory ",
               "replication, or change the common set.")
        
        if (n_groups > 0L) {
          grp          <- sparse_groups[[trt]]
          grp_presence <- group_env_current[grp, candidate_envs]
          pref         <- ifelse(grp_presence == 0L, 1L, 0L)
          load_fraction <- env_load[candidate_envs] /
            k_vec[match(candidate_envs, environments)]
          ord <- if (seed_constrained)
            order(-pref, load_fraction, seed_cost[candidate_envs],
                  candidate_envs)
          else order(-pref, load_fraction, candidate_envs)
          chosen_env   <- candidate_envs[ord][1L]
          group_env_current[grp, chosen_env] <- 1L
        } else {
          load_fraction <- env_load[candidate_envs] /
            k_vec[match(candidate_envs, environments)]
          ord <- if (seed_constrained)
            order(load_fraction, seed_cost[candidate_envs], candidate_envs)
          else order(load_fraction, candidate_envs)
          chosen_env <- candidate_envs[ord][1L]
        }
        
        alloc[trt, chosen_env] <- 1L
        env_load[chosen_env]   <- env_load[chosen_env] + 1L
        line_rep[trt]          <- 1L
        if (seed_constrained)
          seed_consumed[trt] <- seed_consumed[trt] + seed_cost[chosen_env]
      }
      
    } else {
      line_rep          <- stats::setNames(integer(0), character(0))
      group_env_current <- NULL
    }
    
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
    
    # Phase 2: fill remaining capacity. Re-evaluate the most constrained site
    # after every assignment so an early site cannot consume the only lines
    # capable of filling a later, more seed-demanding environment.
    while (any(env_load < k_vec)) {
      active <- environments[env_load < k_vec]
      feasible_by_env <- lapply(active, function(env_name) {
        cand <- sparse_treatments[alloc[sparse_treatments, env_name] == 0L]
        if (seed_constrained)
          cand <- cand[
            seed_consumed[cand] + seed_cost[env_name] <=
              seed_budget[cand] + 1e-8
          ]
        cand
      })
      names(feasible_by_env) <- active
      open_slots <- k_vec[match(active, environments)] - env_load[active]
      feasible_n <- lengths(feasible_by_env)
      if (any(feasible_n == 0L)) {
        failed <- active[which(feasible_n == 0L)[1L]]
        stop("Seed-constrained allocation could not fill environment `", failed,
             "` (filled ", env_load[failed], " of ",
             k_vec[match(failed, environments)],
             " entries). The requested capacities are not feasible with the ",
             "available per-treatment seed budgets.")
      }
      slack <- feasible_n - open_slots
      most_constrained <- active[slack == min(slack)]
      if (seed_constrained && length(most_constrained) > 1L) {
        most_constrained <- most_constrained[
          seed_cost[most_constrained] == max(seed_cost[most_constrained])]
      }
      env_name <- most_constrained[1L]
      candidates <- feasible_by_env[[env_name]]
        
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
        bonus_weight     <- if (n_env_groups_now <
                               min(min_groups_per_environment, n_groups))
          100 else 5
        score <- score + bonus_weight * lacking_bonus
      }
        
      if (n_groups > 0L && force_group_connectivity &&
          !is.null(min_env_per_group)) {
        env_count_group <- rowSums(group_env_current > 0L)
        cand_groups     <- sparse_groups[candidates]
        conn_term       <- 50 * pmax(0, as.numeric(
          min_env_per_group - env_count_group[cand_groups]
        ))
        conn_term[is.na(conn_term)] <- 0
        score <- score + conn_term
      }
        
      if (allocation_method == "equireplicate" && allow_approximate)
        score[deficit[names(score)] < 0] <-
          score[deficit[names(score)] < 0] - 1000
      if (seed_constrained) {
        headroom <- (seed_budget[candidates] - seed_consumed[candidates] -
                       seed_cost[env_name]) / pmax(seed_budget[candidates], 1)
        score <- score + 1e-6 * headroom
      }
        
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
      if (seed_constrained)
        seed_consumed[chosen] <- seed_consumed[chosen] + seed_cost[env_name]
        
      if (n_groups > 0L) {
        gp                              <- sparse_groups[[chosen]]
        group_current[gp]               <- group_current[gp] + 1L
        group_env_current[gp, env_name] <- 1L
      }
    }
  }
  
  # ============================================================
  # 7C. Near-balanced concurrence / connectivity refinement (P1.2/P1.3)
  # ============================================================
  # Swap-based improvement that drives the design toward a near-balanced
  # (regular-graph) structure while preserving the equireplication and equal
  # environment-load margins. Swaps exchange a sparse treatment between two
  # environments, so every treatment keeps its replication count and every
  # environment keeps its size. "env_pair" equalises the number of lines shared
  # by each pair of environments (the connectivity property MET inference depends
  # on); "line_pair" equalises pairwise treatment concurrence; "both" runs
  # env_pair first, then line_pair.
  balance_report <- NULL
  if (balance != "none" && n_sparse >= 2L) {
    if (!is.null(balance_seed)) set.seed(balance_seed)
    S0 <- alloc[sparse_treatments, , drop = FALSE]
    storage.mode(S0) <- "integer"
    kinds <- switch(balance,
                    env_pair  = "env_pair",
                    line_pair = "line_pair",
                    both      = c("env_pair", "line_pair"))
    before <- .balance_metrics(S0)
    S1 <- S0
    for (kd in kinds)
      S1 <- .balance_allocation(
        S1, kind = kd, iter = balance_iter,
        groups = if (n_groups > 0L) sparse_groups else NULL,
        seed_budget = if (seed_constrained)
          seed_budget[sparse_treatments] else NULL,
        environment_cost = if (seed_constrained) seed_cost else NULL)
    after <- .balance_metrics(S1)
    alloc[sparse_treatments, ] <- S1
    if (seed_constrained)
      seed_consumed[sparse_treatments] <-
        as.numeric(S1 %*% seed_cost[colnames(S1)])
    balance_report <- list(
      balance = balance, before = before, after = after,
      constraints_preserved = c(
        margins = TRUE, groups = n_groups > 0L,
        seed_budget = seed_constrained))
  }

  # Correlation-adaptive refinement changes only binary treatment membership.
  # It never adds plots or treats repeated plots of one genotype as additional
  # connectivity.
  pair_refinement_report <- NULL
  if (!is.null(pair_targets)) {
    if (!is.null(balance_seed)) set.seed(balance_seed)
    S0 <- alloc[sparse_treatments, , drop = FALSE]
    storage.mode(S0) <- "integer"
    refined <- .adaptive_pair_refinement(
      S0,
      pair_targets = pair_targets,
      common_offset = n_common,
      aggregate = pair_aggregate,
      cvar_alpha = pair_cvar_alpha,
      iter = balance_iter,
      groups = if (n_groups > 0L) sparse_groups else NULL,
      seed_budget = if (seed_constrained)
        seed_budget[sparse_treatments] else NULL,
      environment_cost = if (seed_constrained) seed_cost else NULL
    )
    alloc[sparse_treatments, ] <- refined$allocation
    if (seed_constrained)
      seed_consumed[sparse_treatments] <-
        as.numeric(refined$allocation %*%
                     seed_cost[colnames(refined$allocation)])
    pair_refinement_report <- list(
      aggregate = pair_aggregate,
      cvar_alpha = pair_cvar_alpha,
      target_se = pair_target_se,
      before_score = refined$before_score,
      after_score = refined$after_score,
      constraints_preserved = c(
        treatment_margins = TRUE,
        environment_margins = TRUE,
        groups = n_groups > 0L,
        seed_budget = seed_constrained
      )
    )
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
  pairwise_connectivity <- NULL
  if (!is.null(pair_targets)) {
    pair_index <- which(upper.tri(overlap_matrix), arr.ind = TRUE)
    achieved <- overlap_matrix[upper.tri(overlap_matrix)]
    pairwise_connectivity <- data.frame(
      Environment1 = rownames(overlap_matrix)[pair_index[, 1L]],
      Environment2 = colnames(overlap_matrix)[pair_index[, 2L]],
      TargetDistinctSharedTreatments = as.numeric(pair_targets),
      AchievedDistinctSharedTreatments = as.numeric(achieved),
      TargetFraction = pmin(1, as.numeric(achieved / pair_targets)),
      stringsAsFactors = FALSE
    )
  }
  
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
    n_groups                 = n_groups,
    pair_aggregate           = if (!is.null(pair_targets))
      pair_aggregate else NA_character_,
    pair_target_se           = if (!is.null(pair_targets))
      pair_target_se else NA_real_
  )

  seed_summary <- if (seed_constrained) data.frame(
    Treatment = treatments,
    SeedAvailable = seed_budget[treatments] + minimum_seed_buffer,
    SeedBuffer = minimum_seed_buffer,
    SeedAllocated = seed_consumed[treatments],
    SeedRemaining = seed_budget[treatments] + minimum_seed_buffer -
      seed_consumed[treatments],
    SeedSpendableRemaining = seed_budget[treatments] -
      seed_consumed[treatments],
    Feasible = seed_consumed[treatments] <= seed_budget[treatments] + 1e-8,
    stringsAsFactors = FALSE
  ) else NULL
  
  list(
    allocation_matrix    = alloc,
    allocation_long      = allocation_long,
    overlap_matrix       = overlap_matrix,
    line_replications    = line_replications,
    environment_sizes    = environment_sizes,
    group_assignment     = group_assignment,
    group_by_environment = group_by_environment,
    group_overlap_matrix = group_overlap_matrix,
    balance_report       = balance_report,
    pair_refinement_report = pair_refinement_report,
    pairwise_connectivity = pairwise_connectivity,
    seed_summary         = seed_summary,
    summary              = summary_out,
    seed_used            = seed_used
  )
}


# Gale-Ryser feasibility for an equireplicate design with (possibly unequal)
# per-environment sparse capacities. The bipartite degree sequence has every
# sparse treatment with degree r (in exactly r environments) and environment e
# with degree k_sparse[e]. It is realizable iff the totals match and, for the
# environment loads sorted in decreasing order, every prefix sum of the top k
# loads does not exceed n_sparse * min(r, k).
.equireplicate_degree_feasible <- function(k_sparse, r, n_sparse) {
  r <- as.integer(r)
  if (r < 1L || r > length(k_sparse)) return(FALSE)
  if (any(k_sparse < 0L) || any(k_sparse > n_sparse)) return(FALSE)
  if (sum(k_sparse) != n_sparse * r) return(FALSE)
  a  <- sort(as.integer(k_sparse), decreasing = TRUE)
  cs <- cumsum(a)
  for (k in seq_along(a)) {
    if (cs[k] > n_sparse * min(r, k)) return(FALSE)
  }
  TRUE
}
