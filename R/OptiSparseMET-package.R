#' OptiSparseMET: Sparse Multi-Environment Trial Design with Flexible Local Field Layout
#'
#' `OptiSparseMET` constructs sparse multi-environment trial designs by
#' addressing treatment allocation across environments and field layout
#' construction within environments as a two-level problem, solved jointly
#' under shared statistical and genetic constraints. The package is intended
#' for breeding programs where the number of candidate lines exceeds what any
#' single environment can accommodate, environments are heterogeneous in
#' capacity and genetic composition, and seed availability imposes practical
#' limits on replication that purely theoretical designs ignore.
#'
#' @details
#' ## Stage 1: across-environment sparse allocation
#'
#' `allocate_sparse_met()` determines which treatments enter which
#' environments. Two allocation strategies are available.
#' `"random_balanced"` implements an M3-type stochastic allocation that
#' approximates balance without requiring exact BIBD parameters to be
#' satisfiable -- appropriate when environment capacities differ or when the
#' trial dimensions do not admit an exact balanced solution.
#' `"balanced_incomplete"` implements an M4-type allocation following BIBD
#' principles at the MET level, enforcing equal replication and approximately
#' uniform pairwise co-occurrence -- appropriate when environments are
#' comparable in size and equal replication is a hard requirement.
#'
#' Allocation can be guided by genetic structure through family labels, genomic
#' relationship clusters (GRM), or pedigree relationship clusters (A), ensuring
#' that genetic groups are distributed across environments rather than
#' concentrated within a subset of them. Common treatments can be forced into
#' all environments to establish direct cross-environment connectivity
#' independent of model assumptions. Feasibility of an exact balanced
#' incomplete allocation can be checked in advance with
#' `check_balanced_incomplete_feasibility()`.
#'
#' The underlying resource identity is:
#'
#' \deqn{N = J \times r = I \times k}
#'
#' where \eqn{J} is the number of treatments, \eqn{r} is the number of
#' environments each treatment enters, \eqn{I} is the number of environments,
#' and \eqn{k} is the number of treatments per environment. This identity makes
#' the tradeoff between coverage and replication depth explicit given fixed
#' total resources \eqn{N}.
#'
#' ## Stage 2: within-environment field design
#'
#' Each environment is designed independently after allocation.
#' `met_prep_famoptg()` constructs augmented, partially replicated (p-rep), or
#' RCBD-type repeated-check designs for environments where block-based local
#' control is appropriate. `met_alpha_rc_stream()` constructs stream-based
#' alpha row-column designs for environments with fixed grid geometry and
#' continuous field-book order, accommodating AR1 and AR1xAR1 spatial models
#' at the analysis stage.
#'
#' Note: `met_prep_famoptg()` and `met_alpha_rc_stream()` are the
#' OptiSparseMET-specific versions of the constructors from the OptiDesign
#' package. They carry the `met_` prefix to avoid namespace conflicts when
#' both packages are loaded simultaneously. Their argument lists and return
#' values are identical to `prep_famoptg()` and `alpha_rc_stream()` in
#' OptiDesign.
#'
#' When replication within an environment is flexible,
#' `assign_replication_by_seed()` evaluates seed availability against the
#' per-plot seed requirement for each treatment and assigns replication roles
#' -- replicated, unreplicated, or excluded -- before the field layout is
#' constructed. Both within-environment design functions accept optional
#' relationship matrices for structure-aware entry arrangement and dispersion
#' optimisation, and can evaluate the resulting design under a mixed-model
#' framework prior to field deployment.
#'
#' Environment-level field books are assembled into a single MET field book
#' with `combine_met_fieldbooks()`, which handles heterogeneous column sets
#' across environments and adds MET-level metadata to every row.
#' The complete two-stage pipeline can also be run end-to-end through
#' `plan_sparse_met_design()`.
#'
#' @section Exported functions:
#'
#' | Function | Role |
#' |---|---|
#' | `allocate_sparse_met()` | Allocate treatments across environments |
#' | `check_balanced_incomplete_feasibility()` | Verify slot identity before allocation |
#' | `derive_allocation_groups()` | Derive genetic group labels for allocation |
#' | `min_k_for_full_coverage()` | Minimum per-environment capacity for full coverage |
#' | `suggest_safe_k()` | Suggest a safe per-environment capacity |
#' | `warn_if_k_too_small()` | Non-fatal capacity pre-flight check |
#' | `assign_replication_by_seed()` | Classify treatments into replication roles |
#' | `met_prep_famoptg()` | Build block-based within-environment layouts |
#' | `met_alpha_rc_stream()` | Build stream-based row-column layouts |
#' | `met_evaluate_famoptg_efficiency()` | Evaluate efficiency of met_prep_famoptg() designs |
#' | `met_evaluate_alpha_efficiency()` | Evaluate efficiency of met_alpha_rc_stream() designs |
#' | `met_optimize_famoptg()` | Optimise met_prep_famoptg() designs (Random Restart) |
#' | `met_optimize_alpha_rc()` | Optimise met_alpha_rc_stream() designs (RS/SA/GA) |
#' | `combine_met_fieldbooks()` | Stack environment field books into one MET field book |
#' | `plan_sparse_met_design()` | Run the full two-stage pipeline |
#'
#' @section Naming convention:
#' The `met_` prefix on all six within-environment design functions
#' (`met_prep_famoptg()`, `met_alpha_rc_stream()`,
#' `met_evaluate_famoptg_efficiency()`, `met_evaluate_alpha_efficiency()`,
#' `met_optimize_famoptg()`, `met_optimize_alpha_rc()`) distinguishes
#' them from their counterparts in the OptiDesign package. Both packages can
#' be loaded simultaneously without namespace conflicts.
#'
#' @section Example data:
#' `OptiSparseMET_example_data` is a synthetic dataset containing treatment
#' vectors, environment names, family metadata, seed inventories, relationship
#' matrices (GRM, A, K), and pre-built argument lists for allocation and local
#' design functions. All components refer to the same set of synthetic lines
#' and environments and can be passed to any package function without
#' modification. Load with:
#'
#' ```r
#' data("OptiSparseMET_example_data", package = "OptiSparseMET")
#' ```
#'
#' @section Dependencies:
#' **Matrix** is used for sparse matrix operations in efficiency evaluation and
#' design construction. **pracma** provides `mod()` for serpentine traversal
#' logic in `met_alpha_rc_stream()` and `met_prep_famoptg()`.
#'
#' @references
#' Montesinos-Lopez, O. A., Mosqueda-Gonzalez, B. A., Salinas-Ruiz, J.,
#' Montesinos-Lopez, A., & Crossa, J. (2023). Sparse multi-trait genomic
#' prediction under balanced incomplete block design. *The Plant Genome*,
#' 16, e20305. \doi{10.1002/tpg2.20305}
#'
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P.,
#' ..., & Moreau, L. (2012). Maximizing the reliability of genomic selection
#' by optimizing the calibration set of reference individuals. *Genetics*,
#' 192(2), 715-728. \doi{10.1534/genetics.112.141473}
#'
#' @keywords internal
"_PACKAGE"
