#' OptiSparseMET: Sparse Multi-Environment Trial Design with Flexible Local Field Layout
#'
#' `OptiSparseMET` constructs sparse multi-environment trial (MET) designs by
#' integrating treatment allocation across environments with field layout
#' within environments. It is intended for breeding programmes in which
#' candidate numbers exceed local field capacity, environments differ in
#' relevance and precision, and finite seed stocks constrain replication. The
#' workflow addresses three linked decisions: (i) which environments and
#' treatments to test; (ii) how to connect and replicate those treatments
#' across the trial network; and (iii) how to arrange their plots within each
#' field.
#'
#' @details
#' ## Stage 1: across-environment sparse allocation
#'
#' `allocate_sparse_met()` determines which treatments enter which
#' environments. Two allocation strategies are available.
#' `"random_balanced"` implements an M3-type stochastic allocation that
#' approximates balance without requiring exact balanced incomplete block
#' design (BIBD) parameters. It is appropriate when environment capacities
#' differ or when the trial dimensions do not admit an exact balanced solution.
#' `"equireplicate"` implements an M4-type allocation at the MET level,
#' enforcing equal distinct-environment incidence across the network while
#' allowing unequal environment sizes when the requested row and column degrees are
#' jointly realisable. Optional balancing improves pairwise co-occurrence
#' without changing either margin. It is not a strict BIBD, which generally
#' cannot exist when treatments greatly outnumber environments.
#'
#' Allocation can be guided by family labels, a genomic relationship matrix
#' (GRM), or a pedigree numerator relationship matrix (A). Additive
#' relationships are the default; dominance is optional and is intended for
#' hybrids. Common treatments establish direct cross-environment connectivity
#' independent of an assumed covariance model. Their number and identities can
#' be selected with `optimize_common_treatments()`; local plot replication is
#' assigned afterwards from remaining seed. Feasibility of an exact equireplicate
#' allocation can be checked before construction with
#' `check_equireplicate_feasibility()`.
#'
#' When seed limits are supplied, allocation draws from one network-wide
#' inventory. A treatment consumes seed only at environments where it is
#' present, but its consumption is accumulated across all such environments.
#' An optional per-treatment reserve is protected; its default is zero grams.
#'
#' The underlying resource identity is:
#'
#' \deqn{N = J \times r = \sum_{e=1}^{I} k_e}
#'
#' where \eqn{J} is the number of treatments, \eqn{r} is the number of
#' environments each treatment enters, \eqn{I} is the number of environments,
#' and \eqn{k_e} is the number of treatments in environment \eqn{e}. Equal
#' site size is the special case \eqn{k_e = k}. This identity makes the trade-off
#' between coverage and replication depth explicit given fixed total resources
#' \eqn{N}.
#'
#' ## Environmental evidence and uncertainty
#'
#' Weather, soil, management, and geography are quality-controlled and retained
#' as separate environmental kernels. Weather-by-soil,
#' weather-by-management, and soil-by-management kernels represent explicit
#' cross-modality hypotheses; they are not interpreted as evidence that the
#' corresponding effects are non-zero. Dedicated variable-level hypotheses
#' can be constructed with `build_variable_interaction_kernels()`. Their tensor
#' features are residualised against their exact parent features, and strong
#' heredity retains the parent-variable kernels. When historical MET responses are
#' available, `calibrate_environment_covariance()` estimates non-negative
#' kernel contributions by blocked validation.
#' `assess_variable_interactions()` evaluates dedicated components separately
#' with multiplicity control before their joint central-covariance assessment.
#' When those responses are
#' unavailable, no user-defined modality weights are required: identity is the
#' central genetic environment covariance and the unweighted kernels define
#' separate maximin sensitivity scenarios.
#'
#' `infer_mega_environments()` estimates the number of descriptive
#' environmental groups by stability-screened hierarchical and
#' kernel-principal-coordinate clustering. It reports a provisional or stable
#' partition only when separation, group size, and reproducibility requirements
#' are met; otherwise it returns one unpartitioned network. Descriptive
#' environmental groups are not labelled genetic mega-environments without
#' historical response evidence.
#'
#' ## Stage 2: within-environment field design
#'
#' Each environment is designed independently after allocation.
#' `met_prep_famoptg()` constructs augmented, partially replicated (p-rep), or
#' randomised complete block design (RCBD)-type repeated-check layouts for
#' environments where block-based local control is appropriate.
#' `met_alpha_rc_stream()` constructs stream-based alpha row-column designs for
#' environments with fixed grid geometry and continuous field-book order,
#' accommodating first-order autoregressive spatial models in one dimension
#' (AR1) or two dimensions (AR1xAR1) at the analysis stage.
#'
#' Note: `met_prep_famoptg()` and `met_alpha_rc_stream()` are the
#' OptiSparseMET-specific versions of the constructors from the OptiDesign
#' package. They carry the `met_` prefix to avoid namespace conflicts when
#' both packages are loaded simultaneously. Their argument lists and return
#' values are identical to `prep_famoptg()` and `alpha_rc_stream()` in
#' OptiDesign.
#'
#' When replication within an environment is flexible, the remaining
#' network-wide seed inventory is evaluated against the local per-plot
#' requirement before assigning replicated, unreplicated, or excluded roles.
#' `optimize_common_treatments()` selects common presence only: every selected
#' anchor occurs in every environment and counts once for pairwise
#' connectivity. Local plot replication is assigned afterwards from the
#' remaining seed inventory and never increases connectivity.
#' Both local design functions accept optional relationship matrices for
#' structure-aware entry arrangement and dispersion optimisation. Statistical
#' information, simulation, and cost summaries use delivered plot counts,
#' including repeated checks and partial replication, rather than incidence
#' alone.
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
#' | `check_equireplicate_feasibility()` | Verify the equireplicate degree sequence before allocation |
#' | `derive_allocation_groups()` | Derive genetic group labels for allocation |
#' | `min_k_for_full_coverage()` | Minimum per-environment capacity for full coverage |
#' | `suggest_safe_k()` | Suggest a safe per-environment capacity |
#' | `warn_if_k_too_small()` | Non-fatal capacity pre-flight check |
#' | `assign_replication_by_seed()` | Classify treatments into replication roles |
#' | `optimize_common_treatments()` | Robustly select the size and identities of the global common core |
#' | `met_prep_famoptg()` | Build block-based within-environment layouts |
#' | `met_alpha_rc_stream()` | Build stream-based row-column layouts |
#' | `met_evaluate_famoptg_efficiency()` | Evaluate efficiency of met_prep_famoptg() designs |
#' | `met_evaluate_alpha_efficiency()` | Evaluate efficiency of met_alpha_rc_stream() designs |
#' | `build_environment_kernels()` | Construct separate environmental main and interaction kernels |
#' | `build_variable_interaction_kernels()` | Construct residualised, auditable kernels for named variable interactions |
#' | `assess_variable_interactions()` | Assess dedicated interactions by blocked validation and multiplicity control |
#' | `calibrate_environment_covariance()` | Calibrate genetic environment covariance from historical MET responses |
#' | `infer_mega_environments()` | Infer stable descriptive environmental groups |
#' | `met_information()` | Evaluate coupled MET prediction-error variance and mean coefficient of determination (CDmean) |
#' | `simulate_met()` | Estimate selection accuracy and gain with Monte Carlo uncertainty |
#' | `make_met_benchmark_designs()` | Construct full, random, M3, and feasible strict-M4 comparators |
#' | `benchmark_met_designs()` | Compare paired performance and operational feasibility across stress scenarios |
#' | `benchmark_environment_models()` | Compare separate and combined environmental covariance hypotheses |
#' | `simulate_environment_model_benchmark()` | Quantify recovery of a known environmental covariance |
#' | `benchmark_environment_missingness()` | Stress-test environmental kernels after covariate masking and imputation |
#' | `summarize_design_stability()` | Quantify common-set, capacity, and allocation stability |
#' | `met_optimize_famoptg()` | Optimise met_prep_famoptg() designs by random restart |
#' | `met_optimize_alpha_rc()` | Optimise met_alpha_rc_stream() designs by random restart, simulated annealing, or genetic algorithm |
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
#' **Matrix** supports sparse matrix operations in efficiency evaluation and
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
