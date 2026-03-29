#' Example data for OptiSparseMET
#'
#' Synthetic dataset for demonstrating and testing sparse multi-environment
#' trial allocation and within-environment field design construction across the
#' functions of `OptiSparseMET`. All components are internally consistent:
#' treatment IDs, family labels, seed inventories, relationship matrices, and
#' environment specifications refer to the same set of synthetic lines and
#' environments, so any combination of components can be passed to the package
#' functions without modification.
#'
#' @format A named list with the following components:
#' \describe{
#'   \item{`treatments`}{Character vector of test treatment IDs. Candidate
#'     lines to be allocated across environments. Does not include check
#'     treatments.}
#'   \item{`environments`}{Character vector of environment names. Used as the
#'     `environments` argument in [allocate_sparse_met()] and as keys in
#'     `env_design_specs`.}
#'   \item{`treatment_info`}{Data frame with columns `Treatment` and `Family`.
#'     Used as the `treatment_info` argument in [allocate_sparse_met()] when
#'     `allocation_group_source = "Family"`.}
#'   \item{`common_treatments`}{Character vector of treatment IDs forced into
#'     all environments before sparse allocation. A subset of `treatments`.}
#'   \item{`seed_info`}{Data frame with columns `Treatment` and
#'     `SeedAvailable`. One row per treatment in `treatments`. Used as the
#'     `seed_available` argument in [assign_replication_by_seed()].}
#'   \item{`seed_required_per_plot`}{Data frame with columns `Environment` and
#'     `SeedRequiredPerPlot`. Per-plot seed requirement for each environment.}
#'   \item{`OptiSparseMET_GRM`}{Numeric matrix. Genomic relationship matrix
#'     with row and column names matching `treatments`. Used when
#'     `allocation_group_source = "GRM"`, `cluster_source = "GRM"`, or
#'     `dispersion_source = "GRM"`.}
#'   \item{`OptiSparseMET_A`}{Numeric matrix. Pedigree relationship matrix
#'     with row and column names matching `treatments`. Used when
#'     `allocation_group_source = "A"`, `cluster_source = "A"`, or
#'     `dispersion_source = "A"`.}
#'   \item{`OptiSparseMET_K`}{Numeric matrix. Relationship matrix for
#'     mixed-model efficiency evaluation and dispersion optimisation. Used
#'     when `prediction_type %in% c("GBLUP", "PBLUP")` or
#'     `dispersion_source = "K"`.}
#'   \item{`sparse_example_args_random_balanced`}{Named list of arguments
#'     ready to pass to [allocate_sparse_met()] with
#'     `allocation_method = "random_balanced"`.}
#'   \item{`sparse_example_args_balanced_incomplete`}{Named list of arguments
#'     ready to pass to [allocate_sparse_met()] with
#'     `allocation_method = "balanced_incomplete"`.}
#'   \item{`env_design_specs`}{Named list with one element per environment.
#'     Each element is a named list specifying the local field design for that
#'     environment. The `design` field in each element is set to either
#'     `"met_prep_famoptg"` or `"met_alpha_rc_stream"` -- the
#'     OptiSparseMET-specific versions of the within-environment constructors.
#'     These carry the `met_` prefix to avoid namespace conflicts when both
#'     OptiSparseMET and OptiDesign are loaded simultaneously.}
#' }
#'
#' @details
#' All data are synthetic and generated solely for use in package examples,
#' unit tests, and vignettes. No real breeding trial data are included.
#' Relationship matrices (`OptiSparseMET_GRM`, `OptiSparseMET_A`,
#' `OptiSparseMET_K`) are positive semi-definite by construction and have row
#' and column names consistent with `treatments`, so they can be used directly
#' without an `id_map`.
#'
#' The `seed_info` and `seed_required_per_plot` components are calibrated so
#' that a range of replication outcomes (fully replicated, partially replicated,
#' and excluded treatments) arise naturally when passed to
#' [assign_replication_by_seed()], making the data useful for demonstrating
#' all three replication modes.
#'
#' The `env_design_specs` component is ready to pass directly to
#' [plan_sparse_met_design()] via the `env_design_specs` argument. Each
#' specification uses `design = "met_prep_famoptg"` or
#' `design = "met_alpha_rc_stream"` to select the local design engine.
#'
#' ## Illustrative workflow
#'
#' ```r
#' data("OptiSparseMET_example_data", package = "OptiSparseMET")
#' x <- OptiSparseMET_example_data
#'
#' # Step 0: check feasibility
#' suggest_safe_k(x$treatments, x$environments, buffer = 3)
#'
#' # Step 1: allocate across environments
#' alloc <- do.call(allocate_sparse_met,
#'   x$sparse_example_args_random_balanced)
#'
#' # Step 2: seed-aware replication for one environment
#' rep_plan <- assign_replication_by_seed(
#'   treatments             = alloc$allocation_long$Treatment[
#'     alloc$allocation_long$Environment == "E1" &
#'     alloc$allocation_long$Assigned == 1L],
#'   seed_available         = x$seed_info,
#'   seed_required_per_plot = x$seed_required_per_plot$SeedRequiredPerPlot[
#'     x$seed_required_per_plot$Environment == "E1"],
#'   replication_mode       = "p_rep"
#' )
#'
#' # Step 3: within-environment design (using met_ prefix function)
#' design_E1 <- met_prep_famoptg(
#'   check_treatments        = x$common_treatments,
#'   check_families          = rep("CHECK", length(x$common_treatments)),
#'   p_rep_treatments        = rep_plan$p_rep_treatments,
#'   p_rep_reps              = rep_plan$p_rep_reps,
#'   p_rep_families          = x$treatment_info$Family[
#'     match(rep_plan$p_rep_treatments, x$treatment_info$Treatment)],
#'   unreplicated_treatments = rep_plan$unreplicated_treatments,
#'   unreplicated_families   = x$treatment_info$Family[
#'     match(rep_plan$unreplicated_treatments, x$treatment_info$Treatment)],
#'   n_blocks = 4L, n_rows = 10L, n_cols = 12L
#' )
#'
#' # Or run the full pipeline end-to-end
#' met_result <- plan_sparse_met_design(
#'   treatments                     = x$treatments,
#'   environments                   = x$environments,
#'   n_test_entries_per_environment = suggest_safe_k(x$treatments, x$environments),
#'   common_treatments              = x$common_treatments,
#'   env_design_specs               = x$env_design_specs,
#'   treatment_info                 = x$treatment_info,
#'   seed_info                      = x$seed_info,
#'   seed_required_per_plot         = x$seed_required_per_plot,
#'   seed                           = 123
#' )
#' ```
#'
#' @source Synthetically generated for package documentation and testing.
"OptiSparseMET_example_data"
