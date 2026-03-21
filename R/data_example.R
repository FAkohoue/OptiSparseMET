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
#'   \item{`treatments`}{Character vector of test treatment IDs. These are the
#'     candidate lines to be allocated across environments. Does not include
#'     check treatments.}
#'   \item{`environments`}{Character vector of environment names. Used as the
#'     `environments` argument in `allocate_sparse_met()` and as keys in
#'     `env_design_specs`.}
#'   \item{`treatment_info`}{Data frame with one row per treatment in
#'     `treatments`. Contains at minimum columns `Treatment` and `Family`.
#'     Used as the `treatment_info` argument in `allocate_sparse_met()` when
#'     `allocation_group_source = "Family"`.}
#'   \item{`common_treatments`}{Character vector of treatment IDs forced into
#'     all environments before sparse allocation. A subset of `treatments`.
#'     Used as the `common_treatments` argument in `allocate_sparse_met()` and
#'     `combine_met_fieldbooks()`.}
#'   \item{`seed_info`}{Data frame with columns `Treatment` and
#'     `SeedAvailable`. Contains one row per treatment in `treatments` and
#'     gives the seed quantity available for each line. Used as the
#'     `seed_available` argument in `assign_replication_by_seed()`.}
#'   \item{`seed_required_per_plot`}{Data frame with columns `Environment` and
#'     `SeedRequiredPerPlot`. Gives the per-plot seed requirement for each
#'     environment, reflecting that seed consumption may differ across
#'     locations. Used as the `seed_required_per_plot` argument in
#'     `assign_replication_by_seed()` after subsetting to the environment of
#'     interest.}
#'   \item{`OptiSparseMET_GRM`}{Numeric matrix. Genomic relationship matrix
#'     with row and column names matching the IDs in `treatments`. Used as the
#'     `GRM` argument in `allocate_sparse_met()` when
#'     `allocation_group_source = "GRM"`, and in `alpha_rc_stream()` or
#'     `prep_famoptg()` when `cluster_source = "GRM"` or
#'     `dispersion_source = "GRM"`.}
#'   \item{`OptiSparseMET_A`}{Numeric matrix. Pedigree-based numerator
#'     relationship matrix with row and column names matching the IDs in
#'     `treatments`. Used as the `A` argument when
#'     `allocation_group_source = "A"`, `cluster_source = "A"`, or
#'     `dispersion_source = "A"`.}
#'   \item{`OptiSparseMET_K`}{Numeric matrix. Relationship matrix intended for
#'     mixed-model efficiency evaluation and dispersion optimization. Used as
#'     the `K` argument in `alpha_rc_stream()` and `prep_famoptg()` when
#'     `prediction_type %in% c("GBLUP", "PBLUP")` or
#'     `dispersion_source = "K"`.}
#'   \item{`sparse_example_args_random_balanced`}{Named list of arguments ready
#'     to pass to `allocate_sparse_met()` with
#'     `allocation_method = "random_balanced"`. All list names correspond
#'     directly to function parameter names.}
#'   \item{`sparse_example_args_balanced_incomplete`}{Named list of arguments
#'     ready to pass to `allocate_sparse_met()` with
#'     `allocation_method = "balanced_incomplete"`. All list names correspond
#'     directly to function parameter names.}
#'   \item{`env_design_specs`}{Named list with one element per environment in
#'     `environments`. Each element is itself a named list of arguments
#'     specifying the local field design for that environment, suitable for
#'     passing to `prep_famoptg()` or `alpha_rc_stream()`.}
#' }
#'
#' @details
#' All data are synthetic and generated solely for use in package examples,
#' unit tests, and vignettes. No real breeding trial data are included.
#' Relationship matrices (`OptiSparseMET_GRM`, `OptiSparseMET_A`,
#' `OptiSparseMET_K`) are positive semi-definite by construction and have
#' row and column names consistent with `treatments`, so they can be used
#' directly without an `id_map`. The `seed_info` and `seed_required_per_plot`
#' components are calibrated so that a range of replication outcomes (fully
#' replicated, partially replicated, and excluded treatments) arise naturally
#' when passed to `assign_replication_by_seed()`, making the data useful for
#' demonstrating all three replication modes.
#'
#' @source Synthetically generated for package documentation and testing.
"OptiSparseMET_example_data"