#' Plan a sparse multi-environment trial design and assemble a combined field book
#'
#' `plan_sparse_met_design()` is the end-to-end workflow wrapper for sparse
#' multi-environment trial (MET) design. It first allocates treatments across
#' environments using [allocate_sparse_met()], then builds a local field design
#' for each environment using the environment-specific specification in
#' `env_design_specs`, and finally stacks all local field books into one
#' combined MET field book using [combine_met_fieldbooks()].
#'
#' Compared with calling the lower-level functions manually,
#' `plan_sparse_met_design()` ensures that:
#' - every environment has a design specification
#' - allocation metadata are propagated into the combined field book
#' - local design metadata are harmonized across heterogeneous design engines
#' - seed-aware replication can be applied environment by environment
#' - local design efficiency outputs, when present as `design$efficiency`,
#'   are automatically captured into `environment_summary` and
#'   `efficiency_summary`
#' - environment specs are validated before design construction
#' - `met_prep_famoptg()` replication mode can be inferred automatically when
#'   not supplied explicitly
#' - a registry-based compatibility layer allows additional future design
#'   engines to be added with minimal changes
#'
#' @param treatments Character vector of treatment IDs.
#' @param environments Character vector of environment names.
#' @param allocation_method Character scalar passed to [allocate_sparse_met()].
#' @param n_test_entries_per_environment Integer scalar or integer vector passed
#'   to [allocate_sparse_met()].
#' @param target_replications Optional positive integer passed to
#'   [allocate_sparse_met()].
#' @param common_treatments Optional character vector of preselected common
#'   treatments. Supplied identities are fixed as present in every environment;
#'   their site-specific replication may differ according to the local design.
#' @param balance,balance_iter,balance_seed Connectivity refinement settings
#'   passed to [allocate_sparse_met()].
#' @param Sigma_E Optional environmental covariance matrix or candidate list
#'   used to allocate distinct shared treatments adaptively across environment
#'   pairs. Replication within a site does not contribute to this connectivity.
#' @param pair_target_se,pair_aggregate,pair_cvar_alpha Correlation-adaptive
#'   pairwise connectivity settings passed to [allocate_sparse_met()].
#' @param env_design_specs Named list with one element per environment. Each
#'   element must contain at least a `design` field selecting the local design
#'   engine (`"met_prep_famoptg"` or `"met_alpha_rc_stream"`). All other
#'   fields are passed through to the corresponding local design function,
#'   except a small set of pipeline-consumed fields used internally for
#'   seed-aware replication. A specification for `met_prep_famoptg` may include
#'   `fixed_treatment_reps`, a named positive-integer vector for a subset of the
#'   allocated treatments. Those exact replication counts are enforced first;
#'   remaining p-rep positions are then assigned by the usual seed-aware rule.
#'   Set `max_extra_replication_plots` when the field has an exact physical
#'   limit: a treatment fixed at three reps consumes two of those extra plots.
#' @param treatment_info Optional data frame with at least `Treatment` and
#'   optionally `Family`.
#' @param seed_info Optional data frame used by [assign_replication_by_seed()]
#'   when seed-aware replication is requested for a local design.
#' @param seed_required_per_plot Optional scalar, named vector, named list, or
#'   data frame with columns `Environment` and `SeedRequiredPerPlot`, used for
#'   network-wide seed accounting and local replication.
#' @param minimum_seed_buffer Non-negative reserve retained for every candidate
#'   treatment after all environments and local replications are planned.
#' @param allocation_group_source Character scalar passed to
#'   [allocate_sparse_met()].
#' @param GRM Optional genomic relationship matrix passed to
#'   [allocate_sparse_met()].
#' @param A Optional pedigree relationship matrix passed to
#'   [allocate_sparse_met()].
#' @param id_map Optional treatment-to-matrix ID map passed to
#'   [allocate_sparse_met()].
#' @param group_method Character scalar passed to [allocate_sparse_met()].
#' @param group_seed Integer seed passed to [allocate_sparse_met()].
#' @param group_attempts Integer passed to [allocate_sparse_met()].
#' @param n_pcs_use Integer or `Inf` passed to [allocate_sparse_met()].
#' @param min_groups_per_environment Optional integer passed to
#'   [allocate_sparse_met()].
#' @param min_env_per_group Optional integer passed to
#'   [allocate_sparse_met()].
#' @param balance_groups_across_env Logical passed to [allocate_sparse_met()].
#' @param force_group_connectivity Logical passed to [allocate_sparse_met()].
#' @param allow_approximate Logical passed to [allocate_sparse_met()].
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A named list with:
#' \describe{
#'   \item{`sparse_allocation`}{Output from [allocate_sparse_met()].}
#'   \item{`environment_designs`}{Named list of local design outputs, one per
#'     environment.}
#'   \item{`combined_field_book`}{Combined MET field book from
#'     [combine_met_fieldbooks()].}
#'   \item{`environment_summary`}{Environment-level summary table including
#'     design metadata, occupied plot counts (`n_total_plots`), grid cells
#'     (`n_field_cells`), treatment counts, and extracted efficiency fields such
#'     as `eff_model`, `eff_A`, `eff_D`, and `eff_mean_PEV` when available.}
#'   \item{`group_environment_summary`}{Group-by-environment summary propagated
#'     from `sparse_allocation$group_by_environment`.}
#'   \item{`efficiency_summary`}{Long-format table of extracted efficiency
#'     outputs with columns `Environment`, `LocalDesign`, `Metric`, `Value`,
#'     and `ValueType`. Empty when no efficiency outputs are detected.}
#'   \item{`seed_ledger`}{When `seed_info` is supplied, one row per treatment
#'     with the original inventory, grams allocated across all final fieldbooks,
#'     remaining inventory, and feasibility flag.}
#'   \item{`summary`}{High-level summary list for the overall MET design.}
#'   \item{`seed_used`}{Seed used internally.}
#' }
#'
#' @export
plan_sparse_met_design <- function(
    treatments,
    environments,
    allocation_method = c("random_balanced", "equireplicate", "M3", "M4"),
    n_test_entries_per_environment,
    target_replications = NULL,
    common_treatments = NULL,
    balance = c("none", "env_pair", "line_pair", "both"),
    balance_iter = 2000L,
    balance_seed = NULL,
    Sigma_E = NULL,
    pair_target_se = 0.15,
    pair_aggregate = c("maximin", "cvar", "mean"),
    pair_cvar_alpha = 0.25,
    env_design_specs,
    treatment_info = NULL,
    seed_info = NULL,
    seed_required_per_plot = NULL,
    minimum_seed_buffer = 0,
    allocation_group_source = c("none", "Family", "GRM", "A"),
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

  allocation_method <- match.arg(allocation_method)
  balance <- match.arg(balance)
  pair_aggregate <- match.arg(pair_aggregate)
  allocation_group_source <- match.arg(allocation_group_source)
  group_method <- match.arg(group_method)

  seed_used <- seed
  if (!is.null(seed_used)) {
    if (!is.numeric(seed_used) || length(seed_used) != 1L ||
        !is.finite(seed_used) || abs(seed_used - round(seed_used)) > 1e-8)
      stop("`seed` must be one finite integer or NULL.")
    seed_used <- as.integer(round(seed_used))
    set.seed(seed_used)
  }

  treatments   <- as.character(treatments)
  environments <- as.character(environments)
  if (anyNA(treatments) || any(!nzchar(treatments)) ||
      anyDuplicated(treatments))
    stop("`treatments` must contain unique, non-missing, non-empty IDs.")
  if (anyNA(environments) || any(!nzchar(environments)) ||
      anyDuplicated(environments))
    stop("`environments` must contain unique, non-missing, non-empty names.")
  if (length(minimum_seed_buffer) != 1L ||
      !is.numeric(minimum_seed_buffer) ||
      !is.finite(minimum_seed_buffer) || minimum_seed_buffer < 0)
    stop("`minimum_seed_buffer` must be one finite non-negative number.")
  if (is.null(seed_info) && minimum_seed_buffer > 0)
    stop("`seed_info` is required when `minimum_seed_buffer` is positive.")

  if (length(treatments) < 1L) {
    stop("`treatments` must contain at least one treatment ID.")
  }
  if (length(environments) < 1L) {
    stop("`environments` must contain at least one environment.")
  }

  if (missing(env_design_specs) || is.null(env_design_specs) || !is.list(env_design_specs)) {
    stop("`env_design_specs` must be a named list with one specification per environment.")
  }
  if (is.null(names(env_design_specs)) || any(!nzchar(names(env_design_specs)))) {
    stop("`env_design_specs` must be a named list, with names matching `environments`.")
  }

  missing_specs <- setdiff(environments, names(env_design_specs))
  if (length(missing_specs) > 0L) {
    stop(
      "Missing environment design specifications for: ",
      paste(missing_specs, collapse = ", ")
    )
  }

  `%||%` <- function(x, y) if (is.null(x)) y else x

  get_family_vec <- function(trt) {
    trt <- as.character(trt)
    if (!is.null(treatment_info) &&
        is.data.frame(treatment_info) &&
        all(c("Treatment", "Family") %in% names(treatment_info))) {
      fam <- treatment_info$Family[match(trt, treatment_info$Treatment)]
      fam[is.na(fam)] <- "UNGROUPED"
      return(as.character(fam))
    }
    rep("UNGROUPED", length(trt))
  }

  resolve_seed_required <- function(seed_required_per_plot, env_name) {
    default_val <- 1

    if (is.null(seed_required_per_plot)) return(default_val)

    if (is.numeric(seed_required_per_plot) && length(seed_required_per_plot) == 1L) {
      if (!is.finite(seed_required_per_plot) || seed_required_per_plot <= 0)
        stop("`seed_required_per_plot` must be strictly positive.")
      return(seed_required_per_plot)
    }

    if (is.numeric(seed_required_per_plot) &&
        length(seed_required_per_plot) > 1L &&
        !is.null(names(seed_required_per_plot))) {
      val <- seed_required_per_plot[[env_name]]
      if (is.null(val) || length(val) != 1L ||
          !is.finite(val) || val <= 0)
        stop("Named `seed_required_per_plot` is missing or invalid for `",
             env_name, "`.")
      return(val)
    }

    if (is.data.frame(seed_required_per_plot)) {
      if (!all(c("Environment", "SeedRequiredPerPlot") %in%
               names(seed_required_per_plot)))
        stop("`seed_required_per_plot` data frame needs columns ",
             "`Environment` and `SeedRequiredPerPlot`.")
      if (anyDuplicated(as.character(seed_required_per_plot$Environment)))
        stop("`seed_required_per_plot` contains duplicate environments.")
      idx <- match(env_name, as.character(seed_required_per_plot$Environment))
      if (is.na(idx))
        stop("`seed_required_per_plot` is missing environment `", env_name, "`.")
      val <- seed_required_per_plot$SeedRequiredPerPlot[idx]
      if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0)
        stop("Invalid SeedRequiredPerPlot for environment `", env_name, "`.")
      return(val)
    }

    if (is.list(seed_required_per_plot)) {
      val <- seed_required_per_plot[[env_name]]
      if (is.null(val) || length(val) != 1L || !is.numeric(val) ||
          !is.finite(val) || val <= 0)
        stop("List `seed_required_per_plot` is missing or invalid for `",
             env_name, "`.")
      return(val)
    }

    stop("`seed_required_per_plot` must be NULL, a positive scalar, a named ",
         "vector/list, or an Environment/SeedRequiredPerPlot data frame.")
  }

  extract_efficiency_wide <- function(local_out) {
    eff <- local_out$efficiency
    if (is.null(eff) || !is.list(eff)) return(NULL)

    out <- list(
      eff_model = if (!is.null(eff$model)) as.character(eff$model) else NA_character_,
      eff_treatment_effect = if (!is.null(eff$treatment_effect)) as.character(eff$treatment_effect) else NA_character_,
      eff_prediction_type = if (!is.null(eff$prediction_type)) as.character(eff$prediction_type) else NA_character_,
      eff_residual_structure_requested = if (!is.null(eff$residual_structure_requested)) as.character(eff$residual_structure_requested) else NA_character_,
      eff_residual_structure_used = if (!is.null(eff$residual_structure_used)) as.character(eff$residual_structure_used) else NA_character_,
      eff_spatial_engine_used = if (!is.null(eff$spatial_engine_used)) as.character(eff$spatial_engine_used) else NA_character_,
      eff_mode = if (!is.null(eff$mode)) as.character(eff$mode) else NA_character_,
      eff_notes = if (!is.null(eff$notes)) as.character(eff$notes) else NA_character_,
      eff_A = if (!is.null(eff$A)) as.numeric(eff$A) else NA_real_,
      eff_D = if (!is.null(eff$D)) as.numeric(eff$D) else NA_real_,
      eff_mean_PEV = if (!is.null(eff$mean_PEV)) as.numeric(eff$mean_PEV) else NA_real_,
      eff_n_lines = if (!is.null(eff$n_lines)) as.numeric(eff$n_lines) else NA_real_,
      has_efficiency = TRUE,
      n_efficiency_metrics = sum(!vapply(eff, is.null, logical(1)))
    )

    as.data.frame(out, stringsAsFactors = FALSE)
  }

  extract_efficiency_long <- function(local_out, env_name, design_name) {
    eff <- local_out$efficiency

    empty <- data.frame(
      Environment = character(0), LocalDesign = character(0),
      Metric = character(0), Value = character(0), ValueType = character(0),
      stringsAsFactors = FALSE
    )

    if (is.null(eff) || !is.list(eff)) return(empty)

    vals <- lapply(names(eff), function(nm) {
      val <- eff[[nm]]
      if (length(val) == 1L) {
        data.frame(
          Environment = env_name, LocalDesign = design_name,
          Metric = nm, Value = as.character(val),
          ValueType = if (is.numeric(val)) "numeric" else "character",
          stringsAsFactors = FALSE
        )
      } else NULL
    })

    vals <- vals[!vapply(vals, is.null, logical(1))]
    if (length(vals) == 0L) return(empty)
    do.call(rbind, vals)
  }

  build_environment_summary_row <- function(
    env_name, design_name, spec, local_out, env_treatments, common_treatments
  ) {
    fb <- local_out$field_book
    n_field_cells <- if (!is.null(fb) && is.data.frame(fb)) nrow(fb) else 0L
    n_total_plots <- if (!is.null(fb) && is.data.frame(fb) &&
                         "Treatment" %in% names(fb)) {
      trt <- as.character(fb$Treatment)
      sum(!is.na(trt) & nzchar(trt))
    } else n_field_cells
    n_unique_treatments <- if (!is.null(fb) && is.data.frame(fb) && "Treatment" %in% names(fb)) {
      trt <- as.character(fb$Treatment)
      length(unique(trt[!is.na(trt) & nzchar(trt)]))
    } else length(unique(env_treatments))

    n_allocated_treatments <- length(env_treatments)
    n_common_present <- sum(env_treatments %in% common_treatments)
    n_sparse_present <- n_allocated_treatments - n_common_present

    data.frame(
      Environment = env_name,
      LocalDesign = design_name,
      ReplicationMode = if (!is.null(spec$replication_mode)) as.character(spec$replication_mode) else NA_character_,
      n_total_plots = as.integer(n_total_plots),
      n_field_cells = as.integer(n_field_cells),
      n_allocated_treatments = as.integer(n_allocated_treatments),
      n_unique_treatments = as.integer(n_unique_treatments),
      n_common_treatments = as.integer(n_common_present),
      n_sparse_treatments = as.integer(n_sparse_present),
      stringsAsFactors = FALSE
    )
  }

  validate_env_spec <- function(spec, env_name) {
    if (is.null(spec) || !is.list(spec))
      stop("Each element of `env_design_specs` must be a list. Problem at environment: ", env_name)
    if (is.null(spec$design))
      stop("Each environment specification must include a `design` field. Missing for environment: ", env_name)

    design_name <- as.character(spec$design)[1L]

    # Validate required fields for built-in engines
    if (identical(design_name, "met_prep_famoptg")) {
      req <- c("check_treatments", "check_families", "n_blocks", "n_rows", "n_cols")
      miss <- req[!req %in% names(spec)]
      if (length(miss) > 0L)
        stop("Environment `", env_name, "` using `met_prep_famoptg` is missing required fields: ",
             paste(miss, collapse = ", "))
      if (length(spec$check_treatments) != length(spec$check_families))
        stop("Environment `", env_name, "`: `check_treatments` and `check_families` must have the same length.")

    } else if (identical(design_name, "met_alpha_rc_stream")) {
      req <- c("check_treatments", "check_families", "n_reps", "n_rows", "n_cols")
      miss <- req[!req %in% names(spec)]
      if (length(miss) > 0L)
        stop("Environment `", env_name, "` using `met_alpha_rc_stream` is missing required fields: ",
             paste(miss, collapse = ", "))
      if (length(spec$check_treatments) != length(spec$check_families))
        stop("Environment `", env_name, "`: `check_treatments` and `check_families` must have the same length.")
    }
    invisible(TRUE)
  }

  infer_replication_mode <- function(spec, env_name, env_treatments) {
    if (!is.null(spec$replication_mode)) return(as.character(spec$replication_mode)[1L])
    desired_replications_used <- spec$desired_replications %||% 2L
    env_seed_req <- resolve_seed_required(seed_required_per_plot, env_name)
    if (!is.null(seed_info) && !is.null(env_seed_req) && env_seed_req > 0 && length(env_treatments) > 0L) {
      return(if (desired_replications_used > 1L) "p_rep" else "augmented")
    }
    "augmented"
  }

  build_prep_famoptg_args <- function(spec, env_name, env_treatments) {
    local_args <- spec
    replication_mode_used <- infer_replication_mode(spec, env_name, env_treatments)
    desired_replications_used <- local_args$desired_replications %||% 2L
    candidate_prep_used <- local_args$candidate_prep %||% NULL
    max_prep_used <- local_args$max_prep %||% NULL
    max_extra_plots_used <- local_args$max_extra_replication_plots %||% NULL
    shortage_action_used <- local_args$shortage_action %||% "downgrade"
    fixed_reps_used <- local_args$fixed_treatment_reps %||% integer(0)

    # Remove pipeline-only arguments
    local_args$design <- NULL
    local_args$replication_mode <- NULL
    local_args$desired_replications <- NULL
    local_args$candidate_prep <- NULL
    local_args$max_prep <- NULL
    local_args$max_extra_replication_plots <- NULL
    local_args$shortage_action <- NULL
    local_args$fixed_treatment_reps <- NULL

    if (length(fixed_reps_used)) {
      if (!is.numeric(fixed_reps_used) || is.null(names(fixed_reps_used)) ||
          any(!nzchar(names(fixed_reps_used))) ||
          anyDuplicated(names(fixed_reps_used)) ||
          any(!is.finite(fixed_reps_used)) ||
          any(abs(fixed_reps_used - round(fixed_reps_used)) > 1e-8) ||
          any(fixed_reps_used < 1L))
        stop("Environment `", env_name, "`: `fixed_treatment_reps` must be ",
             "a uniquely named vector of positive integers.")
      fixed_reps_used <- as.integer(round(fixed_reps_used))
      names(fixed_reps_used) <- names(spec$fixed_treatment_reps)
      missing_fixed <- setdiff(names(fixed_reps_used), env_treatments)
      if (length(missing_fixed))
        stop("Environment `", env_name, "`: fixed-replication treatments ",
             "were not allocated here: ", paste(missing_fixed, collapse = ", "))
      if (!is.null(local_args$n_blocks) &&
          any(fixed_reps_used > as.integer(local_args$n_blocks)))
        stop("Environment `", env_name, "`: a fixed replication exceeds ",
             "`n_blocks`.")
      if (!identical(replication_mode_used, "p_rep"))
        stop("Environment `", env_name, "`: `fixed_treatment_reps` is ",
             "supported only with `replication_mode = \"p_rep\"`.")
    }
    if (!is.null(max_extra_plots_used)) {
      if (!is.numeric(max_extra_plots_used) ||
          length(max_extra_plots_used) != 1L ||
          !is.finite(max_extra_plots_used) ||
          abs(max_extra_plots_used - round(max_extra_plots_used)) > 1e-8 ||
          max_extra_plots_used < 0)
        stop("Environment `", env_name, "`: ",
             "`max_extra_replication_plots` must be one non-negative integer.")
      max_extra_plots_used <- as.integer(round(max_extra_plots_used))
    }

    seed_req <- resolve_seed_required(seed_required_per_plot, env_name)
    if (!is.numeric(seed_req) || length(seed_req) != 1L || is.na(seed_req) || seed_req <= 0) {
      stop("Invalid seed_required_per_plot for environment ", env_name,
           ". Must resolve to a single positive numeric value.")
    }
    fixed_extra_plots_n <- sum(pmax(0L, fixed_reps_used - 1L))
    if (!is.null(max_extra_plots_used) &&
        fixed_extra_plots_n > max_extra_plots_used)
      stop("Environment `", env_name, "`: fixed replication requires ",
           fixed_extra_plots_n, " extra plots, exceeding ",
           "`max_extra_replication_plots = ", max_extra_plots_used, "`.")

    if (replication_mode_used %in% c("p_rep", "rcbd_type") &&
        !is.null(seed_info) && seed_req > 0) {

      fixed_ids <- names(fixed_reps_used)
      fixed_extra <- pmax(0, fixed_reps_used - 1L) * seed_req
      fixed_replicated_n <- sum(fixed_reps_used > 1L)
      if (!is.null(seed_remaining) && length(fixed_ids)) {
        if (any(fixed_extra >
                seed_remaining[fixed_ids] + 1e-8))
          stop("Environment `", env_name, "`: optimized fixed replication ",
               "exceeds the network-wide seed inventory for: ",
               paste(fixed_ids[
                 fixed_extra > seed_remaining[fixed_ids] + 1e-8
               ], collapse = ", "), ".")
        seed_remaining[fixed_ids] <<-
          seed_remaining[fixed_ids] - fixed_extra
      }

      remaining_treatments <- setdiff(env_treatments, fixed_ids)
      local_seed_info <- seed_info
      if (!is.null(seed_remaining)) {
        idx <- match(remaining_treatments, as.character(seed_info$Treatment))
        if (anyNA(idx))
          stop("Seed inventory does not cover every allocated treatment.")
        local_seed_info <- data.frame(
          Treatment = remaining_treatments,
          # One mandatory base plot was already debited during allocation.
          SeedAvailable = seed_req + seed_remaining[remaining_treatments],
          stringsAsFactors = FALSE)
      }
      max_prep_remaining <- if (is.null(max_prep_used)) NULL else
        max(0L, as.integer(max_prep_used) - fixed_replicated_n)
      if (!is.null(max_extra_plots_used)) {
        auto_extra_per_treatment <- max(1L,
          as.integer(desired_replications_used) - 1L)
        by_plot_budget <- floor(
          (max_extra_plots_used - fixed_extra_plots_n) /
            auto_extra_per_treatment
        )
        max_prep_remaining <- if (is.null(max_prep_remaining))
          by_plot_budget else min(max_prep_remaining, by_plot_budget)
      }
      candidate_remaining <- if (is.null(candidate_prep_used))
        remaining_treatments else
        intersect(as.character(candidate_prep_used), remaining_treatments)

      if (length(remaining_treatments)) {
        rep_out_auto <- assign_replication_by_seed(
          treatments = remaining_treatments,
          seed_available = local_seed_info,
          seed_required_per_plot = seed_req,
          replication_mode = replication_mode_used,
          desired_replications = desired_replications_used,
          candidate_prep = candidate_remaining,
          max_prep = max_prep_remaining,
          shortage_action = shortage_action_used,
          check_treatments = local_args$check_treatments
        )
      } else {
        rep_out_auto <- list(
          p_rep_treatments = character(0),
          p_rep_reps = integer(0),
          unreplicated_treatments = character(0),
          excluded_treatments = character(0),
          seed_summary = data.frame(),
          replication_mode = replication_mode_used,
          seed_used = NULL
        )
      }

      forced_rep_ids <- fixed_ids[fixed_reps_used > 1L]
      p_rep_treatments <- c(forced_rep_ids,
                            rep_out_auto$p_rep_treatments)
      p_rep_reps <- c(unname(fixed_reps_used[forced_rep_ids]),
                      rep_out_auto$p_rep_reps)
      unrep_treatments <- c(fixed_ids[fixed_reps_used == 1L],
                            rep_out_auto$unreplicated_treatments)
      ord_rep <- match(p_rep_treatments, env_treatments)
      p_rep_reps <- p_rep_reps[order(ord_rep)]
      p_rep_treatments <- p_rep_treatments[order(ord_rep)]
      unrep_treatments <- env_treatments[
        env_treatments %in% unrep_treatments
      ]

      auto_rep <- rep_out_auto$p_rep_treatments
      if (!is.null(seed_remaining) && length(auto_rep)) {
        extra_cost <- pmax(0, as.numeric(rep_out_auto$p_rep_reps) - 1) *
          seed_req
        seed_remaining[auto_rep] <<-
          seed_remaining[auto_rep] - extra_cost
        if (any(seed_remaining < -1e-8))
          stop("Local replication planning exceeded the network-wide seed ",
               "inventory in environment `", env_name, "`.")
      }
      rep_out <- rep_out_auto
      rep_out$p_rep_treatments <- p_rep_treatments
      rep_out$p_rep_reps <- p_rep_reps
      rep_out$unreplicated_treatments <- unrep_treatments
      rep_out$fixed_treatment_reps <- fixed_reps_used
      attr(local_args, "replication_plan") <- rep_out

    } else if (replication_mode_used == "augmented") {
      p_rep_treatments <- character(0)
      p_rep_reps <- integer(0)
      unrep_treatments <- env_treatments

    } else if (replication_mode_used == "rcbd_type") {
      p_rep_treatments <- env_treatments
      p_rep_reps <- rep(as.integer(desired_replications_used), length(env_treatments))
      unrep_treatments <- character(0)

    } else {
      fixed_ids <- names(fixed_reps_used)
      if (is.null(fixed_ids)) fixed_ids <- character(0)
      p_rep_treatments <- fixed_ids[fixed_reps_used > 1L]
      p_rep_reps <- unname(fixed_reps_used[p_rep_treatments])
      unrep_treatments <- c(fixed_ids[fixed_reps_used == 1L],
                            setdiff(env_treatments, fixed_ids))
    }

    # Assign typed zero-length vectors, never NULL: `$<- NULL` removes a list
    # element and would make the required formal disappear from `do.call()`.
    p_rep_treatments <- p_rep_treatments %||% character(0)
    p_rep_reps <- p_rep_reps %||% integer(0)
    unrep_treatments <- unrep_treatments %||% character(0)
    local_args$p_rep_treatments <- as.character(p_rep_treatments)
    local_args$p_rep_reps <- as.integer(p_rep_reps)
    local_args$p_rep_families <- get_family_vec(p_rep_treatments)
    local_args$unreplicated_treatments <- as.character(unrep_treatments)
    local_args$unreplicated_families <- get_family_vec(unrep_treatments)
    attr(local_args, "replication_mode_used") <- replication_mode_used
    local_args
  }

  build_alpha_rc_args <- function(spec, env_treatments) {
    local_args <- spec
    local_args$design <- NULL
    local_args$entry_treatments <- env_treatments
    local_args$entry_families <- get_family_vec(env_treatments)
    local_args
  }

  build_generic_engine_args <- function(fun, spec, env_treatments) {
    local_args <- spec
    local_args$design <- NULL
    fn_formals <- names(formals(fun))
    trt_aliases <- c("entry_treatments", "treatments", "entries", "genotypes", "lines")
    fam_aliases <- c("entry_families", "families", "entry_family", "family_vec", "line_families")
    trt_name <- trt_aliases[trt_aliases %in% fn_formals]
    fam_name <- fam_aliases[fam_aliases %in% fn_formals]
    if (length(trt_name) > 0L) local_args[[trt_name[1L]]] <- env_treatments
    if (length(fam_name) > 0L) local_args[[fam_name[1L]]] <- get_family_vec(env_treatments)
    local_args
  }

  resolve_design_engine <- function(design_spec) {
    if (is.function(design_spec))
      return(list(name = deparse(substitute(design_spec)), fun = design_spec, type = "custom"))

    design_name <- as.character(design_spec)[1L]

    # Registry of built-in engines — met_ prefix versions only
    builtins <- list(
      met_prep_famoptg    = met_prep_famoptg,
      met_alpha_rc_stream = met_alpha_rc_stream
    )

    if (design_name %in% names(builtins))
      return(list(name = design_name, fun = builtins[[design_name]], type = "builtin"))

    if (exists(design_name, mode = "function", inherits = TRUE))
      return(list(name = design_name, fun = get(design_name, mode = "function"), type = "custom"))

    stop("Unsupported local design engine `", design_name,
         "`. Use `\"met_prep_famoptg\"` or `\"met_alpha_rc_stream\"`.")
  }

  # ============================================================
  # 1. Across-environment allocation
  # ============================================================
  allocation_seed_cost <- NULL
  if (!is.null(seed_info)) {
    if (!is.data.frame(seed_info) ||
        !all(c("Treatment", "SeedAvailable") %in% names(seed_info)))
      stop("`seed_info` must contain `Treatment` and `SeedAvailable` columns.")
    if (anyDuplicated(as.character(seed_info$Treatment)))
      stop("`seed_info` contains duplicate treatment IDs.")
    allocation_seed_cost <- stats::setNames(vapply(
      environments,
      function(env_name) {
        spec <- env_design_specs[[env_name]]
        validate_env_spec(spec, env_name)
        q <- resolve_seed_required(seed_required_per_plot, env_name)
        mandatory_reps <- if (identical(
          as.character(spec$design)[1L], "met_alpha_rc_stream"))
          as.integer(spec$n_reps) else 1L
        q * mandatory_reps
      }, numeric(1)), environments)
  }

  sparse_out <- allocate_sparse_met(
    treatments = treatments,
    environments = environments,
    allocation_method = allocation_method,
    n_test_entries_per_environment = n_test_entries_per_environment,
    target_replications = target_replications,
    common_treatments = common_treatments,
    seed_available = if (!is.null(allocation_seed_cost)) seed_info else NULL,
    seed_required_per_environment = allocation_seed_cost,
    minimum_seed_buffer = minimum_seed_buffer,
    balance = balance,
    balance_iter = balance_iter,
    balance_seed = balance_seed,
    Sigma_E = Sigma_E,
    pair_target_se = pair_target_se,
    pair_aggregate = pair_aggregate,
    pair_cvar_alpha = pair_cvar_alpha,
    allocation_group_source = allocation_group_source,
    treatment_info = treatment_info,
    GRM = GRM,
    A = A,
    id_map = id_map,
    group_method = group_method,
    group_seed = group_seed,
    group_attempts = group_attempts,
    n_pcs_use = n_pcs_use,
    min_groups_per_environment = min_groups_per_environment,
    min_env_per_group = min_env_per_group,
    balance_groups_across_env = balance_groups_across_env,
    force_group_connectivity = force_group_connectivity,
    allow_approximate = allow_approximate,
    seed = seed
  )
  seed_remaining <- if (!is.null(sparse_out$seed_summary))
    stats::setNames(sparse_out$seed_summary$SeedSpendableRemaining,
                    sparse_out$seed_summary$Treatment) else NULL

  # ============================================================
  # 2. Per-environment local design construction
  # ============================================================
  environment_designs <- vector("list", length(environments))
  names(environment_designs) <- environments

  local_designs    <- stats::setNames(rep(NA_character_, length(environments)), environments)
  replication_modes <- stats::setNames(rep(NA_character_, length(environments)), environments)

  env_summary_rows <- list()
  efficiency_rows  <- list()

  for (env_name in environments) {
    spec <- env_design_specs[[env_name]]
    validate_env_spec(spec, env_name)

    engine      <- resolve_design_engine(spec$design)
    design_name <- engine$name
    design_fun  <- engine$fun

    local_designs[env_name] <- design_name

    env_treatments <- rownames(sparse_out$allocation_matrix)[
      sparse_out$allocation_matrix[, env_name] == 1L
    ]

    if (identical(design_name, "met_prep_famoptg")) {
      local_args <- build_prep_famoptg_args(spec, env_name, env_treatments)
      replication_modes[env_name] <- attr(local_args, "replication_mode_used")

    } else if (identical(design_name, "met_alpha_rc_stream")) {
      local_args <- build_alpha_rc_args(spec, env_treatments)
      replication_modes[env_name] <- spec$replication_mode %||% NA_character_

    } else {
      local_args <- build_generic_engine_args(design_fun, spec, env_treatments)
      replication_modes[env_name] <- spec$replication_mode %||% NA_character_
    }

    replication_plan <- attr(local_args, "replication_plan")
    local_out <- do.call(design_fun, local_args)
    if (!is.null(replication_plan))
      local_out$replication_plan <- replication_plan

    if (!is.list(local_out) || !"field_book" %in% names(local_out))
      stop("Local design function for environment ", env_name,
           " must return a list containing `field_book`.")
    if (!is.data.frame(local_out$field_book))
      stop("Local design function for environment ", env_name,
           " returned `field_book`, but it is not a data frame.")

    environment_designs[[env_name]] <- local_out

    env_row <- build_environment_summary_row(
      env_name = env_name, design_name = design_name, spec = spec,
      local_out = local_out, env_treatments = env_treatments,
      common_treatments = common_treatments %||% character(0)
    )

    eff_wide <- extract_efficiency_wide(local_out)
    eff_long <- extract_efficiency_long(local_out, env_name, design_name)

    if (!is.null(eff_wide)) {
      env_row <- cbind(env_row, eff_wide)
    } else {
      env_row$eff_model                      <- NA_character_
      env_row$eff_treatment_effect           <- NA_character_
      env_row$eff_prediction_type            <- NA_character_
      env_row$eff_residual_structure_requested <- NA_character_
      env_row$eff_residual_structure_used    <- NA_character_
      env_row$eff_spatial_engine_used        <- NA_character_
      env_row$eff_mode                       <- NA_character_
      env_row$eff_notes                      <- NA_character_
      env_row$eff_A                          <- NA_real_
      env_row$eff_D                          <- NA_real_
      env_row$eff_mean_PEV                   <- NA_real_
      env_row$eff_n_lines                    <- NA_real_
      env_row$has_efficiency                 <- FALSE
      env_row$n_efficiency_metrics           <- 0L
    }

    env_summary_rows[[env_name]] <- env_row
    if (!is.null(eff_long) && nrow(eff_long) > 0L)
      efficiency_rows[[env_name]] <- eff_long
  }

  # ============================================================
  # 3. Combine local field books
  # ============================================================
  combined_field_book <- combine_met_fieldbooks(
    field_books       = lapply(environment_designs, `[[`, "field_book"),
    local_designs     = local_designs,
    replication_modes = replication_modes,
    sparse_method     = sparse_out$summary$allocation_method,
    common_treatments = common_treatments
  )

  # Reconcile the deliverable fieldbooks against one network-wide inventory.
  # Checks absent from seed_info are excluded because their seed is managed
  # separately from candidate-entry seed.
  seed_ledger <- NULL
  if (!is.null(seed_info)) {
    seed_ids <- as.character(seed_info$Treatment)
    used <- stats::setNames(numeric(length(seed_ids)), seed_ids)
    for (env_name in environments) {
      fb <- environment_designs[[env_name]]$field_book
      q <- resolve_seed_required(seed_required_per_plot, env_name)
      tt <- table(as.character(fb$Treatment))
      ids <- intersect(names(tt), seed_ids)
      used[ids] <- used[ids] + as.numeric(tt[ids]) * q
    }
    available <- stats::setNames(as.numeric(seed_info$SeedAvailable), seed_ids)
    remaining <- available - used
    seed_ledger <- data.frame(
      Treatment = seed_ids,
      SeedAvailable = available[seed_ids],
      MinimumBuffer = rep(minimum_seed_buffer, length(seed_ids)),
      SeedAllocated = used[seed_ids],
      SeedRemaining = remaining[seed_ids],
      Feasible = remaining[seed_ids] >= minimum_seed_buffer - 1e-8,
      stringsAsFactors = FALSE)
    if (any(!seed_ledger$Feasible)) {
      bad <- seed_ledger$Treatment[!seed_ledger$Feasible]
      stop("Final fieldbooks exceed network-wide seed availability for: ",
           paste(bad, collapse = ", "), ".")
    }
  }

  # ============================================================
  # 4. Build summary tables
  # ============================================================
  environment_summary <- do.call(rbind, env_summary_rows)
  rownames(environment_summary) <- NULL

  if (length(efficiency_rows) > 0L) {
    efficiency_summary <- do.call(rbind, efficiency_rows)
    rownames(efficiency_summary) <- NULL
  } else {
    efficiency_summary <- data.frame(
      Environment = character(0), LocalDesign = character(0),
      Metric = character(0), Value = character(0), ValueType = character(0),
      stringsAsFactors = FALSE
    )
  }

  environment_summary <- environment_summary[
    match(environments, environment_summary$Environment), , drop = FALSE
  ]
  rownames(environment_summary) <- NULL

  # ============================================================
  # 5. Overall summary
  # ============================================================
  summary_out <- list(
    n_environments = length(environments),
    allocation_method = sparse_out$summary$allocation_method,
    n_total_test_treatments = length(treatments),
    n_common_treatments = length(common_treatments %||% character(0)),
    n_sparse_treatments = sparse_out$summary$n_sparse_treatments,
    n_total_plots = sum(environment_summary$n_total_plots),
    n_environments_with_efficiency = sum(environment_summary$has_efficiency, na.rm = TRUE)
  )

  list(
    sparse_allocation       = sparse_out,
    environment_designs     = environment_designs,
    combined_field_book     = combined_field_book,
    environment_summary     = environment_summary,
    group_environment_summary = sparse_out$group_by_environment,
    efficiency_summary      = efficiency_summary,
    seed_ledger              = seed_ledger,
    summary                 = summary_out,
    seed_used               = seed_used
  )
}
