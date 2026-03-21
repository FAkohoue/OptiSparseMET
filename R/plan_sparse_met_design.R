standardize_local_field_book <- function(fb) {
  if (!is.data.frame(fb)) {
    stop("Local field book must be a data.frame.")
  }
  
  core_cols <- c("Treatment", "Family", "Gcluster", "Block", "Plot", "Row", "Column")
  
  for (cc in core_cols) {
    if (!cc %in% names(fb)) {
      fb[[cc]] <- NA
    }
  }
  
  extra_cols <- setdiff(names(fb), core_cols)
  fb <- fb[, c(core_cols, extra_cols), drop = FALSE]
  
  fb$Treatment <- as.character(fb$Treatment)
  fb$Family <- as.character(fb$Family)
  fb$Gcluster <- as.character(fb$Gcluster)
  
  fb
}

#' Run the full two-stage sparse MET design pipeline
#'
#' `plan_sparse_met_design()` executes the complete OptiSparseMET pipeline in a
#' single call: across-environment treatment allocation via
#' `allocate_sparse_met()`, followed by within-environment field layout
#' construction via `prep_famoptg()` or `alpha_rc_stream()` for each
#' environment, followed by assembly of the combined MET field book via
#' `combine_met_fieldbooks()`. The function is the recommended entry point when
#' the full pipeline is to be run without manual intervention between stages.
#'
#' @description
#' Each environment is designed independently after allocation, using whichever
#' local design engine is specified in `env_design_specs`. Environments may use
#' different engines and different design parameters. When a `prep_famoptg()`
#' environment operates under `"p_rep"` or `"rcbd_type"` replication mode and
#' both `seed_info` and `seed_required_per_plot` are supplied, replication
#' roles are determined by `assign_replication_by_seed()` before the field
#' layout is constructed. When seed information is absent, replication roles
#' are derived directly from the `replication_mode` and `desired_replications`
#' fields in `env_design_specs`.
#'
#' Genetic structure can guide both the across-environment allocation (via
#' `allocation_group_source`) and the within-environment entry arrangement (via
#' the `cluster_source` field in each environment's specification). These two
#' layers are controlled independently; a GRM-guided allocation does not
#' automatically imply GRM-guided within-environment arrangement, and vice versa.
#'
#' @details
#' ## `env_design_specs` structure
#'
#' `env_design_specs` is a named list with one element per environment. Each
#' element is itself a named list of arguments passed to the local design
#' function. Two fields control the pipeline behaviour for `prep_famoptg()`
#' environments and are consumed by `plan_sparse_met_design()` before the
#' remaining fields are forwarded:
#'
#' - `design`: character scalar, either `"prep_famoptg"` or
#'   `"alpha_rc_stream"`. Required for every environment.
#' - `replication_mode`: character scalar passed to `assign_replication_by_seed()`
#'   when seed information is available; one of `"augmented"`, `"p_rep"`, or
#'   `"rcbd_type"`. Required for `prep_famoptg()` environments; ignored for
#'   `alpha_rc_stream()` environments.
#'
#' Additional fields consumed by `plan_sparse_met_design()` for `prep_famoptg()`
#' environments: `desired_replications`, `candidate_prep`, `max_prep`,
#' `shortage_action`. All remaining fields are forwarded as-is to the local
#' design function, including `check_treatments`, `check_families`, `n_blocks`,
#' `n_rows`, `n_cols`, `cluster_source`, and any other argument accepted by the
#' respective function.
#'
#' For `alpha_rc_stream()` environments, `entry_treatments` and
#' `entry_families` are set internally from the allocation result and must not
#' be included in the specification.
#'
#' ## Seed-aware replication
#'
#' When `seed_info` and `seed_required_per_plot` are both supplied and a
#' `prep_famoptg()` environment uses `"p_rep"` or `"rcbd_type"` mode,
#' `assign_replication_by_seed()` is called for that environment with a
#' per-environment seed derived from the global `seed` argument. The resulting
#' replication roles (`p_rep_treatments`, `unreplicated_treatments`) are then
#' passed to `prep_famoptg()`. For `"augmented"` mode, all assigned treatments
#' are unreplicated regardless of seed information.
#'
#' ## Family vectors
#'
#' Family labels for allocated treatments are resolved from `treatment_info`
#' when supplied. Treatments absent from `treatment_info` receive the label
#' `"GEN"`. For `alpha_rc_stream()` environments, `entry_families` is
#' populated from these resolved labels before the design function is called.
#'
#' @param treatments Character vector of candidate test treatment IDs. Check
#'   treatments should not be included here; they are specified per environment
#'   in `env_design_specs`.
#'
#' @param environments Character vector of environment names. Must have a
#'   corresponding entry in `env_design_specs` for every element.
#'
#' @param allocation_method Character scalar. Sparse allocation strategy
#'   passed to `allocate_sparse_met()`. Accepted values are
#'   `"random_balanced"`, `"balanced_incomplete"`, `"M3"`, and `"M4"`. See
#'   `allocate_sparse_met()` for the distinction between strategies.
#'
#' @param n_test_entries_per_environment Integer scalar or vector. Number of
#'   test treatments per environment, including common treatments. Passed
#'   directly to `allocate_sparse_met()`.
#'
#' @param target_replications Optional integer scalar. Target number of
#'   environments per non-common treatment. Passed to `allocate_sparse_met()`.
#'   If `NULL`, the allocation function infers a value from available slots.
#'
#' @param common_treatments Optional character vector of treatment IDs forced
#'   into all environments before sparse allocation. Passed to
#'   `allocate_sparse_met()` and `combine_met_fieldbooks()`.
#'
#' @param allocation_group_source Character scalar. Source of genetic group
#'   labels used to guide across-environment allocation. One of `"none"`,
#'   `"Family"`, `"GRM"`, `"A"`. Passed to `allocate_sparse_met()`. Does not
#'   affect within-environment arrangement, which is controlled by
#'   `cluster_source` in each environment's specification.
#'
#' @param env_design_specs Named list of environment-specific design
#'   specifications. Must contain one element per environment in
#'   `environments`. Each element is a named list with at minimum a `design`
#'   field (`"prep_famoptg"` or `"alpha_rc_stream"`). See the
#'   \sQuote{env_design_specs structure} section for full details.
#'
#' @param treatment_info Optional data frame with columns `Treatment` and
#'   `Family`. Required when `allocation_group_source = "Family"`. When
#'   supplied under any mode, used to resolve family labels for allocated
#'   treatments before constructing within-environment layouts. Treatments
#'   absent from this data frame receive family label `"GEN"`.
#'
#' @param GRM Optional numeric matrix. Genomic relationship matrix. Required
#'   when `allocation_group_source = "GRM"`. Row and column names must match
#'   treatment IDs or be reachable through `id_map`.
#'
#' @param A Optional numeric matrix. Pedigree-based numerator relationship
#'   matrix. Required when `allocation_group_source = "A"`.
#'
#' @param id_map Optional data frame with columns `Treatment` and `LineID`.
#'   Required only when `allocation_group_source %in% c("GRM", "A")` and
#'   treatment IDs do not match the row names of the relationship matrix.
#'
#' @param group_method Character scalar. Clustering method for matrix-based
#'   allocation grouping. `"kmeans"` or `"hclust"`. Passed to
#'   `allocate_sparse_met()`. Ignored when
#'   `allocation_group_source %in% c("none", "Family")`.
#'
#' @param group_seed Integer. Seed for k-means initialization in allocation
#'   grouping. Active only when `allocation_group_source %in% c("GRM", "A")`
#'   and `group_method = "kmeans"`.
#'
#' @param group_attempts Integer. Number of random restarts for k-means in
#'   allocation grouping.
#'
#' @param n_pcs_use Integer or `Inf`. Number of leading principal components
#'   used in matrix-based allocation grouping. Ignored when
#'   `allocation_group_source %in% c("none", "Family")`.
#'
#' @param min_groups_per_environment Optional integer. Minimum number of
#'   allocation groups represented in each environment, where feasible. Passed
#'   to `allocate_sparse_met()`.
#'
#' @param min_env_per_group Optional integer. Minimum number of environments
#'   in which each allocation group appears, where feasible. Passed to
#'   `allocate_sparse_met()`.
#'
#' @param balance_groups_across_env Logical. If `TRUE`, sparse allocation
#'   preferentially assigns treatments from under-represented groups. Passed
#'   to `allocate_sparse_met()`.
#'
#' @param force_group_connectivity Logical. If `TRUE`, sparse allocation
#'   preferentially assigns treatments from groups that have not yet reached
#'   `min_env_per_group`. Passed to `allocate_sparse_met()`.
#'
#' @param seed_info Optional data frame with columns `Treatment` and
#'   `SeedAvailable`. When supplied together with `seed_required_per_plot`,
#'   enables seed-aware replication assignment for `prep_famoptg()`
#'   environments operating under `"p_rep"` or `"rcbd_type"` mode.
#'
#' @param seed_required_per_plot Optional data frame with columns `Environment`
#'   and `SeedRequiredPerPlot`. Gives the per-plot seed requirement for each
#'   environment. The value for a given environment is extracted by matching on
#'   the `Environment` column before calling `assign_replication_by_seed()`.
#'
#' @param allow_approximate Logical. Passed to `allocate_sparse_met()`. If
#'   `TRUE`, allows approximate balance when exact balanced incomplete
#'   allocation is infeasible. Default `FALSE`.
#'
#' @param seed Optional integer. Global random seed. Used directly in
#'   `allocate_sparse_met()`. For within-environment design functions, a
#'   per-environment seed is derived as `seed + i` where `i` is the
#'   environment index, ensuring reproducibility while avoiding identical
#'   seeds across environments.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`sparse_allocation`}{Full output of `allocate_sparse_met()`,
#'     including the allocation matrix, long-format allocation, overlap
#'     matrix, line replications, and group summaries.}
#'   \item{`environment_designs`}{Named list with one element per environment.
#'     Each element is the full output object returned by `prep_famoptg()` or
#'     `alpha_rc_stream()` for that environment.}
#'   \item{`combined_field_book`}{Data frame. The combined MET field book
#'     produced by `combine_met_fieldbooks()`, with one row per plot across
#'     all environments and MET-level metadata columns prepended. The
#'     `Gcluster` column is preserved from local field books where present.}
#'   \item{`environment_summary`}{Data frame with one row per environment.
#'     Columns: `Environment`, `LocalDesign`, `ReplicationMode`,
#'     `n_assigned_treatments`, `n_common_treatments`,
#'     `n_replicated_treatments`, `n_unreplicated_treatments`,
#'     `n_total_plots`.}
#'   \item{`group_environment_summary`}{Data frame from
#'     `sparse_allocation$group_by_environment`. Summarizes the number of
#'     assigned treatments from each allocation group in each environment.
#'     `NULL` when `allocation_group_source = "none"`.}
#'   \item{`summary`}{Named list with high-level counts: number of
#'     environments, allocation method, allocation group source, total test
#'     treatments, number of common treatments, total environment assignments,
#'     number of allocation groups, and total rows in the combined field
#'     book.}
#'   \item{`seed_used`}{The integer seed passed to the function, or `NULL` if
#'     none was supplied.}
#' }
#'
#' @examples
#' \dontrun{
#' treatments <- paste0("L", sprintf("%03d", 1:120))
#' environments <- c("Env1", "Env2", "Env3")
#'
#' treatment_info <- data.frame(
#'   Treatment = treatments,
#'   Family    = rep(paste0("F", 1:6), each = 20),
#'   stringsAsFactors = FALSE
#' )
#'
#' ## env_design_specs defines the local design engine and parameters
#' ## for each environment independently.
#' env_design_specs <- list(
#'   Env1 = list(
#'     design               = "prep_famoptg",
#'     replication_mode     = "p_rep",
#'     desired_replications = 2,
#'     check_treatments     = c("CHK1", "CHK2"),
#'     check_families       = c("CHECK", "CHECK"),
#'     n_blocks             = 4,
#'     n_rows               = 10,
#'     n_cols               = 10,
#'     cluster_source       = "Family"
#'   ),
#'   Env2 = list(
#'     design           = "alpha_rc_stream",
#'     check_treatments = c("CHK1", "CHK2"),
#'     check_families   = c("CHECK", "CHECK"),
#'     n_reps           = 2,
#'     n_rows           = 10,
#'     n_cols           = 10,
#'     cluster_source   = "Family"
#'   ),
#'   Env3 = list(
#'     design               = "prep_famoptg",
#'     replication_mode     = "augmented",
#'     check_treatments     = c("CHK1", "CHK2"),
#'     check_families       = c("CHECK", "CHECK"),
#'     n_blocks             = 4,
#'     n_rows               = 10,
#'     n_cols               = 10,
#'     cluster_source       = "Family"
#'   )
#' )
#'
#' ## Full pipeline: random balanced allocation with family-guided group
#' ## spreading, followed by environment-specific local designs.
#' out <- plan_sparse_met_design(
#'   treatments                     = treatments,
#'   environments                   = environments,
#'   allocation_method              = "random_balanced",
#'   n_test_entries_per_environment = 40,
#'   target_replications            = 2,
#'   allocation_group_source        = "Family",
#'   env_design_specs               = env_design_specs,
#'   treatment_info                 = treatment_info,
#'   min_groups_per_environment     = 4,
#'   min_env_per_group              = 2,
#'   seed                           = 123
#' )
#'
#' out$summary
#' head(out$combined_field_book)
#' out$environment_summary
#' head(out$group_environment_summary)
#' }
#'
#' @export
#' 
plan_sparse_met_design <- function(
    treatments,
    environments,
    allocation_method = c("random_balanced", "balanced_incomplete", "M3", "M4"),
    n_test_entries_per_environment,
    target_replications = NULL,
    common_treatments = NULL,
    allocation_group_source = c("none", "Family", "GRM", "A"),
    env_design_specs,
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
    seed_info = NULL,
    seed_required_per_plot = NULL,
    allow_approximate = FALSE,
    seed = NULL
) {
  
  allocation_method <- match.arg(allocation_method)
  allocation_group_source <- match.arg(allocation_group_source)
  group_method <- match.arg(group_method)
  
  seed_used <- seed
  if (!is.null(seed_used)) {
    set.seed(seed_used)
  }
  
  if (!is.list(env_design_specs) || is.null(names(env_design_specs))) {
    stop("`env_design_specs` must be a named list with one element per environment.")
  }
  
  missing_specs <- setdiff(environments, names(env_design_specs))
  if (length(missing_specs) > 0) {
    stop("Missing design specifications for environments: ",
         paste(missing_specs, collapse = ", "))
  }
  
  if (!is.null(treatment_info)) {
    if (!is.data.frame(treatment_info) ||
        !all(c("Treatment", "Family") %in% names(treatment_info))) {
      stop("`treatment_info` must be a data frame with columns `Treatment` and `Family`.")
    }
    treatment_info$Treatment <- as.character(treatment_info$Treatment)
    treatment_info$Family <- as.character(treatment_info$Family)
  }
  
  if (!is.null(seed_info)) {
    if (!is.data.frame(seed_info) ||
        !all(c("Treatment", "SeedAvailable") %in% names(seed_info))) {
      stop("`seed_info` must be a data frame with columns `Treatment` and `SeedAvailable`.")
    }
    seed_info$Treatment <- as.character(seed_info$Treatment)
  }
  
  if (!is.null(seed_required_per_plot)) {
    if (!is.data.frame(seed_required_per_plot) ||
        !all(c("Environment", "SeedRequiredPerPlot") %in% names(seed_required_per_plot))) {
      stop("`seed_required_per_plot` must be a data frame with columns `Environment` and `SeedRequiredPerPlot`.")
    }
    seed_required_per_plot$Environment <- as.character(seed_required_per_plot$Environment)
  }
  
  sparse_out <- allocate_sparse_met(
    treatments = treatments,
    environments = environments,
    allocation_method = allocation_method,
    n_test_entries_per_environment = n_test_entries_per_environment,
    target_replications = target_replications,
    common_treatments = common_treatments,
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
    seed = seed_used
  )
  
  alloc_mat <- sparse_out$allocation_matrix
  
  get_family_vec <- function(trt) {
    if (length(trt) == 0) {
      return(character(0))
    }
    if (is.null(treatment_info)) {
      return(rep("GEN", length(trt)))
    }
    fam <- treatment_info$Family[match(trt, treatment_info$Treatment)]
    fam[is.na(fam)] <- "GEN"
    fam
  }
  
  environment_designs <- vector("list", length(environments))
  names(environment_designs) <- environments
  
  fb_list <- vector("list", length(environments))
  names(fb_list) <- environments
  
  local_designs <- setNames(rep(NA_character_, length(environments)), environments)
  replication_modes <- setNames(rep(NA_character_, length(environments)), environments)
  
  for (i in seq_along(environments)) {
    env <- environments[i]
    spec <- env_design_specs[[env]]
    
    if (is.null(spec$design)) {
      stop("Environment '", env, "' is missing `design`.")
    }
    
    design_name <- spec$design
    local_designs[env] <- design_name
    
    env_treatments <- rownames(alloc_mat)[alloc_mat[, env] == 1L]
    
    if (design_name == "prep_famoptg") {
      
      replication_mode <- if (is.null(spec$replication_mode)) "augmented" else spec$replication_mode
      replication_modes[env] <- replication_mode
      
      desired_replications <- if (is.null(spec$desired_replications)) 2L else spec$desired_replications
      candidate_prep <- spec$candidate_prep
      max_prep <- spec$max_prep
      shortage_action <- if (is.null(spec$shortage_action)) "downgrade" else spec$shortage_action
      
      env_seed_req <- NULL
      if (!is.null(seed_required_per_plot)) {
        env_seed_req <- seed_required_per_plot$SeedRequiredPerPlot[
          match(env, seed_required_per_plot$Environment)
        ]
      }
      
      if (replication_mode %in% c("p_rep", "rcbd_type") &&
          !is.null(seed_info) && !is.null(env_seed_req)) {
        
        rep_assign <- assign_replication_by_seed(
          treatments = env_treatments,
          check_treatments = spec$check_treatments,
          seed_available = seed_info,
          seed_required_per_plot = env_seed_req,
          replication_mode = replication_mode,
          desired_replications = desired_replications,
          candidate_prep = candidate_prep,
          max_prep = max_prep,
          shortage_action = shortage_action,
          seed = if (is.null(seed_used)) NULL else seed_used + i
        )
        
        p_rep_treatments <- rep_assign$p_rep_treatments
        p_rep_reps <- rep_assign$p_rep_reps
        unreplicated_treatments <- rep_assign$unreplicated_treatments
        
      } else if (replication_mode == "augmented") {
        
        p_rep_treatments <- character(0)
        p_rep_reps <- integer(0)
        unreplicated_treatments <- env_treatments
        
      } else if (replication_mode == "rcbd_type") {
        
        p_rep_treatments <- env_treatments
        p_rep_reps <- rep(desired_replications, length(env_treatments))
        unreplicated_treatments <- character(0)
        
      } else if (replication_mode == "p_rep") {
        
        if (is.null(candidate_prep)) {
          candidate_prep <- env_treatments
        }
        candidate_prep <- intersect(candidate_prep, env_treatments)
        
        if (!is.null(max_prep)) {
          candidate_prep <- head(candidate_prep, max_prep)
        }
        
        p_rep_treatments <- candidate_prep
        p_rep_reps <- rep(desired_replications, length(p_rep_treatments))
        unreplicated_treatments <- setdiff(env_treatments, p_rep_treatments)
        
      } else {
        stop("Unknown replication_mode in environment '", env, "': ", replication_mode)
      }
      
      p_rep_families <- get_family_vec(p_rep_treatments)
      unrep_families <- get_family_vec(unreplicated_treatments)
      
      local_args <- spec
      local_args$design <- NULL
      local_args$replication_mode <- NULL
      local_args$desired_replications <- NULL
      local_args$candidate_prep <- NULL
      local_args$max_prep <- NULL
      local_args$shortage_action <- NULL
      
      local_args$p_rep_treatments <- if (length(p_rep_treatments)) p_rep_treatments else character(0)
      local_args$p_rep_reps <- if (length(p_rep_reps)) p_rep_reps else integer(0)
      local_args$p_rep_families <- if (length(p_rep_families)) p_rep_families else character(0)
      local_args$unreplicated_treatments <- if (length(unreplicated_treatments)) unreplicated_treatments else character(0)
      local_args$unreplicated_families <- if (length(unrep_families)) unrep_families else character(0)
      
      local_out <- do.call(prep_famoptg, local_args)
    }
    
    if (design_name == "alpha_rc_stream") {
      replication_modes[env] <- NA_character_
      
      local_args <- spec
      local_args$design <- NULL
      local_args$entry_treatments <- env_treatments
      local_args$entry_families <- get_family_vec(env_treatments)
      
      local_out <- do.call(alpha_rc_stream, local_args)
    }
    
    if (!design_name %in% c("prep_famoptg", "alpha_rc_stream")) {
      stop("Unsupported design for environment '", env, "': ", design_name)
    }
    
    environment_designs[[env]] <- local_out
    fb_list[[env]] <- standardize_local_field_book(local_out$field_book)
  }
  
  combined_field_book <- combine_met_fieldbooks(
    field_books = fb_list,
    local_designs = local_designs,
    replication_modes = replication_modes,
    sparse_method = if (allocation_method %in% c("M3", "random_balanced")) "random_balanced" else "balanced_incomplete",
    common_treatments = common_treatments
  )
  
  # ------------------------------------------------------------
  # Build environment-level summary
  # ------------------------------------------------------------
  env_summary_list <- vector("list", length(environments))
  names(env_summary_list) <- environments
  
  for (env in environments) {
    
    fb <- fb_list[[env]]
    spec <- env_design_specs[[env]]
    
    design_name <- local_designs[[env]]
    replication_mode <- replication_modes[[env]]
    
    check_trt_env <- spec$check_treatments
    if (is.null(check_trt_env)) {
      check_trt_env <- character(0)
    }
    
    assigned_trt <- unique(fb$Treatment[!fb$Treatment %in% check_trt_env])
    n_assigned <- length(assigned_trt)
    
    n_common <- if (is.null(common_treatments)) {
      0L
    } else {
      sum(assigned_trt %in% common_treatments)
    }
    
    trt_counts <- table(fb$Treatment)
    trt_counts <- trt_counts[!names(trt_counts) %in% check_trt_env]
    
    n_replicated <- sum(trt_counts > 1)
    n_unreplicated <- sum(trt_counts == 1)
    
    env_summary_list[[env]] <- data.frame(
      Environment = env,
      LocalDesign = design_name,
      ReplicationMode = replication_mode,
      n_assigned_treatments = n_assigned,
      n_common_treatments = n_common,
      n_replicated_treatments = n_replicated,
      n_unreplicated_treatments = n_unreplicated,
      n_total_plots = nrow(fb),
      stringsAsFactors = FALSE
    )
  }
  
  environment_summary <- do.call(rbind, env_summary_list)
  rownames(environment_summary) <- NULL
  
  group_environment_summary <- sparse_out$group_by_environment
  
  summary <- list(
    n_environments = length(environments),
    allocation_method = if (allocation_method %in% c("M3", "random_balanced")) "random_balanced" else "balanced_incomplete",
    allocation_group_source = allocation_group_source,
    n_total_test_treatments = length(treatments),
    n_common_treatments = if (is.null(common_treatments)) 0L else length(intersect(common_treatments, treatments)),
    total_environment_assignments = sum(sparse_out$environment_sizes),
    n_allocation_groups = if (!is.null(sparse_out$summary$n_groups)) sparse_out$summary$n_groups else NA_integer_,
    combined_n_rows = nrow(combined_field_book)
  )
  
  list(
    sparse_allocation = sparse_out,
    environment_designs = environment_designs,
    combined_field_book = combined_field_book,
    environment_summary = environment_summary,
    group_environment_summary = group_environment_summary,
    summary = summary,
    seed_used = seed_used
  )
}