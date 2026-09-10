#' Build a joint allocation-replication-layout evaluator
#'
#' @description
#' Creates the callback consumed by [optimize_design()] to turn each candidate
#' allocation into actual environment fieldbooks, count its realised integer
#' replication, and derive full local treatment-information matrices. This is
#' the package's joint optimisation bridge: allocation moves are accepted or
#' rejected using the plantable layouts they produce, rather than an assumed
#' scalar site efficiency.
#'
#' @param env_design_specs Named local-design specifications, as in
#'   [plan_sparse_met_design()].
#' @param treatment_info,common_treatments,seed_info,seed_required_per_plot
#'   Inputs passed to [plan_sparse_met_design()].
#' @param minimum_seed_buffer Seed-reserve input passed to
#'   [plan_sparse_met_design()].
#' @param information_models Optional named list, one element per environment,
#'   of arguments for [local_treatment_information()] other than `field_book`,
#'   `treatments`, and `check_treatments`. Unspecified environments use the
#'   defaults of that function.
#' @param cost_per_plot One finite non-negative plot cost, or a named value per
#'   environment, returned to the design objective.
#' @param seed Reproducible local-layout seed.
#' @param plan_args Additional named arguments for [plan_sparse_met_design()].
#'   Core candidate-specific arguments cannot be overridden.
#' @return A function suitable for `optimize_design(design_evaluator = ...)`.
#'   Calling it returns `reps`, `local_information`, `fieldbooks`, `metadata`,
#'   `cost_per_plot`, and `feasible`.
#' @export
fieldbook_design_evaluator <- function(
    env_design_specs, treatment_info = NULL, common_treatments = NULL,
    seed_info = NULL, seed_required_per_plot = NULL,
    minimum_seed_buffer = 0, information_models = NULL,
    cost_per_plot = 1, seed = NULL, plan_args = list()) {
  if (!is.list(env_design_specs) || !length(env_design_specs) ||
      is.null(names(env_design_specs)) || any(names(env_design_specs) == "") ||
      anyDuplicated(names(env_design_specs)))
    stop("`env_design_specs` must be a uniquely named non-empty list.")
  if (!is.null(information_models) &&
      (!is.list(information_models) || is.null(names(information_models)) ||
       anyDuplicated(names(information_models)) ||
       !all(names(information_models) %in% names(env_design_specs))))
    stop("`information_models` must be NULL or a named environment subset.")
  if (!is.numeric(cost_per_plot) || !length(cost_per_plot) ||
      any(!is.finite(cost_per_plot)) || any(cost_per_plot < 0))
    stop("`cost_per_plot` must contain finite non-negative values.")
  if (!is.list(plan_args) || (length(plan_args) && is.null(names(plan_args))))
    stop("`plan_args` must be a named list.")
  locked <- c("treatments", "environments", "allocation_matrix",
              "n_test_entries_per_environment", "env_design_specs",
              "treatment_info", "common_treatments", "seed_info",
              "seed_required_per_plot", "minimum_seed_buffer", "seed")
  bad <- intersect(names(plan_args), locked)
  if (length(bad))
    stop("`plan_args` cannot override: ", paste(bad, collapse = ", "), ".")

  function(allocation_matrix) {
    M <- as.matrix(allocation_matrix)
    envs <- colnames(M); lines <- rownames(M)
    if (is.null(envs) || is.null(lines) ||
        !setequal(envs, names(env_design_specs)))
      stop("Candidate allocation environments must match `env_design_specs`.")
    M <- M[lines, envs, drop = FALSE]
    for (e in envs) {
      checks <- as.character(env_design_specs[[e]]$check_treatments %||%
                               character(0))
      if (length(intersect(lines, checks)))
        stop("Candidate treatments must not duplicate fixed check treatments.")
    }
    call_args <- c(list(
      treatments = lines, environments = envs, allocation_matrix = M,
      n_test_entries_per_environment = as.integer(colSums(M)),
      env_design_specs = env_design_specs,
      treatment_info = treatment_info,
      common_treatments = common_treatments,
      seed_info = seed_info,
      seed_required_per_plot = seed_required_per_plot,
      minimum_seed_buffer = minimum_seed_buffer,
      seed = seed), plan_args)
    plan <- do.call(plan_sparse_met_design, call_args)

    reps <- matrix(0, nrow(M), ncol(M), dimnames = dimnames(M))
    local_info <- vector("list", length(envs)); names(local_info) <- envs
    fieldbooks <- vector("list", length(envs)); names(fieldbooks) <- envs
    fixed_overhead <- stats::setNames(numeric(length(envs)), envs)
    for (e in envs) {
      fb <- plan$environment_designs[[e]]$field_book
      fieldbooks[[e]] <- fb
      counts <- table(as.character(fb$Treatment))
      present <- lines[M[, e] == 1L]
      reps[present, e] <- as.numeric(counts[present])
      if (anyNA(reps[present, e]) || any(reps[present, e] <= 0))
        stop("The local fieldbook for `", e,
             "` omitted an allocated treatment.")
      occupied <- sum(!is.na(fb$Treatment) & nzchar(as.character(fb$Treatment)))
      fixed_overhead[e] <- occupied - sum(reps[, e])
      if (fixed_overhead[e] < 0)
        stop("Local plot accounting is inconsistent for `", e, "`.")
      model <- if (is.null(information_models) ||
                   is.null(information_models[[e]])) list() else
        information_models[[e]]
      if (!is.list(model))
        stop("Every information model must be a list.")
      forbidden <- intersect(names(model),
                             c("field_book", "treatments", "check_treatments"))
      if (length(forbidden))
        stop("Information models cannot override field_book, treatments, or checks.")
      checks <- as.character(env_design_specs[[e]]$check_treatments %||%
                               character(0))
      local_info[[e]] <- do.call(local_treatment_information,
        c(list(field_book = fb, treatments = present,
               check_treatments = checks), model))
    }
    list(
      reps = reps,
      local_information = local_info,
      cost_per_plot = cost_per_plot,
      fixed_plot_overhead = fixed_overhead,
      feasible = TRUE,
      fieldbooks = fieldbooks,
      metadata = list(environment_summary = plan$environment_summary,
                      combined_field_book = plan$combined_field_book,
                      seed_ledger = plan$seed_ledger,
                      plan_summary = plan$summary))
  }
}
