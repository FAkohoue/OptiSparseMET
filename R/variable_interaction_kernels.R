#' Build dedicated variable-level environmental interaction kernels
#'
#' @description
#' Constructs auditable kernels for pre-specified interactions among individual
#' weather, soil, management, or geographic variables. This function complements
#' the omnibus modality interactions produced by
#' [add_environment_kernel_interactions()]. An omnibus weather-by-soil kernel
#' represents the complete tensor product of two modality feature spaces. In
#' contrast, a dedicated variable kernel guarantees that a named hypothesis,
#' such as temperature by root-zone bulk density, is represented separately.
#'
#' The procedure has five stages. First, every requested variable is matched to
#' a quality-controlled covariate column. Second, the parent feature blocks are
#' centred and standardised. Third, their row-wise tensor product is formed.
#' Fourth, the tensor-product features are residualised against the intercept
#' and all parent main-effect features by default. Finally, Gaussian
#' bandwidth-ensemble kernels are constructed for the parent and residual
#' interaction features. Residualisation reduces main-effect leakage while
#' preserving a positive-semidefinite kernel because it is applied before the
#' Gaussian kernel is constructed.
#'
#' Strong heredity is the default: every dedicated interaction is returned with
#' the exact parent-variable kernels used to construct it. Strong heredity here
#' is a model-structure rule, not a claim that every parent has a non-zero
#' historical covariance weight. The parent kernels remain in every model in
#' which their interaction is assessed.
#'
#' @param covariates Named list of environment-by-variable matrices or data
#'   frames. The `covariates` component returned by
#'   [build_environment_kernels()] can be supplied directly. Every block must
#'   have the same unique environment row names.
#' @param interactions A long data frame or a named list defining the
#'   hypotheses. The recommended data-frame format has four required columns:
#'   `interaction`, `parent`, `modality`, and `variable`. Each row identifies
#'   one covariate column. Rows with the same `interaction` and `parent` define
#'   one multicolumn parent feature, for example the indicators of a categorical
#'   management factor. Each interaction must contain at least two parents.
#'
#'   The equivalent list format is
#'   `list(temp_x_density = list(temperature = list(modality = "weather",
#'   variables = "mean_temperature"), density = list(modality = "soil",
#'   variables = "bulk_density")))`.
#' @param bandwidth_multipliers Positive multipliers applied to the median
#'   positive pairwise distance. The resulting Gaussian kernels are averaged.
#'   The default `c(0.5, 1, 2)` reduces dependence on one bandwidth in a small
#'   environment network.
#' @param orthogonalise Logical. If `TRUE`, residualise each tensor-product
#'   feature against the intercept and its parent main-effect features.
#' @param hierarchy Either `"strong"` or `"none"`. `"strong"` returns every
#'   dedicated parent kernel and records the parent-child hierarchy.
#'   `"none"` returns only the interaction kernels and should be used only when
#'   equivalent parent main effects are already guaranteed elsewhere.
#' @param max_tensor_columns Maximum number of raw tensor-product columns
#'   permitted for one interaction. This explicit guard prevents an unnoticed
#'   expansion of high-cardinality categorical interactions.
#' @param min_residual_fraction Minimum retained sum-of-squares fraction after
#'   residualisation. A smaller value indicates that the proposed interaction
#'   is almost completely explained by its parent main effects and is rejected
#'   as numerically uninformative.
#'
#' @return A list containing:
#' \describe{
#'   \item{`kernels`}{Dedicated parent and interaction kernels.}
#'   \item{`audit`}{One row per interaction, including feature dimensions,
#'     residual information, bandwidths, effective rank, scope, and hierarchy.}
#'   \item{`variable_ledger`}{One row per requested covariate, documenting its
#'     modality, parent, matched column, and resulting parent kernel.}
#'   \item{`hierarchy`}{A named list mapping each interaction kernel to its
#'     dedicated parent kernels.}
#'   \item{`definitions`}{The normalised long-form interaction specification.}
#'   \item{`environments`}{The realised environment order.}
#' }
#'
#' @section Interpretation:
#' Inclusion in `kernels` means that the named interaction was represented
#' structurally. It does not demonstrate that the interaction explains genetic
#' environment covariance and does not establish causality. Use
#' [assess_variable_interactions()] for blocked historical validation. Without
#' valid historical MET responses, retain each kernel as a separate unweighted
#' maximin sensitivity scenario.
#'
#' @examples
#' env <- paste0("E", 1:8)
#' covariates <- list(
#'   weather = matrix(
#'     c(seq(22, 29), seq(300, 650, length.out = 8)),
#'     nrow = 8,
#'     dimnames = list(env, c("mean_temperature", "rainfall"))
#'   ),
#'   soil = matrix(
#'     c(seq(1.1, 1.5, length.out = 8), seq(5.2, 6.6, length.out = 8)),
#'     nrow = 8,
#'     dimnames = list(env, c("bulk_density", "pH"))
#'   )
#' )
#' hypotheses <- data.frame(
#'   interaction = c("temperature_x_density", "temperature_x_density"),
#'   parent = c("temperature", "density"),
#'   modality = c("weather", "soil"),
#'   variable = c("mean_temperature", "bulk_density")
#' )
#' out <- build_variable_interaction_kernels(covariates, hypotheses)
#' names(out$kernels)
#' out$audit
#' @export
build_variable_interaction_kernels <- function(
    covariates, interactions,
    bandwidth_multipliers = c(0.5, 1, 2),
    orthogonalise = TRUE,
    hierarchy = c("strong", "none"),
    max_tensor_columns = 200L,
    min_residual_fraction = 1e-6) {
  hierarchy <- match.arg(hierarchy)
  if (!is.logical(orthogonalise) || length(orthogonalise) != 1L ||
      is.na(orthogonalise))
    stop("`orthogonalise` must be TRUE or FALSE.")
  if (!is.numeric(bandwidth_multipliers) ||
      any(!is.finite(bandwidth_multipliers)) ||
      any(bandwidth_multipliers <= 0))
    stop("`bandwidth_multipliers` must be finite and positive.")
  bandwidth_multipliers <- sort(unique(as.numeric(bandwidth_multipliers)))
  if (!is.numeric(max_tensor_columns) || length(max_tensor_columns) != 1L ||
      !is.finite(max_tensor_columns) || max_tensor_columns < 1)
    stop("`max_tensor_columns` must be a positive integer.")
  max_tensor_columns <- as.integer(max_tensor_columns)
  if (!is.numeric(min_residual_fraction) ||
      length(min_residual_fraction) != 1L ||
      !is.finite(min_residual_fraction) ||
      min_residual_fraction <= 0 || min_residual_fraction >= 1)
    stop("`min_residual_fraction` must be one number in (0, 1).")

  C <- .validate_variable_covariates(covariates)
  envs <- rownames(C[[1L]])
  specification <- .normalise_variable_interactions(interactions)
  unknown_modalities <- setdiff(unique(specification$modality), names(C))
  if (length(unknown_modalities))
    stop("Interaction definitions name absent modalities: ",
         paste(unknown_modalities, collapse = ", "), ".")

  interaction_order <- unique(specification$interaction)
  kernels <- list()
  audit_rows <- list()
  ledger_rows <- list()
  hierarchy_map <- list()
  parent_cache <- list()
  parent_names <- character()

  for (interaction_label in interaction_order) {
    definition <- specification[
      specification$interaction == interaction_label, , drop = FALSE
    ]
    parent_order <- unique(definition$parent)
    if (length(parent_order) < 2L)
      stop("Interaction '", interaction_label,
           "' must contain at least two distinct parents.")

    parent_features <- list()
    realised_parent_names <- character(length(parent_order))
    parent_modalities <- character(length(parent_order))
    parent_variables <- vector("list", length(parent_order))

    for (j in seq_along(parent_order)) {
      parent_label <- parent_order[j]
      part <- definition[definition$parent == parent_label, , drop = FALSE]
      modalities <- unique(part$modality)
      if (length(modalities) != 1L)
        stop("Parent '", parent_label, "' in interaction '",
             interaction_label, "' must belong to one modality.")
      modality <- modalities[[1L]]
      variables <- unique(part$variable)
      absent <- setdiff(variables, colnames(C[[modality]]))
      if (length(absent))
        stop("Interaction '", interaction_label, "', parent '",
             parent_label, "' requests variables absent after environmental ",
             "quality control in modality '", modality, "': ",
             paste(absent, collapse = ", "), ".")
      X <- C[[modality]][, variables, drop = FALSE]
      Z <- .standardise_interaction_features(
        X, context = paste0("parent '", parent_label,
                            "' of interaction '", interaction_label, "'")
      )
      parent_features[[parent_label]] <- Z
      parent_modalities[j] <- modality
      parent_variables[[j]] <- variables

      cache_key <- paste(
        modality, paste(sort(variables), collapse = "\034"), sep = "\035"
      )
      kernel_name <- paste0(
        "variable_main__", .kernel_safe_id(modality), "__",
        .kernel_safe_id(paste(sort(variables), collapse = "_plus_"))
      )
      if (cache_key %in% names(parent_cache)) {
        kernel_name <- parent_cache[[cache_key]]
      } else {
        if (kernel_name %in% parent_names)
          stop("Two distinct parent definitions produce the same safe kernel ",
               "name '", kernel_name, "'. Rename the colliding variables.")
        parent_kernel <- .gaussian_feature_ensemble(
          Z, bandwidth_multipliers,
          context = paste0("parent '", parent_label, "'")
        )
        attr(parent_kernel, "kernel_role") <- "variable_main"
        attr(parent_kernel, "variable_modality") <- modality
        attr(parent_kernel, "variable_columns") <- variables
        attr(parent_kernel, "bandwidths") <-
          attr(parent_kernel, "realised_bandwidths")
        parent_cache[[cache_key]] <- kernel_name
        parent_names <- c(parent_names, kernel_name)
        if (hierarchy == "strong")
          kernels[[kernel_name]] <- parent_kernel
      }
      realised_parent_names[j] <- kernel_name
      ledger_rows[[length(ledger_rows) + 1L]] <- data.frame(
        interaction = interaction_label,
        parent = parent_label,
        modality = modality,
        variable = variables,
        matched_after_qc = TRUE,
        parent_kernel = kernel_name,
        stringsAsFactors = FALSE
      )
    }
    if (anyDuplicated(realised_parent_names))
      stop("Interaction '", interaction_label,
           "' defines two parents with the same modality-variable feature ",
           "block. Merge the duplicate parents or revise the hypothesis.")

    raw_column_count <- prod(vapply(parent_features, ncol, integer(1)))
    if (!is.finite(raw_column_count) ||
        raw_column_count > max_tensor_columns)
      stop("Interaction '", interaction_label, "' expands to ",
           format(raw_column_count, scientific = FALSE),
           " tensor columns, exceeding `max_tensor_columns = ",
           max_tensor_columns, "`.")
    tensor <- .row_tensor_product(parent_features)
    tensor <- .standardise_interaction_features(
      tensor, context = paste0("raw tensor for interaction '",
                               interaction_label, "'"),
      allow_drop = TRUE
    )
    retained_before_residualisation <- ncol(tensor)

    residual_fraction <- 1
    if (orthogonalise) {
      parent_design <- do.call(cbind, parent_features)
      parent_design <- cbind(intercept = 1, parent_design)
      parent_qr <- qr(parent_design)
      residual <- qr.resid(parent_qr, tensor)
      numerator <- sum(residual^2)
      denominator <- sum(tensor^2)
      residual_fraction <- if (denominator > 0) numerator / denominator else 0
      if (!is.finite(residual_fraction) ||
          residual_fraction < min_residual_fraction)
        stop("Interaction '", interaction_label,
             "' retains only ", signif(residual_fraction, 4),
             " of its tensor-feature information after parent-effect ",
             "residualisation, below `min_residual_fraction = ",
             min_residual_fraction, "`. The interaction is aliased with its ",
             "parents or unsupported by the environment sample.")
      tensor <- .standardise_interaction_features(
        residual, context = paste0("residual tensor for interaction '",
                                   interaction_label, "'"),
        allow_drop = TRUE
      )
    }

    interaction_kernel <- .gaussian_feature_ensemble(
      tensor, bandwidth_multipliers,
      context = paste0("interaction '", interaction_label, "'")
    )
    interaction_name <- paste0(
      "variable_x__", .kernel_safe_id(interaction_label)
    )
    if (interaction_name %in% c(names(kernels), parent_names))
      stop("Interaction names are not unique after safe-name conversion: '",
           interaction_name, "'.")
    attr(interaction_kernel, "kernel_role") <- "interaction"
    attr(interaction_kernel, "interaction_class") <- "variable_level"
    attr(interaction_kernel, "interaction_label") <- interaction_label
    attr(interaction_kernel, "interaction_parents") <-
      realised_parent_names
    attr(interaction_kernel, "interaction_parent_labels") <- parent_order
    attr(interaction_kernel, "interaction_mode") <-
      if (orthogonalise) "residualised_variable_tensor" else
        "variable_tensor"
    attr(interaction_kernel, "interaction_scope") <-
      if (length(unique(parent_modalities)) == 1L)
        "within_modality" else "cross_modality"
    attr(interaction_kernel, "strong_heredity") <- hierarchy == "strong"
    attr(interaction_kernel, "bandwidths") <-
      attr(interaction_kernel, "realised_bandwidths")
    kernels[[interaction_name]] <- interaction_kernel
    hierarchy_map[[interaction_name]] <- realised_parent_names

    eigenvalues <- eigen(
      interaction_kernel, symmetric = TRUE, only.values = TRUE
    )$values
    rank_tolerance <- max(abs(eigenvalues)) *
      max(nrow(interaction_kernel), 1L) * .Machine$double.eps^0.5
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      interaction = interaction_label,
      kernel = interaction_name,
      order = length(parent_order),
      scope = attr(interaction_kernel, "interaction_scope"),
      modalities = paste(parent_modalities, collapse = " x "),
      variables = paste(
        vapply(parent_variables, paste, collapse = " + ",
               FUN.VALUE = character(1)),
        collapse = " x "
      ),
      parent_kernels = paste(realised_parent_names, collapse = ";"),
      raw_tensor_columns = raw_column_count,
      informative_tensor_columns = retained_before_residualisation,
      final_tensor_columns = ncol(tensor),
      orthogonalised = orthogonalise,
      residual_information_fraction = residual_fraction,
      effective_kernel_rank = sum(eigenvalues > rank_tolerance),
      bandwidths = paste(
        signif(attr(interaction_kernel, "realised_bandwidths"), 6),
        collapse = ";"
      ),
      hierarchy = hierarchy,
      status = "constructed_structural_hypothesis",
      stringsAsFactors = FALSE
    )
  }

  list(
    kernels = kernels,
    audit = do.call(rbind, audit_rows),
    variable_ledger = do.call(rbind, ledger_rows),
    hierarchy = hierarchy_map,
    definitions = specification,
    environments = envs
  )
}

#' Assess dedicated environmental interactions with historical MET evidence
#'
#' @description
#' Evaluates each dedicated variable-level interaction separately against a
#' common main-effect model. The comparison uses the blocked validation
#' machinery in [calibrate_environment_covariance()]. Leave-target-out folds
#' are used when several year-specific targets are available; otherwise
#' leave-environment-out folds are used. One-sided paired sign tests determine
#' whether adding an interaction reduces prediction error consistently.
#' Multiplicity is controlled across the assessed interactions by Holm's method
#' by default.
#'
#' Without valid historical multi-environment trial (MET) responses, no
#' interaction is declared supported or unsupported. Every interaction is
#' labelled `structural_uncertainty` and should remain a separate, unweighted
#' maximin scenario.
#'
#' @param kernels Named list containing main and dedicated interaction kernels,
#'   usually the combined output of [build_environment_kernels()] and
#'   [build_variable_interaction_kernels()].
#' @param target Optional empirical environment genetic-correlation matrix or
#'   list of year-specific matrices.
#' @param historical Optional long data frame of historical adjusted genetic
#'   values. Supply only one of `target` and `historical`.
#' @param genotype_col,environment_col,value_col Column names in `historical`.
#' @param year_col Optional year column. When present, a correlation target is
#'   constructed for each year and leave-year-out validation is preferred.
#' @param interaction_kernels Optional character vector naming the dedicated
#'   interactions. By default, kernels carrying
#'   `interaction_class = "variable_level"` are selected.
#' @param alpha Family-wise significance level for interaction support.
#' @param p_adjust Multiplicity adjustment passed to [stats::p.adjust()].
#'   Holm's procedure is the default. Use `"none"` only for a pre-registered
#'   single hypothesis or an explicitly exploratory analysis.
#' @param n_boot Number of pair-bootstrap calibration draws used to report
#'   component inclusion stability.
#' @param seed Optional reproducibility seed.
#'
#' @return A list with an interaction-level `evidence` table, character vectors
#'   `supported` and `structural_uncertainty`, the fitted single-interaction
#'   models, and the historical targets used for validation. A supported
#'   interaction has positive mean blocked-validation improvement, an adjusted
#'   one-sided sign-test probability below `alpha`, and a positive fitted
#'   covariance contribution.
#'
#' @section Interpretation:
#' Statistical support means that the dedicated kernel improved prediction of
#' historical genetic environment correlations under the available blocked
#' validation. It is not causal evidence. An unsupported interaction is
#' excluded from the central calibrated covariance but remains a legitimate
#' structural sensitivity scenario.
#'
#' @examples
#' \dontrun{
#' evidence <- assess_variable_interactions(
#'   kernels = c(environment_kernels, variable_kernels),
#'   historical = adjusted_historical_met,
#'   year_col = "year",
#'   n_boot = 500,
#'   seed = 2026
#' )
#' evidence$evidence
#' }
#' @export
assess_variable_interactions <- function(
    kernels, target = NULL, historical = NULL,
    genotype_col = "genotype", environment_col = "environment",
    value_col = "value", year_col = NULL,
    interaction_kernels = NULL,
    alpha = 0.05, p_adjust = "holm",
    n_boot = 200L, seed = NULL) {
  K <- .validate_environment_kernels(kernels)
  if (!is.null(target) && !is.null(historical))
    stop("Supply only one of `target` and `historical`.")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("`alpha` must be one finite number in (0, 1).")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      is.na(p_adjust) || !p_adjust %in% stats::p.adjust.methods)
    stop("`p_adjust` must be one of `stats::p.adjust.methods`.")
  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) ||
      n_boot < 0)
    stop("`n_boot` must be a non-negative integer.")
  n_boot <- as.integer(n_boot)

  original <- kernels
  if (is.null(interaction_kernels)) {
    is_variable <- vapply(original, function(M)
      identical(attr(M, "interaction_class"), "variable_level"), logical(1))
    interaction_kernels <- names(original)[is_variable]
  }
  interaction_kernels <- as.character(interaction_kernels)
  if (!length(interaction_kernels))
    stop("No dedicated variable-level interaction kernels were identified.")
  if (anyNA(interaction_kernels) || any(!nzchar(interaction_kernels)) ||
      anyDuplicated(interaction_kernels) ||
      any(!interaction_kernels %in% names(K)))
    stop("`interaction_kernels` must uniquely name kernels in `kernels`.")

  all_interactions <- names(K)[vapply(names(K), function(nm) {
    identical(attr(original[[nm]], "kernel_role"), "interaction") ||
      grepl("_x_", nm, fixed = TRUE)
  }, logical(1))]
  main_names <- setdiff(names(K), all_interactions)
  if (!length(main_names))
    stop("At least one main-effect kernel is required.")

  envs <- rownames(K[[1L]])
  targets <- target
  if (!is.null(historical)) {
    if (!is.null(year_col)) {
      if (!year_col %in% names(historical))
        stop("`historical` has no `", year_col, "` column.")
      historical_split <- split(historical, historical[[year_col]])
      targets <- lapply(
        historical_split, .historical_environment_correlation,
        envs = envs, genotype_col = genotype_col,
        environment_col = environment_col, value_col = value_col
      )
    } else {
      targets <- .historical_environment_correlation(
        historical, envs, genotype_col, environment_col, value_col
      )
    }
  }

  if (is.null(targets)) {
    evidence <- do.call(rbind, lapply(interaction_kernels, function(nm) {
      data.frame(
        interaction = nm,
        label = .interaction_label(original[[nm]], nm),
        scope = .interaction_attribute(original[[nm]], "interaction_scope"),
        parent_kernels = paste(
          attr(original[[nm]], "interaction_parents"), collapse = ";"
        ),
        validation_basis = NA_character_,
        mean_rmse_improvement = NA_real_,
        positive_folds = NA_integer_,
        non_tied_folds = NA_integer_,
        sign_test_p = NA_real_,
        adjusted_p = NA_real_,
        fitted_weight = NA_real_,
        bootstrap_inclusion = NA_real_,
        supported = NA,
        decision = "structural_uncertainty",
        reason = paste0(
          "No historical MET target was available; no response-derived ",
          "interaction decision was made."
        ),
        stringsAsFactors = FALSE
      )
    }))
    return(list(
      evidence = evidence,
      supported = character(),
      structural_uncertainty = interaction_kernels,
      fits = NULL,
      target = NULL
    ))
  }

  fits <- vector("list", length(interaction_kernels))
  names(fits) <- interaction_kernels
  rows <- vector("list", length(interaction_kernels))
  for (i in seq_along(interaction_kernels)) {
    nm <- interaction_kernels[i]
    subset_names <- unique(c(main_names, nm))
    fit_seed <- if (is.null(seed)) NULL else as.integer(seed) + i - 1L
    fit <- calibrate_environment_covariance(
      kernels = K[subset_names],
      target = targets,
      interaction_kernels = nm,
      interaction_policy = "include",
      n_boot = n_boot,
      seed = fit_seed
    )
    fits[[nm]] <- fit
    paired <- .paired_interaction_validation(fit$cross_validation)
    fitted_weight <- unname(fit$weights[nm])
    bootstrap_inclusion <- if (is.null(fit$bootstrap_weights)) NA_real_ else
      mean(fit$bootstrap_weights[, nm] > sqrt(.Machine$double.eps))
    rows[[i]] <- data.frame(
      interaction = nm,
      label = .interaction_label(original[[nm]], nm),
      scope = .interaction_attribute(original[[nm]], "interaction_scope"),
      parent_kernels = paste(
        attr(original[[nm]], "interaction_parents"), collapse = ";"
      ),
      validation_basis = paired$basis,
      mean_rmse_improvement = paired$mean_improvement,
      positive_folds = paired$positive_folds,
      non_tied_folds = paired$non_tied_folds,
      sign_test_p = paired$sign_p,
      adjusted_p = NA_real_,
      fitted_weight = fitted_weight,
      bootstrap_inclusion = bootstrap_inclusion,
      supported = FALSE,
      decision = "retain_as_structural_uncertainty",
      reason = "Historical support has not yet been evaluated.",
      stringsAsFactors = FALSE
    )
  }
  evidence <- do.call(rbind, rows)
  finite_p <- is.finite(evidence$sign_test_p)
  evidence$adjusted_p[finite_p] <- stats::p.adjust(
    evidence$sign_test_p[finite_p], method = p_adjust
  )
  evidence$supported <- with(
    evidence,
    is.finite(mean_rmse_improvement) & mean_rmse_improvement > 0 &
      is.finite(adjusted_p) & adjusted_p < alpha &
      is.finite(fitted_weight) &
      fitted_weight > sqrt(.Machine$double.eps)
  )
  evidence$decision[evidence$supported] <-
    "eligible_for_central_covariance"
  evidence$reason <- ifelse(
    evidence$supported,
    paste0(
      "The interaction improved blocked prediction, passed the adjusted ",
      "paired sign test, and received a positive covariance contribution."
    ),
    paste0(
      "The interaction did not satisfy every response-evidence requirement; ",
      "exclude it from the central covariance and retain it as a structural ",
      "uncertainty scenario."
    )
  )
  rownames(evidence) <- NULL
  list(
    evidence = evidence,
    supported = evidence$interaction[evidence$supported],
    structural_uncertainty =
      evidence$interaction[!evidence$supported],
    fits = fits,
    target = targets
  )
}

.validate_variable_covariates <- function(covariates) {
  if (!is.list(covariates) || !length(covariates) ||
      is.null(names(covariates)) || anyNA(names(covariates)) ||
      any(!nzchar(names(covariates))) || anyDuplicated(names(covariates)))
    stop("`covariates` must be a non-empty uniquely named list.")
  out <- lapply(covariates, function(x) {
    d <- as.data.frame(
      x, check.names = FALSE, stringsAsFactors = FALSE
    )
    if (any(!vapply(d, is.numeric, logical(1))))
      stop("Dedicated interaction covariates must already be numeric. ",
           "Use the quality-controlled `covariates` output from ",
           "`build_environment_kernels()` so categorical variables are ",
           "expanded consistently.")
    M <- data.matrix(d)
    rownames(M) <- rownames(x)
    M
  })
  ref <- rownames(out[[1L]])
  if (is.null(ref) || anyNA(ref) || any(!nzchar(ref)) || anyDuplicated(ref))
    stop("Every covariate block needs unique, non-missing environment row names.")
  for (nm in names(out)) {
    M <- out[[nm]]
    if (is.null(rownames(M)) || !setequal(ref, rownames(M)))
      stop("Covariate block '", nm,
           "' does not cover the same environments as the first block.")
    M <- M[ref, , drop = FALSE]
    if (!ncol(M) || is.null(colnames(M)) || anyNA(colnames(M)) ||
        any(!nzchar(colnames(M))) || anyDuplicated(colnames(M)))
      stop("Covariate block '", nm,
           "' needs unique, non-empty variable column names.")
    if (any(!is.finite(M)))
      stop("Covariate block '", nm,
           "' contains non-finite values after quality control.")
    out[[nm]] <- M
  }
  out
}

.normalise_variable_interactions <- function(interactions) {
  if (is.data.frame(interactions)) {
    specification <- as.data.frame(
      interactions, check.names = FALSE, stringsAsFactors = FALSE
    )
    required <- c("interaction", "parent", "modality", "variable")
    if (!all(required %in% names(specification)))
      stop("Data-frame `interactions` needs columns: ",
           paste(required, collapse = ", "), ".")
    specification <- specification[, required, drop = FALSE]
  } else if (is.list(interactions) && length(interactions)) {
    if (is.null(names(interactions)) || anyNA(names(interactions)) ||
        any(!nzchar(names(interactions))) ||
        anyDuplicated(names(interactions)))
      stop("List `interactions` must have unique, non-empty names.")
    rows <- list()
    for (interaction_name in names(interactions)) {
      parents <- interactions[[interaction_name]]
      if (!is.list(parents) || length(parents) < 2L ||
          is.null(names(parents)) || anyNA(names(parents)) ||
          any(!nzchar(names(parents))) || anyDuplicated(names(parents)))
        stop("Interaction '", interaction_name,
             "' must be a named list of at least two parents.")
      for (parent_name in names(parents)) {
        parent <- parents[[parent_name]]
        if (is.character(parent) && !is.null(names(parent)) &&
            length(unique(names(parent))) == 1L) {
          modality <- unique(names(parent))
          variables <- unname(parent)
        } else {
          if (!is.list(parent) ||
              !"modality" %in% names(parent) ||
              !any(c("variable", "variables") %in% names(parent)))
            stop("Parent '", parent_name, "' in interaction '",
                 interaction_name, "' needs `modality` and `variables`.")
          modality <- as.character(parent$modality)
          variables <- if (!is.null(parent$variables))
            as.character(parent$variables) else
              as.character(parent$variable)
        }
        if (length(modality) != 1L)
          stop("Every interaction parent must name one modality.")
        rows[[length(rows) + 1L]] <- data.frame(
          interaction = interaction_name,
          parent = parent_name,
          modality = modality,
          variable = variables,
          stringsAsFactors = FALSE
        )
      }
    }
    specification <- do.call(rbind, rows)
  } else {
    stop("`interactions` must be a non-empty data frame or named list.")
  }
  for (nm in names(specification))
    specification[[nm]] <- as.character(specification[[nm]])
  if (!nrow(specification) || anyNA(specification) ||
      any(!nzchar(as.matrix(specification))))
    stop("Interaction definitions cannot contain missing or empty values.")
  duplicate_rows <- duplicated(specification)
  if (any(duplicate_rows))
    specification <- specification[!duplicate_rows, , drop = FALSE]
  parent_counts <- vapply(
    split(specification$parent, specification$interaction),
    function(z) length(unique(z)), integer(1)
  )
  if (any(parent_counts < 2L))
    stop("Every interaction must define at least two distinct parents.")
  rownames(specification) <- NULL
  specification
}

.standardise_interaction_features <- function(
    X, context, allow_drop = FALSE) {
  X <- as.matrix(X)
  if (!ncol(X)) stop("No columns were supplied for ", context, ".")
  scales <- apply(X, 2L, stats::sd)
  informative <- is.finite(scales) & scales > sqrt(.Machine$double.eps)
  if (!all(informative) && !allow_drop)
    stop("Non-informative variables remain in ", context, ": ",
         paste(colnames(X)[!informative], collapse = ", "), ".")
  X <- X[, informative, drop = FALSE]
  if (!ncol(X))
    stop("No informative feature remains in ", context, ".")
  Z <- scale(X)
  if (any(!is.finite(Z)))
    stop("Non-finite standardised values were produced for ", context, ".")
  Z <- as.matrix(Z)
  attr(Z, "scaled:center") <- NULL
  attr(Z, "scaled:scale") <- NULL
  Z
}

.row_tensor_product <- function(features) {
  out <- matrix(
    1, nrow(features[[1L]]), 1L,
    dimnames = list(rownames(features[[1L]]), "tensor")
  )
  for (parent_name in names(features)) {
    P <- features[[parent_name]]
    next_out <- matrix(
      NA_real_, nrow(out), ncol(out) * ncol(P),
      dimnames = list(rownames(out), NULL)
    )
    next_names <- character(ncol(next_out))
    k <- 0L
    for (a in seq_len(ncol(out))) {
      for (b in seq_len(ncol(P))) {
        k <- k + 1L
        next_out[, k] <- out[, a] * P[, b]
        next_names[k] <- paste(colnames(out)[a], parent_name,
                               colnames(P)[b], sep = ":")
      }
    }
    colnames(next_out) <- next_names
    out <- next_out
  }
  out
}

.gaussian_feature_ensemble <- function(
    Z, bandwidth_multipliers, context) {
  Z <- as.matrix(Z)
  distances <- as.matrix(stats::dist(Z))
  positive <- distances[upper.tri(distances) & distances > 0]
  base_bandwidth <- stats::median(positive, na.rm = TRUE)
  if (!is.finite(base_bandwidth) || base_bandwidth <= 0)
    stop("All environment feature profiles are identical for ", context, ".")
  bandwidths <- base_bandwidth * bandwidth_multipliers
  squared <- distances^2
  kernels <- lapply(bandwidths, function(h)
    exp(-squared / (2 * h^2)))
  out <- Reduce(`+`, kernels) / length(kernels)
  dimnames(out) <- list(rownames(Z), rownames(Z))
  out <- .normalise_environment_kernel(.bend_pd(out))
  attr(out, "realised_bandwidths") <- bandwidths
  out
}

.kernel_safe_id <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) stop("A kernel label contains no usable characters.")
  out
}

.paired_interaction_validation <- function(cv) {
  empty <- list(
    basis = NA_character_, mean_improvement = NA_real_,
    positive_folds = 0L, non_tied_folds = 0L, sign_p = NA_real_
  )
  if (is.null(cv) || !nrow(cv)) return(empty)
  basis <- if (any(cv$validation == "leave_target_out"))
    "leave_target_out" else if (any(cv$validation == "leave_environment_out"))
      "leave_environment_out" else return(empty)
  main <- cv[
    cv$model == "main_only" & cv$validation == basis,
    c("held_out", "rmse"), drop = FALSE
  ]
  interaction <- cv[
    cv$model == "main_plus_interactions" & cv$validation == basis,
    c("held_out", "rmse"), drop = FALSE
  ]
  paired <- merge(
    main, interaction, by = "held_out",
    suffixes = c("_main", "_interaction")
  )
  improvement <- paired$rmse_main - paired$rmse_interaction
  improvement <- improvement[is.finite(improvement)]
  if (!length(improvement)) return(empty)
  tolerance <- sqrt(.Machine$double.eps)
  signs <- improvement[abs(improvement) > tolerance]
  non_tied <- length(signs)
  positive <- sum(signs > 0)
  list(
    basis = basis,
    mean_improvement = mean(improvement),
    positive_folds = positive,
    non_tied_folds = non_tied,
    sign_p = if (non_tied)
      stats::pbinom(positive - 1L, non_tied, 0.5,
                    lower.tail = FALSE) else NA_real_
  )
}

.interaction_label <- function(K, fallback) {
  label <- attr(K, "interaction_label")
  if (is.null(label) || !length(label)) fallback else as.character(label)[1L]
}

.interaction_attribute <- function(K, attribute) {
  value <- attr(K, attribute)
  if (is.null(value) || !length(value)) NA_character_ else
    as.character(value)[1L]
}
