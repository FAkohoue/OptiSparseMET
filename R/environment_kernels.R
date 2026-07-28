#' Build separate environmental kernels by information modality
#'
#' @description
#' Constructs weather, soil, management, and geographic kernels separately.
#' This prevents a high-dimensional source (usually weather) from receiving
#' accidental extra weight merely because it has more columns. Continuous
#' blocks use an ensemble of Gaussian bandwidths by default. Categorical
#' management is one-hot encoded within its own block. Every block is aligned,
#' quality controlled, imputed explicitly, and normalised before combination.
#'
#' @param weather,soil,management,geography Optional environment-by-variable
#'   matrices/data frames. Supply row names or an `environment` column.
#'   `geography` needs latitude and longitude.
#' @param environments Optional required environment order. By default the union
#'   of environment names across supplied blocks is used.
#' @param kernels Optional named character vector choosing `"gaussian"` or
#'   `"correlation"` per non-geographic block.
#' @param bandwidth_multipliers Positive Gaussian bandwidth multipliers around
#'   the median pairwise distance. Averaging several bandwidths reduces
#'   sensitivity in small environment networks.
#' @param min_coverage Minimum source coverage passed to
#'   [qc_environmental_data()].
#' @param missing_action,impute Missing-data policy passed to
#'   [qc_environmental_data()].
#' @param include_interactions Logical. When `TRUE`, add all available pairwise
#'   weather-by-soil, weather-by-management, and soil-by-management kernels.
#' @param interaction_terms Optional named list defining an exact set of
#'   interactions, for example
#'   `list(weather_x_soil_x_management = c("weather", "soil", "management"))`.
#'   Supplying this argument also enables interaction construction.
#' @param interaction_mode `"anova"` (default) constructs functional-ANOVA
#'   interactions from separately centred parent kernels, reducing main-effect
#'   leakage. `"product"` retains the legacy raw Hadamard product.
#' @param variable_interactions Optional dedicated variable-level hypotheses
#'   passed to [build_variable_interaction_kernels()]. Use these for a short,
#'   pre-specified set of interactions that must be represented and audited
#'   separately, such as temperature by root-zone bulk density. The broad
#'   modality interactions controlled by `include_interactions` remain distinct
#'   omnibus uncertainty structures.
#' @param variable_interaction_control Optional named list of advanced arguments
#'   passed to [build_variable_interaction_kernels()], excluding `covariates`
#'   and `interactions`. This keeps the main interface compact while allowing
#'   control of bandwidths, residualisation, hierarchy, and tensor dimension.
#'
#' @return A list with `kernels`, aligned `covariates`, `audit`, `provenance`,
#'   `imputation`, realised `bandwidths`, and dedicated variable-interaction
#'   diagnostics when requested.
#' @export
build_environment_kernels <- function(
    weather = NULL, soil = NULL, management = NULL, geography = NULL,
    environments = NULL, kernels = NULL,
    bandwidth_multipliers = c(0.5, 1, 2),
    min_coverage = 0.80,
    missing_action = c("impute", "warn", "error"),
    impute = c("median", "none"),
    include_interactions = FALSE,
    interaction_terms = NULL,
    interaction_mode = c("anova", "product"),
    variable_interactions = NULL,
    variable_interaction_control = list()) {
  missing_action <- match.arg(missing_action)
  impute <- match.arg(impute)
  interaction_mode <- match.arg(interaction_mode)
  if (!is.logical(include_interactions) || length(include_interactions) != 1L ||
      is.na(include_interactions))
    stop("`include_interactions` must be TRUE or FALSE.")
  if (!is.list(variable_interaction_control))
    stop("`variable_interaction_control` must be a named list.")
  if (length(variable_interaction_control) &&
      (is.null(names(variable_interaction_control)) ||
       anyNA(names(variable_interaction_control)) ||
       any(!nzchar(names(variable_interaction_control))) ||
       anyDuplicated(names(variable_interaction_control))))
    stop("`variable_interaction_control` must have unique, non-empty names.")
  protected_control <- intersect(
    names(variable_interaction_control), c("covariates", "interactions")
  )
  if (length(protected_control))
    stop("Do not supply ", paste(protected_control, collapse = ", "),
         " in `variable_interaction_control`.")
  blocks0 <- Filter(Negate(is.null), list(
    weather = weather, soil = soil, management = management,
    geography = geography
  ))
  if (!length(blocks0))
    stop("Supply at least one environmental block.")
  if (any(!is.finite(bandwidth_multipliers)) ||
      any(bandwidth_multipliers <= 0))
    stop("`bandwidth_multipliers` must be finite and positive.")
  bandwidth_multipliers <- sort(unique(as.numeric(bandwidth_multipliers)))

  prepared <- lapply(names(blocks0), function(nm)
    .prepare_environment_block(blocks0[[nm]], nm))
  names(prepared) <- names(blocks0)
  discovered <- unique(unlist(lapply(prepared, function(z)
    as.character(z$environment)), use.names = FALSE))
  if (is.null(environments)) environments <- discovered
  environments <- as.character(environments)
  if (!length(environments) || anyNA(environments) || any(!nzchar(environments)) ||
      anyDuplicated(environments))
    stop("`environments` must be unique, non-missing, non-empty identifiers.")
  outside <- setdiff(discovered, environments)
  if (length(outside))
    stop("Environmental blocks contain unknown environments: ",
         paste(outside, collapse = ", "), ".")

  default_kernels <- c(weather = "gaussian", soil = "gaussian",
                       management = "gaussian")
  if (is.null(kernels)) kernels <- default_kernels
  if (is.null(names(kernels)))
    stop("`kernels` must be a named vector by modality.")

  qcs <- list()
  Xs <- list()
  Ks <- list()
  bws <- list()
  for (nm in names(prepared)) {
    z <- prepared[[nm]]
    if (nm == "management") {
      mz <- .onehot(z[, setdiff(names(z), "environment"), drop = FALSE])
      z <- data.frame(environment = z$environment, mz,
                      check.names = FALSE, stringsAsFactors = FALSE)
    }
    absent <- setdiff(environments, z$environment)
    if (length(absent)) {
      filler <- z[rep(1L, length(absent)), , drop = FALSE]
      filler[] <- NA
      filler$environment <- absent
      z <- rbind(z, filler)
    }
    z <- z[match(environments, z$environment), , drop = FALSE]

    if (nm == "geography") {
      need <- match(c("latitude", "longitude"), tolower(names(z)))
      if (anyNA(need))
        stop("`geography` needs latitude and longitude columns.")
      q <- qc_environmental_data(
        z[, c("environment", names(z)[need]), drop = FALSE],
        environments = environments, min_coverage = min_coverage,
        missing_action = missing_action, impute = impute,
        add_missing_indicators = FALSE, source = nm
      )
      X <- data.matrix(q$data[, -1L, drop = FALSE])
      rownames(X) <- q$data$environment
      K <- build_environment_relationship(
        X, source = "enviromic", kernel = "gaussian", geo = TRUE
      )
      Ks[[nm]] <- .normalise_environment_kernel(K)
      Xs[[nm]] <- X
      qcs[[nm]] <- q
      bws[[nm]] <- attr(K, "bandwidth")
      next
    }

    q <- qc_environmental_data(
      z, environments = environments, min_coverage = min_coverage,
      missing_action = missing_action, impute = impute,
      add_missing_indicators = TRUE, source = nm
    )
    Xdf <- q$data[, setdiff(names(q$data), "environment"), drop = FALSE]
    Xdf <- .onehot(Xdf)
    X <- data.matrix(Xdf)
    rownames(X) <- q$data$environment
    informative <- vapply(seq_len(ncol(X)), function(j) {
      x <- X[, j]
      sum(is.finite(x)) > 1L && stats::sd(x, na.rm = TRUE) > 0
    }, logical(1))
    if (!any(informative)) {
      warning("Block '", nm, "' has no variable information; omitting it.",
              call. = FALSE)
      qcs[[nm]] <- q
      Xs[[nm]] <- X
      next
    }
    X <- X[, informative, drop = FALSE]
    which_kernel <- if (nm %in% names(kernels)) kernels[[nm]] else "gaussian"
    if (!which_kernel %in% c("gaussian", "correlation"))
      stop("Kernel for '", nm, "' must be gaussian or correlation.")
    if (which_kernel == "gaussian") {
      Z <- scale(X)
      if (any(!is.finite(Z)))
        stop("Internal non-finite values remained after QC in block '",
             nm, "'.")
      dd <- as.matrix(stats::dist(Z))
      base_bw <- stats::median(dd[upper.tri(dd) & dd > 0], na.rm = TRUE)
      if (!is.finite(base_bw) || base_bw <= 0) base_bw <- 1
      hs <- base_bw * bandwidth_multipliers
      Kset <- lapply(hs, function(h)
        build_environment_relationship(
          X, source = "enviromic", kernel = "gaussian", bandwidth = h
        ))
      K <- Reduce(`+`, Kset) / length(Kset)
      bws[[nm]] <- hs
    } else {
      K <- build_environment_relationship(
        X, source = "enviromic", kernel = "correlation"
      )
      bws[[nm]] <- NA_real_
    }
    Ks[[nm]] <- .normalise_environment_kernel(K)
    Xs[[nm]] <- X
    qcs[[nm]] <- q
  }
  if (!length(Ks))
    stop("No environmental block retained usable information.")

  if (include_interactions || length(interaction_terms))
    Ks <- add_environment_kernel_interactions(
      Ks, interactions = interaction_terms, mode = interaction_mode
    )
  variable_result <- NULL
  if (!is.null(variable_interactions)) {
    variable_result <- do.call(
      build_variable_interaction_kernels,
      c(
        list(covariates = Xs, interactions = variable_interactions),
        variable_interaction_control
      )
    )
    collisions <- intersect(names(Ks), names(variable_result$kernels))
    if (length(collisions))
      stop("Dedicated variable kernels duplicate existing kernel names: ",
           paste(collisions, collapse = ", "), ".")
    Ks <- c(Ks, variable_result$kernels)
  }
  block_diagnostics <- do.call(rbind, lapply(names(Xs), function(nm) {
    X <- Xs[[nm]]
    data.frame(
      block = nm, n_environments = nrow(X), n_variables = ncol(X),
      effective_rank = qr(scale(X, center = TRUE, scale = FALSE))$rank,
      stringsAsFactors = FALSE
    )
  }))
  kernel_agreement <- .environment_kernel_agreement(Ks)
  list(
    kernels = Ks,
    covariates = Xs,
    audit = .bind_qc_component(qcs, "audit"),
    provenance = .bind_qc_component(qcs, "provenance"),
    imputation = .bind_qc_component(qcs, "imputation"),
    bandwidths = bws,
    block_diagnostics = block_diagnostics,
    kernel_agreement = kernel_agreement,
    variable_interaction_audit = if (is.null(variable_result))
      NULL else variable_result$audit,
    variable_ledger = if (is.null(variable_result))
      NULL else variable_result$variable_ledger,
    kernel_hierarchy = if (is.null(variable_result))
      NULL else variable_result$hierarchy,
    variable_interaction_definitions = if (is.null(variable_result))
      NULL else variable_result$definitions,
    environments = environments
  )
}

#' Add functional-ANOVA interaction kernels
#'
#' Constructs separate positive-semidefinite interaction kernels without
#' combining or weighting them. In the default `"anova"` mode, each parent
#' kernel is centred as \eqn{H K H}, and the interaction is their Hadamard
#' product. This is the finite-sample tensor-product/functional-ANOVA
#' construction: it represents joint deviations after removing the constant
#' baseline from each main-effect feature space. It reduces, but cannot
#' eliminate, dependence between main-effect and interaction hypotheses.
#'
#' When `interactions = NULL`, every available pair among weather, soil, and
#' management is generated. Explicit named definitions may include higher-order
#' hypotheses. Geography is never included automatically but can be named
#' explicitly.
#'
#' @param kernels Named list of aligned positive-semidefinite kernels.
#' @param interactions Optional named list whose elements contain at least two
#'   parent-kernel names. `NULL` generates all available pairwise interactions
#'   among weather, soil, and management.
#' @param mode `"anova"` for centred functional-ANOVA products or `"product"`
#'   for legacy raw Hadamard products.
#'
#' @return The original named kernel list augmented with separately normalised
#'   interaction kernels. Each interaction records its parents and mode as
#'   attributes.
#' @export
add_environment_kernel_interactions <- function(
    kernels, interactions = NULL, mode = c("anova", "product")) {
  mode <- match.arg(mode)
  K <- .validate_environment_kernels(kernels)

  if (is.null(interactions)) {
    available <- intersect(c("weather", "soil", "management"), names(K))
    if (length(available) < 2L) return(K)
    cmb <- utils::combn(available, 2L, simplify = FALSE)
    interactions <- stats::setNames(
      cmb, vapply(cmb, function(z) paste(z, collapse = "_x_"),
                  character(1))
    )
  }
  if (!is.list(interactions) || !length(interactions))
    stop("`interactions` must be NULL or a non-empty list.")

  interaction_names <- names(interactions)
  if (is.null(interaction_names)) {
    interaction_names <- vapply(interactions, function(z)
      paste(as.character(z), collapse = "_x_"), character(1))
  } else {
    blank <- is.na(interaction_names) | !nzchar(interaction_names)
    interaction_names[blank] <- vapply(interactions[blank], function(z)
      paste(as.character(z), collapse = "_x_"), character(1))
  }
  if (anyDuplicated(interaction_names))
    stop("Interaction names must be unique.")
  if (any(interaction_names %in% names(K)))
    stop("Interaction names must not overwrite existing kernels: ",
         paste(intersect(interaction_names, names(K)), collapse = ", "), ".")

  n <- nrow(K[[1L]])
  H <- diag(n) - matrix(1 / n, n, n)
  for (i in seq_along(interactions)) {
    parts <- as.character(interactions[[i]])
    if (length(parts) < 2L || anyNA(parts) || any(!nzchar(parts)) ||
        anyDuplicated(parts) || any(!parts %in% names(K)))
      stop("Every interaction must name at least two distinct existing kernels.")
    parents <- K[parts]
    if (mode == "anova") {
      parents <- lapply(parents, function(M) {
        Mc <- H %*% M %*% H
        dimnames(Mc) <- dimnames(M)
        (Mc + t(Mc)) / 2
      })
    }
    Ki <- Reduce(`*`, parents)
    Ki <- (Ki + t(Ki)) / 2
    scale_ref <- max(1, max(abs(Ki)))
    if (max(abs(Ki)) <= 1e-12 ||
        any(diag(Ki) <= 1e-12 * scale_ref))
      stop("Interaction `", interaction_names[i],
           "` is degenerate after ", mode,
           " construction; revise or remove its parent kernels.")
    Ki <- .normalise_environment_kernel(.bend_pd(Ki))
    attr(Ki, "interaction_parents") <- parts
    attr(Ki, "interaction_mode") <- mode
    K[[interaction_names[i]]] <- Ki
  }
  K
}

#' Combine modality-specific environmental kernels
#'
#' Forms one positive-semidefinite environmental relationship matrix from
#' aligned modality-specific kernels. When more than one kernel is supplied,
#' the function requires explicit, non-negative weights; it does not interpret
#' equal weights as a scientific default. Users without empirical weights
#' should retain the kernels as separate uncertainty scenarios or use
#' [consensus_environment_kernels()] for descriptive integration.
#'
#' @param kernels Named list of aligned square kernels.
#' @param weights Non-negative named weights, required when more than one kernel
#'   is supplied. They are normalised to sum to one. Use
#'   [consensus_environment_kernels()] for descriptive integration without
#'   empirically justified weights.
#' @param identity_weight Optional shrinkage weight on an identity kernel.
#' @param interactions Optional named list of products such as
#'   `list(weather_management = c("weather", "management"))`.
#' @param interaction_mode Interaction construction passed to
#'   [add_environment_kernel_interactions()].
#'
#' @return A normalised positive-semidefinite environment relationship matrix
#'   carrying the realised weights as an attribute.
#' @export
combine_environment_kernels <- function(kernels, weights = NULL,
                                        identity_weight = 0,
                                        interactions = NULL,
                                        interaction_mode = c("anova",
                                                             "product")) {
  interaction_mode <- match.arg(interaction_mode)
  K <- .validate_environment_kernels(kernels)
  if (!is.numeric(identity_weight) || length(identity_weight) != 1L ||
      !is.finite(identity_weight) || identity_weight < 0 ||
      identity_weight >= 1)
    stop("`identity_weight` must be one finite number in [0, 1).")
  if (length(interactions))
    K <- add_environment_kernel_interactions(
      K, interactions = interactions, mode = interaction_mode
    )
  if (is.null(weights)) {
    if (length(K) > 1L)
      stop("`weights` is required for a multi-kernel weighted combination. ",
           "Use `consensus_environment_kernels()` for weight-free descriptive ",
           "integration.")
    weights <- stats::setNames(1, names(K))
  } else {
    if (is.null(names(weights)))
      stop("`weights` must be named by kernel.")
    miss <- setdiff(names(K), names(weights))
    extra <- setdiff(names(weights), names(K))
    if (length(miss))
      stop("`weights` is missing kernels: ", paste(miss, collapse = ", "), ".")
    if (length(extra))
      stop("`weights` names unknown kernels: ", paste(extra, collapse = ", "), ".")
    weights <- weights[names(K)]
    if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0)
      stop("Kernel weights must be finite, non-negative, and not all zero.")
    weights <- weights / sum(weights)
  }
  weights <- weights * (1 - identity_weight)
  out <- Reduce(`+`, Map(`*`, K, as.list(weights)))
  if (identity_weight > 0)
    out <- out + diag(identity_weight, nrow(out))
  out <- .normalise_environment_kernel(.bend_pd(out))
  attr(out, "weights") <- c(weights, identity = identity_weight)
  out
}

#' Weight-free consensus of environmental kernels
#'
#' @description
#' Constructs a descriptive, multimodal environmental similarity without
#' eliciting or estimating modality weights. Each kernel is converted to
#' within-kernel off-diagonal ranks, the ranks are combined entrywise by their
#' median, and the result is projected to a valid correlation-like matrix. This
#' robust consensus is intended for descriptive clustering and environmental
#' representativeness only. It is not an estimate of genetic covariance.
#'
#' @param kernels Named list of aligned square environmental kernels.
#'
#' @return A positive-semidefinite, unit-diagonal descriptive environmental
#'   similarity matrix. The `integration` attribute records the weight-free
#'   rule used.
#' @export
consensus_environment_kernels <- function(kernels) {
  K <- .validate_environment_kernels(kernels)
  n <- nrow(K[[1L]])
  envs <- rownames(K[[1L]])
  if (length(K) == 1L) {
    out <- K[[1L]]
    attr(out, "integration") <- "single_kernel_no_weights"
    return(out)
  }

  upper <- upper.tri(K[[1L]])
  ranked <- lapply(K, function(M) {
    z <- M[upper]
    if (length(unique(z)) <= 1L) {
      r <- rep(0.5, length(z))
    } else {
      r <- (rank(z, ties.method = "average") - 1) / (length(z) - 1)
    }
    R <- matrix(0, n, n, dimnames = list(envs, envs))
    R[upper] <- r
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
    diag(R) <- 1
    R
  })
  offdiag <- do.call(cbind, lapply(ranked, function(M) M[upper]))
  med <- apply(offdiag, 1L, stats::median)
  out <- matrix(0, n, n, dimnames = list(envs, envs))
  out[upper] <- med
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  diag(out) <- 1
  out <- .normalise_environment_kernel(.bend_pd(out))
  attr(out, "integration") <- "entrywise_median_of_within_kernel_ranks"
  out
}

#' Calibrate environmental covariance from historical MET information
#'
#' @description
#' Learns non-negative, sum-to-one weights for modality kernels against an
#' empirical genetic-correlation target. Calibration can use a target matrix, a
#' list of year-specific target matrices, or historical genotype-by-environment
#' values. Leave-environment-out predictions and bootstrap weight uncertainty
#' are returned. Without a historical target, weights are not defined:
#' `Sigma_E` is the identity matrix and all main-effect and interaction kernels
#' are returned as separate, unweighted sensitivity scenarios. Lack of
#' response evidence is not treated as absence of a physical interaction.
#' Environmental similarity is never converted into an assumed genetic
#' covariance by arbitrary prior weights.
#'
#' @param kernels Named list of aligned environmental kernels.
#' @param target Optional empirical environment genetic-correlation matrix or
#'   list of such matrices.
#' @param historical Optional long data frame containing genotype, environment,
#'   and genetic-value columns.
#' @param genotype_col,environment_col,value_col Column names in `historical`.
#' @param year_col Optional year column. When supplied, a target is built per
#'   year and leave-year-out validation is added.
#' @param prior_weights Optional named modality weights used only as a weak
#'   ridge target when a historical target is supplied. They are ignored when
#'   `target` and `historical` are both `NULL`.
#' @param identity_weight Prior identity/shrinkage weight in `[0,1)`, used only
#'   to construct the ridge target when historical MET evidence is supplied.
#' @param ridge Non-negative penalty toward `prior_weights`. Set to zero for
#'   wholly data-calibrated weights.
#' @param interaction_kernels Optional character vector naming interaction
#'   kernels. `NULL` identifies names containing `"_x_"`; `character(0)`
#'   declares that none are interactions.
#' @param eligible_interactions Optional character vector restricting which
#'   interaction kernels may enter the central covariance. Kernels excluded by
#'   this screen remain in `candidates` as structural uncertainty scenarios.
#'   Supply the `supported` component from
#'   [assess_variable_interactions()] after component-wise blocked validation.
#'   `NULL` makes every identified interaction eligible for the aggregate
#'   evidence gate.
#' @param interaction_policy `"evidence"` (default) activates interaction
#'   components only when a paired blocked-validation sign test favours the
#'   main-plus-interaction model over the main-only model. `"exclude"` fixes all
#'   interaction weights at zero in the central covariance. `"include"` is an
#'   explicit user override for the central covariance.
#'   Leave-year/target-out folds are preferred when available; otherwise
#'   leave-environment-out folds are used.
#'   This policy never removes interaction kernels from the sensitivity set:
#'   absence of response evidence is not evidence that physical interactions
#'   are absent. Without a historical target, interaction kernels remain
#'   separate, unweighted uncertainty scenarios and never enter the central
#'   covariance.
#' @param interaction_alpha Significance level for the one-sided paired sign
#'   test used by `interaction_policy = "evidence"`.
#' @param n_boot Number of pair-bootstrap calibration draws.
#' @param seed Optional seed for bootstrap reproducibility.
#'
#' @return A list with `Sigma_E`, weights (or `NULL` without historical MET
#'   evidence), status, interaction status/evidence, fitted target, diagnostics,
#'   blocked validation, bootstrap draws, and sensitivity candidates.
#' @export
calibrate_environment_covariance <- function(
    kernels, target = NULL, historical = NULL,
    genotype_col = "genotype", environment_col = "environment",
    value_col = "value", year_col = NULL,
    prior_weights = NULL, identity_weight = 0.10,
    ridge = 0, interaction_kernels = NULL,
    eligible_interactions = NULL,
    interaction_policy = c("evidence", "exclude", "include"),
    interaction_alpha = 0.05,
    n_boot = 200L, seed = NULL) {
  interaction_policy <- match.arg(interaction_policy)
  interaction_attributes <- names(kernels)[vapply(kernels, function(M)
    identical(attr(M, "kernel_role"), "interaction"), logical(1))]
  K <- .validate_environment_kernels(kernels)
  envs <- rownames(K[[1L]])
  if (is.null(interaction_kernels)) {
    interaction_kernels <- union(
      names(K)[grepl("_x_", names(K), fixed = TRUE)],
      interaction_attributes
    )
  } else {
    interaction_kernels <- as.character(interaction_kernels)
    if (anyNA(interaction_kernels) || any(!nzchar(interaction_kernels)) ||
        anyDuplicated(interaction_kernels) ||
        any(!interaction_kernels %in% names(K)))
      stop("`interaction_kernels` must uniquely name kernels in `kernels`.")
  }
  main_kernels <- setdiff(names(K), interaction_kernels)
  if (!length(main_kernels))
    stop("At least one non-interaction kernel is required.")
  if (is.null(eligible_interactions)) {
    eligible_interactions <- interaction_kernels
  } else {
    eligible_interactions <- as.character(eligible_interactions)
    if (anyNA(eligible_interactions) ||
        any(!nzchar(eligible_interactions)) ||
        anyDuplicated(eligible_interactions) ||
        any(!eligible_interactions %in% interaction_kernels))
      stop("`eligible_interactions` must uniquely name identified ",
           "interaction kernels.")
  }
  if (!is.null(historical)) {
    if (!is.null(target))
      stop("Supply only one of `target` and `historical`.")
    if (!is.null(year_col)) {
      if (!year_col %in% names(historical))
        stop("`historical` has no `", year_col, "` column.")
      hs <- split(historical, historical[[year_col]])
      target <- lapply(hs, .historical_environment_correlation,
                       envs = envs, genotype_col = genotype_col,
                       environment_col = environment_col,
                       value_col = value_col)
    } else {
      target <- .historical_environment_correlation(
        historical, envs, genotype_col, environment_col, value_col
      )
    }
  }
  if (!is.numeric(ridge) || length(ridge) != 1L ||
      !is.finite(ridge) || ridge < 0)
    stop("`ridge` must be a finite non-negative scalar.")
  if (!is.numeric(identity_weight) || length(identity_weight) != 1L ||
      !is.finite(identity_weight) || identity_weight < 0 ||
      identity_weight >= 1)
    stop("`identity_weight` must be one finite number in [0, 1).")
  if (!is.numeric(interaction_alpha) || length(interaction_alpha) != 1L ||
      !is.finite(interaction_alpha) || interaction_alpha <= 0 ||
      interaction_alpha >= 1)
    stop("`interaction_alpha` must be one finite number in (0, 1).")
  if (!is.numeric(n_boot) || length(n_boot) != 1L ||
      !is.finite(n_boot) || n_boot < 0)
    stop("`n_boot` must be a non-negative integer.")
  n_boot <- as.integer(n_boot)

  kernel_names <- c(names(K), "identity")
  identity_kernel <- diag(length(envs))
  dimnames(identity_kernel) <- list(envs, envs)
  Kall <- c(K, list(identity = identity_kernel))

  if (is.null(target)) {
    candidates <- c(
      list(independent = Kall$identity),
      stats::setNames(
        lapply(K, .normalise_environment_kernel),
        paste0("kernel_", names(K))
      )
    )
    n_interactions <- length(interaction_kernels)
    interaction_evidence <- data.frame(
      interaction = interaction_kernels,
      policy = rep(interaction_policy, n_interactions),
      model_activated = rep(FALSE, n_interactions),
      component_active = rep(FALSE, n_interactions),
      final_weight = rep(NA_real_, n_interactions),
      bootstrap_inclusion = rep(NA_real_, n_interactions),
      cv_mean_rmse_improvement = rep(NA_real_, n_interactions),
      cv_positive_folds = rep(NA_integer_, n_interactions),
      cv_non_tied_folds = rep(NA_integer_, n_interactions),
      sign_test_p = rep(NA_real_, n_interactions),
      validation_basis = rep(NA_character_, n_interactions),
      reason = rep(
        paste0(
          "Retained as an unweighted structural uncertainty scenario; no ",
          "historical MET target was available for central weighting."
        ),
        n_interactions
      ),
      stringsAsFactors = FALSE
    )
    return(list(
      Sigma_E = Kall$identity, weights = NULL,
      status = "no_historical_met",
      interaction_status = if (!length(interaction_kernels))
        "no_interaction_kernels" else "unweighted_interaction_uncertainty",
      interaction_evidence = interaction_evidence,
      target = NULL, fitted = NULL,
      diagnostics = data.frame(
        metric = c("historical_target_available",
                   "kernel_weights_estimable",
                   "central_covariance_is_identity",
                   "calibration_rmse"),
        value = c(0, 0, 1, NA_real_), stringsAsFactors = FALSE
      ),
      cross_validation = NULL, bootstrap_weights = NULL,
      candidates = candidates
    ))
  }

  if (is.null(prior_weights)) {
    prior_weights <- c(
      stats::setNames(rep((1 - identity_weight) / length(K), length(K)),
                      names(K)),
      identity = identity_weight
    )
  } else {
    if (is.null(names(prior_weights)) ||
        any(!kernel_names %in% names(prior_weights)))
      stop("`prior_weights` must name every kernel and optionally identity.")
    prior_weights <- prior_weights[kernel_names]
    if (any(!is.finite(prior_weights)) || any(prior_weights < 0) ||
        sum(prior_weights) <= 0)
      stop("`prior_weights` must be finite, non-negative, and not all zero.")
    prior_weights <- prior_weights / sum(prior_weights)
  }

  targets <- if (is.list(target) && !is.matrix(target)) target else list(target)
  targets <- lapply(targets, function(T)
    .align_environment_target(T, envs))
  target_mean <- .mean_environment_targets(targets)
  idx <- upper.tri(target_mean) & is.finite(target_mean)
  main_names <- c(main_kernels, "identity")
  p_main <- length(main_names)
  if (sum(idx) < p_main)
    stop("Historical target has too few finite environment pairs to calibrate ",
         p_main, " main-effect weights.")
  X <- do.call(cbind, lapply(Kall, function(M) M[idx]))
  y <- target_mean[idx]
  prior_main <- prior_weights[main_names]
  prior_main <- prior_main / sum(prior_main)
  Xmain <- X[, main_names, drop = FALSE]
  fit_main <- .fit_kernel_simplex(Xmain, y, prior_main, ridge)
  names(fit_main) <- main_names

  full_names <- c(main_names, eligible_interactions)
  prior_full <- prior_weights[full_names]
  if (sum(prior_full) <= 0)
    stop("`prior_weights` assigns no mass to the main and eligible ",
         "interaction kernels.")
  prior_full <- prior_full / sum(prior_full)
  full_estimable <- sum(idx) >= length(full_names)
  fit_full <- NULL
  if (length(eligible_interactions) && full_estimable &&
      !identical(interaction_policy, "exclude")) {
    fit_full <- .fit_kernel_simplex(
      X[, full_names, drop = FALSE], y, prior_full, ridge
    )
    names(fit_full) <- full_names
  } else if (identical(interaction_policy, "include") &&
             length(eligible_interactions)) {
    stop("Historical target has too few finite environment pairs to fit the ",
         "explicitly included interaction model.")
  }

  cv_main <- .cross_validate_environment_kernels(
    Kall[main_names], target_mean, targets, envs,
    prior_main, ridge, model = "main_only"
  )
  cv_full <- if (!is.null(fit_full))
    .cross_validate_environment_kernels(
      Kall[full_names], target_mean, targets, envs,
      prior_full, ridge, model = "main_plus_interactions"
    ) else NULL
  cv <- if (is.null(cv_main)) cv_full else if (is.null(cv_full))
    cv_main else rbind(cv_main, cv_full)

  cv_improvement <- NA_real_
  sign_p <- NA_real_
  positive_folds <- 0L
  non_tied_folds <- 0L
  validation_basis <- NA_character_
  interaction_accepted <- FALSE
  evidence_reason <- "No interaction kernels were supplied."
  if (length(interaction_kernels) && !length(eligible_interactions)) {
    evidence_reason <- paste0(
      "All interactions were screened out of the central covariance; their ",
      "kernels remain structural uncertainty scenarios."
    )
  } else if (length(eligible_interactions)) {
    if (identical(interaction_policy, "exclude")) {
      evidence_reason <- paste0(
        "Interactions were excluded from the central covariance by policy; ",
        "their kernels remain structural uncertainty scenarios."
      )
    } else if (identical(interaction_policy, "include")) {
      interaction_accepted <- !is.null(fit_full)
      evidence_reason <- "Interactions were included by explicit user policy."
    } else if (is.null(fit_full)) {
      evidence_reason <- paste0(
        "The full interaction model was not estimable from the available ",
        "historical environment pairs; interaction kernels remain structural ",
        "uncertainty scenarios."
      )
    } else if (is.null(cv_main) || is.null(cv_full)) {
      evidence_reason <- paste0(
        "Blocked validation could not compare the main-only and interaction ",
        "models; interaction kernels remain structural uncertainty scenarios."
      )
    } else {
      validation_basis <- if (
        any(cv_main$validation == "leave_target_out") &&
        any(cv_full$validation == "leave_target_out")
      ) "leave_target_out" else "leave_environment_out"
      paired <- merge(
        cv_main[cv_main$validation == validation_basis,
                c("validation", "held_out", "rmse")],
        cv_full[cv_full$validation == validation_basis,
                c("validation", "held_out", "rmse")],
        by = c("validation", "held_out"),
        suffixes = c("_main", "_interaction")
      )
      improvement <- paired$rmse_main - paired$rmse_interaction
      improvement <- improvement[is.finite(improvement)]
      cv_improvement <- if (length(improvement))
        mean(improvement) else NA_real_
      tol <- sqrt(.Machine$double.eps)
      signs <- improvement[abs(improvement) > tol]
      non_tied_folds <- length(signs)
      positive_folds <- sum(signs > 0)
      sign_p <- if (non_tied_folds)
        stats::pbinom(positive_folds - 1L, non_tied_folds, 0.5,
                      lower.tail = FALSE) else NA_real_
      positive_weight <- sum(fit_full[eligible_interactions]) > tol
      interaction_accepted <- is.finite(cv_improvement) &&
        cv_improvement > 0 && is.finite(sign_p) &&
        sign_p < interaction_alpha && positive_weight
      evidence_reason <- if (interaction_accepted) {
        "Blocked validation significantly favoured main plus interactions."
      } else {
        paste0(
          "Interactions did not pass the paired blocked-validation evidence ",
          "gate; their central weights were fixed at zero, but their kernels ",
          "remain structural uncertainty scenarios."
        )
      }
    }
  }

  active_names <- if (interaction_accepted) full_names else main_names
  fit_active <- if (interaction_accepted) fit_full else fit_main
  fit <- stats::setNames(rep(0, length(kernel_names)), kernel_names)
  fit[active_names] <- fit_active
  Sigma <- Reduce(`+`, Map(`*`, Kall, as.list(fit)))
  Sigma <- .normalise_environment_kernel(.bend_pd(Sigma))
  fitted <- matrix(NA_real_, length(envs), length(envs),
                   dimnames = list(envs, envs))
  fitted[idx] <- as.numeric(X %*% fit)
  fitted[lower.tri(fitted)] <- t(fitted)[lower.tri(fitted)]
  diag(fitted) <- 1

  boot <- NULL
  if (n_boot > 0L) {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }, add = TRUE)
      set.seed(seed)
    }
    boot_active <- t(replicate(n_boot, {
      ii <- sample.int(length(y), length(y), replace = TRUE)
      .fit_kernel_simplex(
        X[ii, active_names, drop = FALSE], y[ii],
        if (interaction_accepted) prior_full else prior_main, ridge
      )
    }))
    boot <- matrix(0, n_boot, length(kernel_names),
                   dimnames = list(NULL, kernel_names))
    boot[, active_names] <- boot_active
    colnames(boot) <- kernel_names
  }
  rmse <- sqrt(mean((y - as.numeric(X %*% fit))^2))
  Xactive <- X[, active_names, drop = FALSE]
  design_rank <- qr(Xactive)$rank
  design_condition <- tryCatch(kappa(Xactive), error = function(e) Inf)
  diagnostics <- data.frame(
    metric = c("historical_target_available", "n_target_matrices",
               "n_environment_pairs", "calibration_rmse",
               "kernel_design_rank", "kernel_weight_count",
               "weights_fully_identifiable", "kernel_design_condition",
               "interaction_kernel_count", "eligible_interaction_count",
               "interaction_model_activated",
               "interaction_cv_mean_rmse_improvement",
               "interaction_cv_positive_folds",
               "interaction_cv_non_tied_folds",
               "interaction_sign_test_p"),
    value = c(1, length(targets), length(y), rmse,
              design_rank, length(active_names),
              as.numeric(design_rank == length(active_names)),
              design_condition, length(interaction_kernels),
              length(eligible_interactions),
              as.numeric(interaction_accepted), cv_improvement,
              positive_folds, non_tied_folds, sign_p),
    stringsAsFactors = FALSE
  )
  bootstrap_inclusion <- if (!is.null(boot) && length(interaction_kernels))
    colMeans(boot[, interaction_kernels, drop = FALSE] >
               sqrt(.Machine$double.eps)) else
    stats::setNames(rep(NA_real_, length(interaction_kernels)),
                    interaction_kernels)
  n_interactions <- length(interaction_kernels)
  eligible_flag <- interaction_kernels %in% eligible_interactions
  interaction_evidence <- data.frame(
    interaction = interaction_kernels,
    policy = rep(interaction_policy, n_interactions),
    model_activated = unname(interaction_accepted & eligible_flag),
    component_active = unname(
      interaction_accepted & eligible_flag &
        fit[interaction_kernels] > sqrt(.Machine$double.eps)
    ),
    final_weight = unname(fit[interaction_kernels]),
    bootstrap_inclusion = unname(bootstrap_inclusion[interaction_kernels]),
    cv_mean_rmse_improvement = ifelse(
      eligible_flag, cv_improvement, NA_real_
    ),
    cv_positive_folds = ifelse(
      eligible_flag, positive_folds, NA_integer_
    ),
    cv_non_tied_folds = ifelse(
      eligible_flag, non_tied_folds, NA_integer_
    ),
    sign_test_p = ifelse(eligible_flag, sign_p, NA_real_),
    validation_basis = ifelse(
      eligible_flag, validation_basis, NA_character_
    ),
    reason = ifelse(
      interaction_kernels %in% eligible_interactions,
      evidence_reason,
      paste0(
        "Excluded from the central covariance by component-wise screening; ",
        "retained as a structural uncertainty scenario."
      )
    ),
    stringsAsFactors = FALSE
  )
  candidates <- c(
    list(calibrated = Sigma, independent = Kall$identity),
    stats::setNames(
      lapply(K, .normalise_environment_kernel),
      paste0("kernel_", names(K))
    )
  )
  if (!is.null(boot)) {
    probs <- c(0.10, 0.50, 0.90)
    for (q in probs) {
      wq <- apply(boot, 2L, stats::quantile, probs = q, na.rm = TRUE)
      wq <- pmax(0, wq); wq <- wq / sum(wq)
      candidates[[paste0("bootstrap_q", formatC(100 * q, width = 2,
                                                flag = "0"))]] <-
        .normalise_environment_kernel(.bend_pd(
          Reduce(`+`, Map(`*`, Kall, as.list(wq)))
        ))
    }
  }
  list(
    Sigma_E = Sigma, weights = fit, status = "historically_calibrated",
    interaction_status = if (!length(interaction_kernels))
      "no_interaction_kernels" else if (!length(eligible_interactions))
        "interactions_screened_out" else if (interaction_accepted)
        "interaction_model_activated" else "interaction_model_rejected",
    interaction_evidence = interaction_evidence,
    target = target_mean, fitted = fitted, diagnostics = diagnostics,
    cross_validation = cv, bootstrap_weights = boot,
    candidates = candidates
  )
}

.cross_validate_environment_kernels <- function(
    Kall, target_mean, targets, envs, prior_weights, ridge, model) {
  p <- length(Kall)
  cv_environment <- do.call(rbind, lapply(envs, function(held) {
    train_env <- setdiff(envs, held)
    train_idx <- outer(envs, envs, function(a, b)
      a %in% train_env & b %in% train_env) &
      upper.tri(target_mean) & is.finite(target_mean)
    test_idx <- outer(envs, envs, function(a, b)
      xor(a == held, b == held)) &
      upper.tri(target_mean) & is.finite(target_mean)
    if (sum(train_idx) < p || !any(test_idx)) return(NULL)
    Xtr <- do.call(cbind, lapply(Kall, function(M) M[train_idx]))
    w <- .fit_kernel_simplex(Xtr, target_mean[train_idx],
                             prior_weights, ridge)
    Xte <- do.call(cbind, lapply(Kall, function(M) M[test_idx]))
    pred <- as.numeric(Xte %*% w)
    obs <- target_mean[test_idx]
    data.frame(
      model = model, validation = "leave_environment_out",
      held_out = held, n_pairs = length(obs),
      rmse = sqrt(mean((obs - pred)^2)),
      correlation = if (length(obs) > 2L)
        suppressWarnings(stats::cor(obs, pred)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  cv_target <- NULL
  if (length(targets) > 1L) {
    cv_target <- do.call(rbind, lapply(seq_along(targets), function(k) {
      train_target <- .mean_environment_targets(targets[-k])
      tr <- upper.tri(train_target) & is.finite(train_target)
      te <- upper.tri(targets[[k]]) & is.finite(targets[[k]])
      if (sum(tr) < p || !any(te)) return(NULL)
      Xtr <- do.call(cbind, lapply(Kall, function(M) M[tr]))
      wt <- .fit_kernel_simplex(Xtr, train_target[tr],
                                prior_weights, ridge)
      Xte <- do.call(cbind, lapply(Kall, function(M) M[te]))
      pred <- as.numeric(Xte %*% wt)
      obs <- targets[[k]][te]
      data.frame(
        model = model, validation = "leave_target_out",
        held_out = if (is.null(names(targets)))
          as.character(k) else names(targets)[k],
        n_pairs = length(obs),
        rmse = sqrt(mean((obs - pred)^2)),
        correlation = if (length(obs) > 2L)
          suppressWarnings(stats::cor(obs, pred)) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }
  if (is.null(cv_environment)) cv_target else if (is.null(cv_target))
    cv_environment else rbind(cv_environment, cv_target)
}

.prepare_environment_block <- function(x, source) {
  d <- as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"environment" %in% names(d)) {
    if (is.null(rownames(d)))
      stop("Block '", source, "' needs row names or an environment column.")
    d$environment <- rownames(d)
    d <- d[, c("environment", setdiff(names(d), "environment")), drop = FALSE]
  }
  d$environment <- as.character(d$environment)
  if (anyDuplicated(d$environment))
    stop("Block '", source, "' has duplicate environment rows.")
  d
}

.bind_qc_component <- function(x, component) {
  z <- lapply(x, `[[`, component)
  z <- z[!vapply(z, is.null, logical(1))]
  if (!length(z)) return(NULL)
  do.call(rbind, z)
}

.normalise_environment_kernel <- function(K, eps = 1e-10) {
  dn <- dimnames(K)
  K <- (as.matrix(K) + t(as.matrix(K))) / 2
  dd <- sqrt(pmax(diag(K), eps))
  K <- K / tcrossprod(dd)
  diag(K) <- 1
  dimnames(K) <- dn
  K
}

.validate_environment_kernels <- function(kernels) {
  if (!is.list(kernels) || !length(kernels) || is.null(names(kernels)) ||
      any(!nzchar(names(kernels))) || anyDuplicated(names(kernels)))
    stop("`kernels` must be a non-empty uniquely named list.")
  K <- lapply(kernels, as.matrix)
  ref <- rownames(K[[1L]])
  if (is.null(ref) || anyDuplicated(ref))
    stop("Kernels need unique environment row names.")
  for (nm in names(K)) {
    M <- K[[nm]]
    if (nrow(M) != ncol(M) || is.null(colnames(M)) ||
        !setequal(rownames(M), colnames(M)) ||
        !setequal(ref, rownames(M)))
      stop("Kernel '", nm, "' is not an aligned named square matrix.")
    M <- M[ref, ref, drop = FALSE]
    if (any(!is.finite(M)))
      stop("Kernel '", nm, "' contains non-finite values.")
    if (max(abs(M - t(M))) > 1e-8)
      stop("Kernel '", nm, "' is not symmetric.")
    ee <- eigen((M + t(M)) / 2, symmetric = TRUE,
                only.values = TRUE)$values
    if (min(ee) < -1e-7 * max(1, max(abs(ee))))
      stop("Kernel '", nm, "' is not positive semidefinite.")
    K[[nm]] <- .normalise_environment_kernel(.bend_pd(M))
  }
  K
}

.historical_environment_correlation <- function(d, envs, genotype_col,
                                                environment_col, value_col) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  need <- c(genotype_col, environment_col, value_col)
  if (!all(need %in% names(d)))
    stop("`historical` needs columns: ", paste(need, collapse = ", "), ".")
  d <- d[is.finite(as.numeric(d[[value_col]])), need, drop = FALSE]
  d[[value_col]] <- as.numeric(d[[value_col]])
  av <- stats::aggregate(
    d[[value_col]],
    list(genotype = as.character(d[[genotype_col]]),
         environment = as.character(d[[environment_col]])),
    mean
  )
  gen <- unique(av$genotype)
  M <- matrix(NA_real_, length(gen), length(envs),
              dimnames = list(gen, envs))
  ii <- match(av$genotype, gen)
  jj <- match(av$environment, envs)
  keep <- !is.na(jj)
  M[cbind(ii[keep], jj[keep])] <- av$x[keep]
  D <- suppressWarnings(stats::cor(M, use = "pairwise.complete.obs"))
  dimnames(D) <- list(envs, envs)
  diag(D) <- 1
  D
}

.align_environment_target <- function(T, envs) {
  T <- as.matrix(T)
  if (nrow(T) != ncol(T) || is.null(rownames(T)) || is.null(colnames(T)) ||
      !all(envs %in% rownames(T)) || !all(envs %in% colnames(T)))
    stop("Every target must be a named square matrix covering all kernels' ",
         "environments.")
  T <- T[envs, envs, drop = FALSE]
  T <- (T + t(T)) / 2
  diag(T) <- 1
  T
}

.mean_environment_targets <- function(targets) {
  arr <- simplify2array(targets)
  out <- apply(arr, c(1, 2), function(z) {
    z <- z[is.finite(z)]
    if (length(z)) mean(z) else NA_real_
  })
  dimnames(out) <- dimnames(targets[[1L]])
  out
}

.fit_kernel_simplex <- function(X, y, prior, ridge = 0,
                                maxit = 5000L, tol = 1e-10) {
  X <- as.matrix(X); y <- as.numeric(y); prior <- as.numeric(prior)
  if (ncol(X) != length(prior))
    stop("Internal kernel calibration dimension mismatch.")
  L <- 2 * sum(X^2) + 2 * ridge
  step <- if (is.finite(L) && L > 0) 1 / L else 0.01
  ## With no ridge penalty, the prior must not influence a data-calibrated fit.
  ## A uniform start is a deterministic, label-invariant tie breaker if kernels
  ## are collinear and the constrained optimum is not unique.
  w <- if (ridge > 0) .project_simplex(prior) else
    rep(1 / ncol(X), ncol(X))
  for (iter in seq_len(maxit)) {
    grad <- 2 * crossprod(X, X %*% w - y) + 2 * ridge * (w - prior)
    next_w <- .project_simplex(w - step * as.numeric(grad))
    if (max(abs(next_w - w)) < tol) {
      w <- next_w
      break
    }
    w <- next_w
  }
  w
}

.project_simplex <- function(v) {
  v <- as.numeric(v)
  u <- sort(v, decreasing = TRUE)
  css <- cumsum(u) - 1
  rho <- max(which(u - css / seq_along(u) > 0))
  theta <- css[rho] / rho
  pmax(v - theta, 0)
}

.environment_kernel_agreement <- function(kernels) {
  if (length(kernels) < 2L)
    return(data.frame(kernel_1 = character(), kernel_2 = character(),
                      off_diagonal_correlation = numeric(),
                      stringsAsFactors = FALSE))
  nm <- names(kernels)
  rows <- list()
  for (i in seq_len(length(kernels) - 1L)) {
    for (j in (i + 1L):length(kernels)) {
      a <- kernels[[i]][upper.tri(kernels[[i]])]
      b <- kernels[[j]][upper.tri(kernels[[j]])]
      rows[[length(rows) + 1L]] <- data.frame(
        kernel_1 = nm[i], kernel_2 = nm[j],
        off_diagonal_correlation =
          suppressWarnings(stats::cor(a, b, use = "pairwise.complete.obs")),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
