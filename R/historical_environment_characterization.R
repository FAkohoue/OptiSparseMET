#' Historical multimodal environmental characterization
#'
#' Runs the planning-stage weather, soil, kernel, consensus, descriptive-strata,
#' and strict-validation steps as one audited workflow. Weather years remain
#' separate evidence replicates, while soil remains a distinct evidence block.
#' Without a historical genetic-response target, multimodal descriptive
#' integration uses [consensus_environment_kernels()] and no arbitrary modality
#' weights are estimated.
#'
#' @param sites,years,window,daily_by_year,pars,envirotype Inputs passed to
#'   [historical_envirotype()]. Offline workflows may supply `daily_by_year` and
#'   leave `sites` unset when soil is supplied separately.
#' @param soil Optional already-fetched SoilGrids output or environment-by-soil
#'   feature matrix/data frame.
#' @param fetch_soil Logical; retrieve SoilGrids data from `sites`.
#' @param history_control Additional named arguments for
#'   [historical_envirotype()].
#' @param stability_control Additional named arguments for
#'   [assess_envirotype_stability()].
#' @param soil_control Additional named arguments for [fetch_soilgrids()].
#' @param soil_profile_control Additional named arguments for
#'   [soil_profile_features()]. Multi-depth SoilGrids columns are summarized
#'   automatically; already summarized soil features are used as supplied.
#' @param kernel_control Additional named arguments for
#'   [build_environment_kernels()], including `redundancy` and
#'   `redundancy_control`.
#' @param inference_control Additional named arguments for environmental
#'   partition inference. `n_boot` and thresholds can be set here.
#' @param covariance_control Additional named arguments for
#'   [calibrate_environment_covariance()].
#' @param historical_target Optional named environment relationship/covariance
#'   matrix, or list of matrices, derived from historical MET responses. When
#'   supplied it is passed to [calibrate_environment_covariance()].
#' @param seed Reproducibility seed used by partition inference.
#'
#' @return A list containing historical weather summaries and stability, raw and
#'   profile soil data, modality kernels and redundancy audits, descriptive
#'   consensus, weather-soil agreement, candidate environmental strata, strict
#'   mega-environment validation, covariance calibration, and provenance/QC.
#' @export
historical_environment_characterization <- function(
    sites = NULL, years = NULL, window = NULL, daily_by_year = NULL,
    pars = NULL, envirotype = list(), soil = NULL, fetch_soil = FALSE,
    history_control = list(), stability_control = list(), soil_control = list(),
    soil_profile_control = list(), kernel_control = list(),
    inference_control = list(), covariance_control = list(),
    historical_target = NULL, seed = 1L) {
  controls <- list(
    history_control = history_control, stability_control = stability_control,
    soil_control = soil_control, soil_profile_control = soil_profile_control,
    kernel_control = kernel_control, inference_control = inference_control,
    covariance_control = covariance_control
  )
  for (nm in names(controls)) {
    if (!is.list(controls[[nm]]) ||
        (length(controls[[nm]]) &&
         (is.null(names(controls[[nm]])) ||
          any(!nzchar(names(controls[[nm]]))) ||
          anyDuplicated(names(controls[[nm]])))))
      stop("`", nm, "` must be a uniquely named list.")
  }
  if (!is.logical(fetch_soil) || length(fetch_soil) != 1L || is.na(fetch_soil))
    stop("`fetch_soil` must be TRUE or FALSE.")

  .reject_characterization_overrides(
    history_control,
    c("sites", "years", "window", "daily_by_year", "pars", "envirotype")
  )
  history <- do.call(historical_envirotype, c(list(
    sites = sites, years = years, window = window,
    daily_by_year = daily_by_year, pars = pars, envirotype = envirotype
  ), history_control))
  if (is.null(history)) return(NULL)

  stability <- NULL
  shared_environments <- Reduce(intersect, lapply(history$per_year, rownames))
  if (length(history$per_year) >= 2L && length(shared_environments) >= 3L) {
    .reject_characterization_overrides(stability_control, "hist")
    stability <- do.call(
      assess_envirotype_stability,
      c(list(hist = history), stability_control)
    )
  }

  soil_raw <- soil
  if (is.null(soil_raw) && fetch_soil) {
    if (is.null(sites))
      stop("`sites` is required when `fetch_soil = TRUE`.")
    .reject_characterization_overrides(soil_control, "sites")
    soil_raw <- do.call(fetch_soilgrids, c(list(sites = sites), soil_control))
  }
  soil_features <- .characterization_soil_features(
    soil_raw, soil_profile_control
  )

  .reject_characterization_overrides(
    kernel_control, c("weather", "soil", "environments")
  )
  kernel_result <- do.call(build_environment_kernels, c(list(
    weather = history$typical, soil = soil_features,
    environments = rownames(history$typical)
  ), kernel_control))

  integration_kernels <- kernel_result$kernels[
    intersect(c("weather", "soil"), names(kernel_result$kernels))
  ]
  if (!is.null(stability) && "weather" %in% names(integration_kernels)) {
    weather_consensus <- stability$consensus_D
    envs <- rownames(integration_kernels$weather)
    integration_kernels$weather <- .normalise_environment_kernel(
      weather_consensus[envs, envs, drop = FALSE]
    )
  }
  descriptive_consensus <- if (length(integration_kernels) == 1L)
    integration_kernels[[1L]] else
      consensus_environment_kernels(integration_kernels)
  kernel_agreement <- .environment_kernel_agreement(integration_kernels)

  relationships <- list()
  relationship_groups <- character()
  if (!is.null(stability)) {
    relationships <- stability$per_year_D
    names(relationships) <- paste0("weather_", names(relationships))
    relationship_groups <- stats::setNames(
      rep("weather", length(relationships)), names(relationships)
    )
  }
  if ("soil" %in% names(integration_kernels)) {
    relationships$soil <- integration_kernels$soil
    relationship_groups <- c(relationship_groups, soil = "soil")
  }

  .reject_characterization_overrides(
    inference_control,
    c("D", "relationships", "relationship_groups", "seed", "mode",
      "min_cluster_size")
  )
  descriptive_control <- inference_control
  descriptive_min_size <- descriptive_control$descriptive_min_cluster_size %||%
    1L
  strict_min_size <- descriptive_control$strict_min_cluster_size %||% 2L
  descriptive_control$descriptive_min_cluster_size <- NULL
  descriptive_control$strict_min_cluster_size <- NULL
  common_inference <- list(
    D = descriptive_consensus,
    relationships = if (length(relationships)) relationships else NULL,
    relationship_groups = if (length(relationships)) relationship_groups else
      NULL,
    seed = seed
  )
  strata <- do.call(infer_environmental_strata, c(
    common_inference, list(min_cluster_size = descriptive_min_size),
    descriptive_control
  ))
  strict_validation <- do.call(infer_mega_environments, c(
    common_inference,
    list(min_cluster_size = strict_min_size, mode = "mega_environment"),
    descriptive_control
  ))

  .reject_characterization_overrides(
    covariance_control, c("kernels", "target", "seed")
  )
  covariance <- do.call(calibrate_environment_covariance, c(list(
    kernels = kernel_result$kernels, target = historical_target, seed = seed
  ), covariance_control))
  list(
    historical_weather = history,
    weather_stability = stability,
    soil_raw = soil_raw,
    soil_features = soil_features,
    environmental_kernels = kernel_result,
    integration_kernels = integration_kernels,
    descriptive_consensus = descriptive_consensus,
    kernel_agreement = kernel_agreement,
    environmental_strata = strata,
    mega_environment_validation = strict_validation,
    covariance_calibration = covariance,
    provenance = list(
      weather = history$weather_provenance,
      weather_qc = history$daily_qc_provenance,
      soil = if (is.null(soil_raw)) NULL else attr(soil_raw, "provenance"),
      soil_request = if (is.null(soil_raw)) NULL else
        attr(soil_raw, "soilgrids_request")
    )
  )
}


.reject_characterization_overrides <- function(control, protected) {
  duplicate <- intersect(names(control), protected)
  if (length(duplicate))
    stop("Control list must not override: ", paste(duplicate, collapse = ", "),
         ".")
  invisible(TRUE)
}


.characterization_soil_features <- function(soil, control) {
  if (is.null(soil)) return(NULL)
  d <- as.data.frame(soil, check.names = FALSE, stringsAsFactors = FALSE)
  if ("environment" %in% names(d)) {
    env <- as.character(d$environment)
    d$environment <- NULL
    rownames(d) <- env
  }
  is_profile <- grepl(
    "__(0-5cm|5-15cm|15-30cm|30-60cm|60-100cm|100-200cm|0-30cm)__",
    names(d)
  )
  if (any(is_profile)) {
    .reject_characterization_overrides(control, "soil")
    return(do.call(soil_profile_features, c(list(soil = d), control)))
  }
  data.matrix(d)
}
