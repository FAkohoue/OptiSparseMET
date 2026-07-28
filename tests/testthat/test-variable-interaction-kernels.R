test_that("dedicated interactions represent within- and cross-modality terms", {
  env <- paste0("E", 1:8)
  weather <- cbind(
    mean_temperature = c(22, 24, 23, 27, 26, 29, 25, 28),
    rainfall = c(610, 430, 720, 310, 520, 260, 650, 390)
  )
  soil <- cbind(
    bulk_density = c(1.12, 1.34, 1.21, 1.48, 1.28, 1.41, 1.17, 1.37),
    pH = c(5.3, 6.1, 5.8, 6.5, 5.6, 6.2, 5.1, 5.9)
  )
  rownames(weather) <- rownames(soil) <- env
  definitions <- data.frame(
    interaction = c(
      "temperature_x_rainfall", "temperature_x_rainfall",
      "temperature_x_density", "temperature_x_density"
    ),
    parent = c("temperature", "rainfall", "temperature", "density"),
    modality = c("weather", "weather", "weather", "soil"),
    variable = c(
      "mean_temperature", "rainfall",
      "mean_temperature", "bulk_density"
    )
  )

  out <- build_variable_interaction_kernels(
    list(weather = weather, soil = soil),
    definitions
  )

  expect_true(all(c(
    "variable_x__temperature_x_rainfall",
    "variable_x__temperature_x_density"
  ) %in% names(out$kernels)))
  expect_equal(
    out$audit$scope[
      out$audit$interaction == "temperature_x_rainfall"
    ],
    "within_modality"
  )
  expect_equal(
    out$audit$scope[
      out$audit$interaction == "temperature_x_density"
    ],
    "cross_modality"
  )
  expect_true(all(out$audit$orthogonalised))
  expect_true(all(out$audit$residual_information_fraction > 0))
  expect_true(all(out$audit$hierarchy == "strong"))
  expect_true(all(out$variable_ledger$matched_after_qc))
  expect_true(all(vapply(out$kernels, function(K)
    min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) >= -1e-7,
    logical(1))))
})

test_that("strong heredity returns exact parent kernels and metadata", {
  env <- paste0("E", 1:7)
  weather <- cbind(
    temperature = c(20, 24, 21, 27, 23, 29, 25)
  )
  soil <- cbind(
    density = c(1.1, 1.4, 1.2, 1.5, 1.3, 1.45, 1.25)
  )
  rownames(weather) <- rownames(soil) <- env
  definition <- list(
    temperature_density = list(
      temperature = list(
        modality = "weather", variables = "temperature"
      ),
      density = list(modality = "soil", variables = "density")
    )
  )

  out <- build_variable_interaction_kernels(
    list(weather = weather, soil = soil), definition
  )
  interaction <- out$kernels$variable_x__temperature_density
  parents <- attr(interaction, "interaction_parents")
  expect_length(parents, 2L)
  expect_true(all(parents %in% names(out$kernels)))
  expect_equal(attr(interaction, "interaction_class"), "variable_level")
  expect_equal(
    attr(interaction, "interaction_mode"),
    "residualised_variable_tensor"
  )
  expect_true(isTRUE(attr(interaction, "strong_heredity")))
  expect_equal(out$hierarchy$variable_x__temperature_density, parents)
})

test_that("interaction audit rejects absent and aliased variables", {
  env <- paste0("E", 1:6)
  weather <- cbind(temperature = c(20, 21, 23, 22, 25, 24))
  soil <- cbind(density = c(1.1, 1.3, 1.2, 1.5, 1.4, 1.6))
  rownames(weather) <- rownames(soil) <- env

  absent <- data.frame(
    interaction = c("temperature_density", "temperature_density"),
    parent = c("temperature", "density"),
    modality = c("weather", "soil"),
    variable = c("temperature", "bulk_density")
  )
  expect_error(
    build_variable_interaction_kernels(
      list(weather = weather, soil = soil), absent
    ),
    "absent after environmental quality control"
  )

  too_large <- data.frame(
    interaction = c("temperature_density", "temperature_density"),
    parent = c("temperature", "density"),
    modality = c("weather", "soil"),
    variable = c("temperature", "density")
  )
  expect_error(
    build_variable_interaction_kernels(
      list(weather = weather, soil = soil), too_large,
      max_tensor_columns = 0
    ),
    "positive integer"
  )
})

test_that("build_environment_kernels integrates dedicated hypotheses", {
  env <- paste0("E", 1:8)
  weather <- data.frame(
    environment = env,
    temperature = c(22, 25, 23, 28, 24, 29, 26, 27),
    rainfall = c(500, 420, 610, 300, 570, 260, 450, 390)
  )
  soil <- data.frame(
    environment = env,
    density = c(1.1, 1.4, 1.2, 1.5, 1.3, 1.45, 1.25, 1.36)
  )
  definition <- data.frame(
    interaction = c("temperature_density", "temperature_density"),
    parent = c("temperature", "density"),
    modality = c("weather", "soil"),
    variable = c("temperature", "density")
  )

  out <- build_environment_kernels(
    weather = weather,
    soil = soil,
    include_interactions = TRUE,
    variable_interactions = definition
  )

  expect_true(all(c(
    "weather_x_soil",
    "variable_x__temperature_density"
  ) %in% names(out$kernels)))
  expect_equal(nrow(out$variable_interaction_audit), 1L)
  expect_equal(nrow(out$variable_ledger), 2L)
  expect_true(
    "variable_x__temperature_density" %in% names(out$kernel_hierarchy)
  )
})

test_that("variable interactions remain unweighted without historical MET", {
  env <- paste0("E", 1:7)
  weather <- cbind(temperature = c(20, 24, 21, 27, 23, 29, 25))
  soil <- cbind(density = c(1.1, 1.4, 1.2, 1.5, 1.3, 1.45, 1.25))
  rownames(weather) <- rownames(soil) <- env
  definition <- data.frame(
    interaction = c("temperature_density", "temperature_density"),
    parent = c("temperature", "density"),
    modality = c("weather", "soil"),
    variable = c("temperature", "density")
  )
  built <- build_variable_interaction_kernels(
    list(weather = weather, soil = soil), definition
  )

  evidence <- assess_variable_interactions(built$kernels, n_boot = 0L)
  expect_equal(evidence$evidence$decision, "structural_uncertainty")
  expect_true(is.na(evidence$evidence$supported))
  expect_length(evidence$supported, 0L)
  expect_equal(
    evidence$structural_uncertainty,
    "variable_x__temperature_density"
  )

  covariance <- calibrate_environment_covariance(
    built$kernels, n_boot = 0L
  )
  expect_equal(covariance$status, "no_historical_met")
  expect_null(covariance$weights)
  expect_true(
    "kernel_variable_x__temperature_density" %in%
      names(covariance$candidates)
  )
})

test_that("component screening controls central interaction eligibility", {
  env <- paste0("E", 1:8)
  weather <- cbind(
    temperature = c(20, 24, 21, 27, 23, 29, 25, 26)
  )
  soil <- cbind(
    density = c(1.1, 1.4, 1.2, 1.5, 1.3, 1.45, 1.25, 1.36)
  )
  rownames(weather) <- rownames(soil) <- env
  definition <- data.frame(
    interaction = c("temperature_density", "temperature_density"),
    parent = c("temperature", "density"),
    modality = c("weather", "soil"),
    variable = c("temperature", "density")
  )
  built <- build_variable_interaction_kernels(
    list(weather = weather, soil = soil), definition
  )
  target <- built$kernels[[built$hierarchy[[1L]][1L]]]

  evidence <- assess_variable_interactions(
    built$kernels, target = target, n_boot = 0L
  )
  expect_true(all(c(
    "mean_rmse_improvement", "adjusted_p", "fitted_weight",
    "bootstrap_inclusion", "decision"
  ) %in% names(evidence$evidence)))

  central <- calibrate_environment_covariance(
    built$kernels,
    target = target,
    eligible_interactions = character(),
    n_boot = 0L
  )
  expect_equal(central$interaction_status, "interactions_screened_out")
  expect_equal(
    unname(central$weights["variable_x__temperature_density"]),
    0
  )
  expect_true(
    "kernel_variable_x__temperature_density" %in%
      names(central$candidates)
  )
})
