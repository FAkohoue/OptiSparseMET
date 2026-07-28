test_that("environmental QC records ranges, coverage, and imputation", {
  d <- data.frame(
    environment = rep(c("E1", "E2"), each = 3),
    date = rep(as.Date("2024-01-01") + 0:2, 2),
    T2M = c(20, NA, 22, 200, 25, 26),
    RAIN = c(0, 1, 2, 0, NA, 3)
  )
  q <- qc_environmental_data(
    d, environments = c("E1", "E2"), date_col = "date",
    ranges = list(T2M = c(-20, 60), RAIN = c(0, 500)),
    min_coverage = 0.5, missing_action = "impute", impute = "linear",
    source = "weather"
  )
  expect_false(anyNA(q$data$T2M))
  expect_true(any(q$audit$n_invalid > 0))
  expect_true(nrow(q$imputation) >= 2L)
  expect_true(all(c("T2M__was_missing", "RAIN__was_missing") %in%
                    names(q$data)))
})

test_that("environmental QC rejects duplicate environment-date keys", {
  d <- data.frame(environment = c("E1", "E1"),
                  date = as.Date(c("2024-01-01", "2024-01-01")),
                  T2M = c(20, 21))
  expect_error(qc_environmental_data(d, date_col = "date"), "duplicate")
})

test_that("supplied covariates are filled cellwise from fetched covariates", {
  supplied <- matrix(c(1, NA, 3, 4), 2,
                     dimnames = list(c("E1", "E2"), c("a", "b")))
  fetched <- matrix(c(10, 20, 30, 40, 50, 60), 2,
                    dimnames = list(c("E1", "E2"), c("a", "b", "c")))
  z <- .merge_cov_blocks(supplied, fetched)
  expect_equal(z["E2", "a"], 20)
  expect_equal(z["E1", "a"], 1)
  expect_true("c" %in% names(z))
})

test_that("standalone environment relationship does not silently impute", {
  X <- matrix(c(1, NA, 3, 4, 5, 6), 3,
              dimnames = list(paste0("E", 1:3), c("a", "b")))
  expect_error(
    build_environment_relationship(X, source = "enviromic"),
    "missing values"
  )
  D <- build_environment_relationship(
    X, source = "enviromic", kernel = "gaussian", missing = "median"
  )
  expect_equal(diag(D), rep(1, 3), tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("modality kernels are separate, normalised, and PSD", {
  env <- paste0("E", 1:6)
  weather <- matrix(rnorm(6 * 30), 6, dimnames = list(env, NULL))
  soil <- matrix(rnorm(6 * 3), 6, dimnames = list(env, NULL))
  management <- data.frame(
    environment = env,
    irrigation = rep(c("rainfed", "irrigated"), 3),
    nitrogen = c(80, 100, 90, 110, 85, 105)
  )
  geo <- data.frame(environment = env, latitude = 1:6,
                    longitude = seq(-80, -75))
  out <- build_environment_kernels(
    weather = weather, soil = soil, management = management,
    geography = geo, environments = env, include_interactions = TRUE
  )
  expect_true(all(c("weather", "soil", "management", "geography",
                    "weather_x_soil", "weather_x_management",
                    "soil_x_management") %in%
                    names(out$kernels)))
  expect_true(all(vapply(
    out$kernels[c("weather_x_soil", "weather_x_management",
                  "soil_x_management")],
    function(K) identical(attr(K, "interaction_mode"), "anova"),
    logical(1)
  )))
  for (K in out$kernels) {
    expect_equal(diag(K), rep(1, 6), tolerance = 1e-8,
                 ignore_attr = TRUE)
    expect_gte(min(eigen(K, symmetric = TRUE,
                         only.values = TRUE)$values), -1e-7)
  }
  expect_equal(out$block_diagnostics$n_variables[
    out$block_diagnostics$block == "weather"], 30)

  custom <- build_environment_kernels(
    weather = weather, soil = soil, management = management,
    environments = env,
    interaction_terms = list(
      weather_x_soil_x_management =
        c("weather", "soil", "management")
    ),
    interaction_mode = "product"
  )
  expect_true("weather_x_soil_x_management" %in% names(custom$kernels))
  expect_false("weather_x_soil" %in% names(custom$kernels))
  expect_equal(
    attr(custom$kernels$weather_x_soil_x_management, "interaction_mode"),
    "product"
  )
})

test_that("ANOVA interactions are centred and higher orders are explicit", {
  env <- paste0("E", 1:6)
  K_weather <- outer(1:6, 1:6, function(i, j) exp(-abs(i - j) / 2))
  K_soil <- outer(c(1, 3, 2, 6, 4, 5), c(1, 3, 2, 6, 4, 5),
                  function(i, j) exp(-abs(i - j) / 2))
  K_management <- outer(rep(1:2, 3), rep(1:2, 3),
                        function(i, j) ifelse(i == j, 1, 0.25))
  dimnames(K_weather) <- dimnames(K_soil) <-
    dimnames(K_management) <- list(env, env)
  main <- list(weather = K_weather, soil = K_soil,
               management = K_management)

  pairwise <- add_environment_kernel_interactions(main)
  expect_true(all(c("weather_x_soil", "weather_x_management",
                    "soil_x_management") %in% names(pairwise)))
  expect_equal(attr(pairwise$weather_x_soil, "interaction_parents"),
               c("weather", "soil"))
  expect_equal(attr(pairwise$weather_x_soil, "interaction_mode"), "anova")

  raw <- add_environment_kernel_interactions(
    main, interactions = list(weather_x_soil = c("weather", "soil")),
    mode = "product"
  )
  expect_false(isTRUE(all.equal(
    unname(pairwise$weather_x_soil),
    unname(raw$weather_x_soil),
    tolerance = 1e-7
  )))

  triple <- add_environment_kernel_interactions(
    main,
    interactions = list(
      weather_x_soil_x_management =
        c("weather", "soil", "management")
    )
  )
  expect_true("weather_x_soil_x_management" %in% names(triple))
  expect_equal(
    attr(triple$weather_x_soil_x_management, "interaction_parents"),
    c("weather", "soil", "management")
  )
  expect_gte(min(eigen(triple$weather_x_soil_x_management,
                       symmetric = TRUE, only.values = TRUE)$values), -1e-7)
})

test_that("interaction definitions cannot overwrite or name absent parents", {
  env <- paste0("E", 1:4)
  K <- diag(4); dimnames(K) <- list(env, env)
  expect_error(
    add_environment_kernel_interactions(
      list(weather = K, soil = K),
      interactions = list(weather = c("weather", "soil"))
    ),
    "must not overwrite"
  )
  expect_error(
    add_environment_kernel_interactions(
      list(weather = K, soil = K),
      interactions = list(weather_x_management =
                            c("weather", "management"))
    ),
    "existing kernels"
  )
})

test_that("kernel combination honours modality weights rather than column count", {
  env <- paste0("E", 1:4)
  K1 <- diag(4); dimnames(K1) <- list(env, env)
  K2 <- matrix(0.5, 4, 4); diag(K2) <- 1
  dimnames(K2) <- list(env, env)
  D <- combine_environment_kernels(
    list(weather = K1, soil = K2),
    weights = c(weather = 0.25, soil = 0.75)
  )
  expect_equal(D[1, 2], 0.375, tolerance = 1e-8)
  expect_equal(attr(D, "weights")[c("weather", "soil")],
                c(weather = 0.25, soil = 0.75))
})

test_that("multi-kernel combination never assigns implicit equal weights", {
  env <- paste0("E", 1:3)
  K1 <- K2 <- diag(3)
  dimnames(K1) <- dimnames(K2) <- list(env, env)
  expect_error(
    combine_environment_kernels(list(weather = K1, soil = K2)),
    "weights.*required"
  )
})

test_that("descriptive kernel consensus requires no modality weights", {
  env <- paste0("E", 1:5)
  K1 <- outer(1:5, 1:5, function(i, j) exp(-abs(i - j)))
  K2 <- outer(c(1, 4, 2, 5, 3), c(1, 4, 2, 5, 3),
              function(i, j) exp(-abs(i - j)))
  dimnames(K1) <- dimnames(K2) <- list(env, env)
  D <- consensus_environment_kernels(list(weather = K1, soil = K2))
  expect_equal(diag(D), rep(1, 5), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_gte(min(eigen(D, symmetric = TRUE,
                       only.values = TRUE)$values), -1e-7)
  expect_equal(attr(D, "integration"),
               "entrywise_median_of_within_kernel_ranks")
  expect_null(attr(D, "weights"))
})

test_that("historical correlations calibrate non-negative simplex weights", {
  env <- paste0("E", 1:5)
  K_weather <- outer(1:5, 1:5, function(i, j) exp(-abs(i - j)))
  K_soil <- outer(c(1, 3, 2, 5, 4), c(1, 3, 2, 5, 4),
                  function(i, j) exp(-abs(i - j)))
  dimnames(K_weather) <- dimnames(K_soil) <- list(env, env)
  target <- 0.8 * K_weather + 0.1 * K_soil + 0.1 * diag(5)
  dimnames(target) <- list(env, env)
  fit <- calibrate_environment_covariance(
    list(weather = K_weather, soil = K_soil), target = target,
    ridge = 0, n_boot = 20, seed = 1
  )
  expect_equal(fit$status, "historically_calibrated")
  expect_equal(sum(fit$weights), 1, tolerance = 1e-8)
  expect_true(all(fit$weights >= 0))
  expect_gt(fit$weights["weather"], fit$weights["soil"])
  expect_true(all(c("validation", "held_out", "rmse") %in%
                    names(fit$cross_validation)))
})

test_that("no-history covariance is identity and does not define weights", {
  env <- paste0("E", 1:4)
  K <- matrix(0.7, 4, 4); diag(K) <- 1
  dimnames(K) <- list(env, env)
  out1 <- calibrate_environment_covariance(
    list(weather = K), n_boot = 0
  )
  out2 <- calibrate_environment_covariance(
    list(weather = K),
    prior_weights = c(weather = 0.01, identity = 0.99),
    identity_weight = 0.99, n_boot = 0
  )
  expect_equal(out1$status, "no_historical_met")
  expect_equal(out1$Sigma_E, diag(4), ignore_attr = TRUE)
  expect_equal(out1$Sigma_E, out2$Sigma_E)
  expect_null(out1$weights)
  expect_true(all(c("independent", "kernel_weather") %in%
                    names(out1$candidates)))
})

test_that("interactions remain uncertainty scenarios without historical MET", {
  env <- paste0("E", 1:6)
  K_weather <- outer(1:6, 1:6, function(i, j) exp(-abs(i - j) / 2))
  K_soil <- outer(c(1, 3, 2, 6, 4, 5), c(1, 3, 2, 6, 4, 5),
                  function(i, j) exp(-abs(i - j) / 2))
  dimnames(K_weather) <- dimnames(K_soil) <- list(env, env)
  kernels <- add_environment_kernel_interactions(
    list(weather = K_weather, soil = K_soil)
  )

  out <- calibrate_environment_covariance(kernels, n_boot = 0L)
  expect_equal(out$interaction_status,
               "unweighted_interaction_uncertainty")
  expect_true("kernel_weather_x_soil" %in% names(out$candidates))
  expect_null(out$weights)
  expect_false(any(out$interaction_evidence$model_activated))
  expect_false(any(out$interaction_evidence$component_active))
  expect_match(out$interaction_evidence$reason,
               "unweighted structural uncertainty")

  sensitivity <- calibrate_environment_covariance(
    kernels, interaction_policy = "exclude", n_boot = 0L
  )
  expect_true("kernel_weather_x_soil" %in% names(sensitivity$candidates))
  expect_false(any(sensitivity$interaction_evidence$model_activated))
})

test_that("historical calibration rejects unsupported interactions", {
  env <- paste0("E", 1:7)
  K_weather <- outer(1:7, 1:7, function(i, j) exp(-abs(i - j) / 2))
  K_soil <- outer(c(1, 3, 2, 7, 4, 6, 5),
                  c(1, 3, 2, 7, 4, 6, 5),
                  function(i, j) exp(-abs(i - j) / 2))
  dimnames(K_weather) <- dimnames(K_soil) <- list(env, env)
  kernels <- add_environment_kernel_interactions(
    list(weather = K_weather, soil = K_soil)
  )

  out <- calibrate_environment_covariance(
    kernels, target = K_weather, ridge = 0,
    interaction_policy = "evidence", interaction_alpha = 1e-12,
    n_boot = 0L
  )
  expect_equal(out$interaction_status, "interaction_model_rejected")
  expect_equal(unname(out$weights["weather_x_soil"]), 0,
               tolerance = 1e-10)
  expect_true("kernel_weather_x_soil" %in% names(out$candidates))
  expect_true(all(c("model", "validation", "held_out", "rmse") %in%
                    names(out$cross_validation)))
  expect_match(out$interaction_evidence$reason,
               "did not pass|could not compare")
})

test_that("interaction inclusion is an explicit historical-data override", {
  env <- paste0("E", 1:7)
  K_weather <- outer(1:7, 1:7, function(i, j) exp(-abs(i - j) / 2))
  K_soil <- outer(c(1, 3, 2, 7, 4, 6, 5),
                  c(1, 3, 2, 7, 4, 6, 5),
                  function(i, j) exp(-abs(i - j) / 2))
  dimnames(K_weather) <- dimnames(K_soil) <- list(env, env)
  kernels <- add_environment_kernel_interactions(
    list(weather = K_weather, soil = K_soil)
  )
  out <- calibrate_environment_covariance(
    kernels, target = K_weather, ridge = 0,
    interaction_policy = "include", n_boot = 0L
  )
  expect_equal(out$interaction_status, "interaction_model_activated")
  expect_true("kernel_weather_x_soil" %in% names(out$candidates))
  expect_true(all(out$interaction_evidence$model_activated))
})

test_that("robust scenarios can carry alternative environment covariances", {
  env <- c("E1", "E2")
  A <- diag(2); B <- matrix(c(1, 0.5, 0.5, 1), 2)
  dimnames(A) <- dimnames(B) <- list(env, env)
  sc <- robust_scenarios(
    sigma_e2 = c(1, 2), sigmaE_shrink = c(0, 0.5),
    Sigma_E_candidates = list(independent = A, correlated = B)
  )
  expect_length(sc, 8L)
  expect_true(all(vapply(sc, function(z) !is.null(z$Sigma_E), logical(1))))
  expect_equal(sum(vapply(sc, `[[`, numeric(1), "prob")), 1)
})

test_that("actual response dates are extracted and ordered", {
  d <- data.frame(YEAR = c(2024, 2024), MO = c(1, 1), DY = c(2, 1))
  expect_equal(.power_response_dates(d),
               as.Date(c("2024-01-02", "2024-01-01")))
})

test_that("phenology and rice stress features use observed stages", {
  dates <- as.Date("2024-01-01") + 0:9
  d <- data.frame(
    environment = "E1", date = dates,
    T2M = rep(28, 10), T2M_MAX = c(rep(34, 5), rep(38, 5)),
    T2M_MIN = c(rep(22, 5), rep(26, 5)),
    RH2M = 70, PRECTOTCORR = c(rep(0, 5), rep(25, 5)),
    EVPTRNS = 4, GWETROOT = c(rep(0.2, 5), rep(0.5, 5))
  )
  ph <- data.frame(
    environment = c("E1", "E1"), stage = c("veg", "flower"),
    start_date = dates[c(1, 6)], end_date = dates[c(5, 10)]
  )
  W <- rice_weather_features(d, phenology = ph, stats = c("mean", "sum"))
  expect_true(all(c("VPD_mean_veg", "heat_days_flower",
                    "hot_nights_flower", "dry_spell_longest_veg",
                    "root_water_deficit_severity_veg") %in% colnames(W)))
  expect_equal(W["E1", "heat_days_flower"], 5)
  expect_equal(W["E1", "hot_nights_flower"], 5)
})

test_that("weather calibration gives station values precedence", {
  dl <- data.frame(environment = "E1", day = 0:5,
                   T2M = 20:25, PRECTOTCORR = 1:6)
  st <- data.frame(environment = "E1", day = 0:5,
                   T2M = 22:27, PRECTOTCORR = 2 * (1:6))
  z <- calibrate_weather_series(dl, st, min_overlap = 3)
  expect_equal(z$data$T2M, st$T2M)
  expect_equal(z$data$PRECTOTCORR, st$PRECTOTCORR)
  expect_true(all(z$diagnostics$n_overlap == 6))
})

test_that("soil profiles are thickness weighted and include uncertainty", {
  env <- c("E1", "E2")
  mk <- function(property, depth, q, stat = "mean")
    paste(property, depth, q, stat, sep = "__")
  S <- data.frame(row.names = env, check.names = FALSE)
  for (p in c("sand", "silt", "clay", "wv0033", "wv1500")) {
    for (d in c("0-5cm", "5-15cm", "15-30cm")) {
      for (q in c("Q0.05", "Q0.5", "Q0.95"))
        S[[mk(p, d, q)]] <- if (p == "wv0033") c(30, 35) else
          if (p == "wv1500") c(12, 15) else
            c(sand = 50, silt = 30, clay = 20)[p] +
              c(0, 2) + c(`Q0.05` = -2, `Q0.5` = 0, `Q0.95` = 2)[q]
    }
  }
  F <- soil_profile_features(S, root_depth_cm = 30)
  expect_true(all(c("available_water_capacity_root",
                    "texture_ilr_sand_vs_silt",
                    "clay_root_prediction_width") %in% colnames(F)))
  expect_equal(unname(F[, "available_water_capacity_root"]), c(18, 20))
  expect_true(all(F[, "clay_root_prediction_width"] > 0))
})
