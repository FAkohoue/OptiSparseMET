# Tests for configurable weather parameters/aggregation, the supplied+fetched
# merge logic, and the extended variable catalogue. The network fetch itself is
# covered (opt-in) in test-network-fetch.R; here the pure logic is tested offline.

# ---- weather aggregation ----------------------------------------------------

test_that(".aggregate_weather applies per-parameter statistics and friendly names", {
  d <- data.frame(T2M = c(20, 22, 24), PRECTOTCORR = c(5, 0, 10),
                  WS2M = c(1, 2, 3), EVPTRNS = c(2, 2, 2))
  stat_map <- list(T2M = "mean", PRECTOTCORR = "sum", WS2M = "max",
                   EVPTRNS = "sum")
  agg <- .aggregate_weather(d, stat_map)
  expect_equal(agg[["mean_temp"]], 22)          # T2M -> friendly name, mean
  expect_equal(agg[["total_precip"]], 15)        # PRECTOTCORR -> friendly, sum
  expect_equal(agg[["WS2M"]], 3)                 # custom param keeps code, max
  expect_equal(agg[["EVPTRNS"]], 6)              # sum
})

test_that(".aggregate_weather returns NA for a requested but absent parameter", {
  d <- data.frame(T2M = c(20, 22))
  agg <- .aggregate_weather(d, list(T2M = "mean", WS2M = "mean"))
  expect_equal(agg[["mean_temp"]], 21)
  expect_true(is.na(agg[["WS2M"]]))
})

# ---- supplied + fetched merge ----------------------------------------------

test_that(".merge_cov_blocks keeps user columns and adds only the missing ones", {
  own <- matrix(c(21, 24, 400, 300), 2, 2,
                dimnames = list(c("E1", "E2"), c("mean_temp", "total_precip")))
  fetched <- matrix(c(20, 25, 2, 3, 0.4, 0.5), 2, 3,
                    dimnames = list(c("E1", "E2"),
                                    c("mean_temp", "WS2M", "GWETROOT")))
  m <- .merge_cov_blocks(own, fetched)
  expect_true(all(c("mean_temp", "total_precip", "WS2M", "GWETROOT") %in% colnames(m)))
  expect_equal(rownames(m), c("E1", "E2"))
  expect_equal(unname(m[, "mean_temp"]),
               unname(own[, "mean_temp"]))              # user's value wins
  expect_equal(unname(m[, "WS2M"]),
               unname(fetched[, "WS2M"]))                # missing var filled
})

test_that(".merge_cov_blocks handles NULL inputs", {
  a <- matrix(1, 1, 1, dimnames = list("E1", "x"))
  expect_identical(.merge_cov_blocks(NULL, a), a)
  expect_identical(.merge_cov_blocks(a, NULL), a)
})

test_that("build_enviromic_covariates keeps supplied weather columns", {
  sites <- data.frame(environment = c("E1", "E2", "E3"))
  wx <- data.frame(environment = sites$environment,
                   my_temp = c(21, 24, 27), my_rain = c(400, 300, 250))
  X <- build_enviromic_covariates(sites, weather = wx, standardize = FALSE)
  expect_true(all(c("my_temp", "my_rain") %in% colnames(X)))
})

# ---- extended catalogue -----------------------------------------------------

test_that("the weather catalogue documents the additional POWER parameters", {
  wc <- enviromic_variable_catalog("weather")
  expect_true(all(c("mean_temp", "WS2M", "GWETROOT", "EVPTRNS", "T2MDEW") %in%
                    wc$variable))
  expect_true(all(c("variable", "description", "units") %in% names(wc)))
})
