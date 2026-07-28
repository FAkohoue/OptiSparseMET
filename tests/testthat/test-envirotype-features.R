# Tests for temporal envirotyping features (windows, variability/stress/GDD,
# functional representation). All offline: the engine operates on a daily data
# frame, so no network is involved.

mk_daily <- function() {
  set.seed(1)
  mk <- function(e, len, base) data.frame(
    environment = e, day = 0:(len - 1),
    T2M = base + sin((0:(len - 1)) / 10) * 3 + rnorm(len, 0, 0.3),
    T2M_MAX = base + 6 + rnorm(len, 0, 0.3),
    PRECTOTCORR = pmax(0, rnorm(len, 3, 3)),
    stringsAsFactors = FALSE)
  rbind(mk("E1", 90, 24), mk("E2", 80, 29))    # E2 hotter, shorter season
}

test_that("named stage windows produce per-stage statistics", {
  d <- mk_daily()
  W <- envirotype_features(d, windows = list(veg = c(0, 40), flower = c(40, 70)),
                           stats = c("mean", "sd"))
  expect_true(all(c("T2M_mean_veg", "T2M_sd_veg", "T2M_mean_flower",
                    "PRECTOTCORR_mean_flower") %in% colnames(W)))
  expect_equal(rownames(W), c("E1", "E2"))
  # E2 is hotter than E1 in every stage
  expect_gt(W["E2", "T2M_mean_veg"], W["E1", "T2M_mean_veg"])
  # sd captures within-stage variation (non-zero, finite)
  expect_true(all(is.finite(W[, "T2M_sd_veg"])) && all(W[, "T2M_sd_veg"] > 0))
})

test_that("integer K splits the season into K equal windows", {
  d <- mk_daily()
  W <- envirotype_features(d, windows = 3, stats = "mean")
  expect_true(all(paste0("T2M_mean_w", 1:3) %in% colnames(W)))
})

test_that("stress-day counts and growing degree days are computed per window", {
  d <- mk_daily()
  W <- envirotype_features(
    d, windows = 2,
    thresholds = list(heat_days = list(par = "T2M_MAX", above = 33)),
    gdd = list(gdd = list(par = "T2M", base = 10)))
  expect_true(any(grepl("^heat_days_w", colnames(W))))
  expect_true(any(grepl("^gdd_w", colnames(W))))
  # hotter E2 accumulates more GDD-per-window and more heat days than cooler E1
  expect_gte(sum(W["E2", grep("^heat_days_w", colnames(W))]),
             sum(W["E1", grep("^heat_days_w", colnames(W))]))
  expect_gt(sum(W["E2", grep("^gdd_w", colnames(W))]),
            sum(W["E1", grep("^gdd_w", colnames(W))]))
})

test_that("functional representation returns df coefficients per parameter", {
  d <- mk_daily()
  W <- envirotype_features(d, windows = 1,
                           functional = list(temp_shape = list(par = "T2M", df = 4)))
  fcols <- grep("^temp_shape_f", colnames(W), value = TRUE)
  expect_length(fcols, 4L)
  expect_true(all(is.finite(W[, fcols])))
})

test_that("day index is derived from row order when absent", {
  d <- mk_daily(); d$day <- NULL
  W <- NULL
  expect_warning(
    W <- envirotype_features(d, windows = 2, stats = "mean"),
    "deriving day from row order"
  )
  expect_equal(nrow(W), 2L)
  expect_equal(rownames(W), c("E1", "E2"))
  expect_true(all(is.finite(W[, "T2M_mean_w1"])))
})

test_that("envirotype features feed the environment relationship matrix", {
  d <- mk_daily()
  W <- envirotype_features(d, windows = 2, stats = c("mean", "sd"))
  D <- build_environment_relationship(W, source = "enviromic", kernel = "gaussian")
  expect_equal(dim(D), c(2L, 2L))
  expect_equal(diag(D), rep(1, 2), tolerance = 1e-8, ignore_attr = TRUE)
})
