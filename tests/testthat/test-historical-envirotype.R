# Tests for planning-stage historical envirotyping: across-year aggregation into
# a typical profile plus interannual variability. Offline (daily_by_year path).

test_that(".aggregate_across_years computes typical mean and variability SD", {
  m1 <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("E1", "E2"), c("a", "b")))
  m2 <- matrix(c(3, 4, 5, 6), 2, dimnames = list(c("E1", "E2"), c("a", "b")))
  m3 <- matrix(c(2, 3, 4, 5), 2, dimnames = list(c("E1", "E2"), c("a", "b")))
  agg <- .aggregate_across_years(list(m1, m2, m3), "sd")
  expect_equal(agg$typical["E1", "a"], mean(c(1, 3, 2)))          # 2
  expect_equal(agg$variability["E1", "a"], stats::sd(c(1, 3, 2)))
  expect_equal(dim(agg$typical), c(2L, 2L))
})

test_that("cv variability divides SD by the typical mean", {
  m1 <- matrix(c(2, 10), 1, dimnames = list("E1", c("a", "b")))
  m2 <- matrix(c(4, 10), 1, dimnames = list("E1", c("a", "b")))
  agg <- .aggregate_across_years(list(m1, m2), "cv")
  expect_equal(agg$variability["E1", "a"], stats::sd(c(2, 4)) / 3)  # mean=3
  expect_equal(agg$variability["E1", "b"], 0)                        # constant -> cv 0
})

mk_year <- function(base) rbind(
  data.frame(environment = "E1", day = 0:59,
             T2M = base + rnorm(60, 0, 0.5), T2M_MAX = base + 6 + rnorm(60, 0, 0.5),
             stringsAsFactors = FALSE),
  data.frame(environment = "E2", day = 0:59,
             T2M = base + 4 + rnorm(60, 0, 0.5), T2M_MAX = base + 10 + rnorm(60, 0, 0.5),
             stringsAsFactors = FALSE))

test_that("historical_envirotype summarises several past seasons (offline)", {
  set.seed(1)
  hist <- historical_envirotype(
    years = 2020:2022,
    daily_by_year = list(`2020` = mk_year(24), `2021` = mk_year(27),
                         `2022` = mk_year(22)),
    envirotype = list(windows = 2, stats = c("mean", "sd")))
  expect_equal(hist$n_years, 3L)
  expect_equal(rownames(hist$typical), c("E1", "E2"))
  expect_true("T2M_mean_w1" %in% colnames(hist$typical))
  # E2 is the hotter site in the typical profile
  expect_gt(hist$typical["E2", "T2M_mean_w1"], hist$typical["E1", "T2M_mean_w1"])
  # interannual variability is defined (3 seasons) and non-negative
  expect_true(all(is.finite(hist$variability[, "T2M_mean_w1"])))
  expect_true(all(hist$variability[, "T2M_mean_w1"] >= 0))
})

test_that("the typical profile plugs into the environment relationship", {
  set.seed(2)
  hist <- historical_envirotype(
    years = 2019:2021,
    daily_by_year = list(`2019` = mk_year(23), `2020` = mk_year(25),
                         `2021` = mk_year(24)),
    envirotype = list(windows = 2, stats = "mean"))
  D <- build_environment_relationship(hist$typical, source = "enviromic",
                                      kernel = "gaussian")
  expect_equal(dim(D), c(2L, 2L))
})

test_that("assess_envirotype_stability builds per-year D and reports stability", {
  mk3 <- function(base) rbind(
    data.frame(environment = "E1", day = 0:49,
               T2M = base + rnorm(50, 0, 0.3), T2M_MAX = base + 6 + rnorm(50, 0, 0.3)),
    data.frame(environment = "E2", day = 0:49,
               T2M = base + 3 + rnorm(50, 0, 0.3), T2M_MAX = base + 9 + rnorm(50, 0, 0.3)),
    data.frame(environment = "E3", day = 0:49,
               T2M = base - 2 + rnorm(50, 0, 0.3), T2M_MAX = base + 3 + rnorm(50, 0, 0.3)))
  set.seed(1)
  hist <- historical_envirotype(
    years = 2019:2021,
    daily_by_year = list(`2019` = mk3(24), `2020` = mk3(26), `2021` = mk3(22)),
    envirotype = list(windows = 2, stats = c("mean", "sd")))
  st <- assess_envirotype_stability(hist, kernel = "gaussian")
  expect_length(st$per_year_D, 3L)
  expect_equal(dim(st$stability), c(3L, 3L))
  expect_equal(diag(st$stability), rep(1, 3), ignore_attr = TRUE)
  expect_true(is.finite(st$mean_stability))
  # consensus D = mean of per-year matrices, same shape, symmetric
  expect_equal(dim(st$consensus_D), dim(st$per_year_D[[1]]))
  expect_equal(st$consensus_D, t(st$consensus_D), tolerance = 1e-8)
  # `combined` binds typical + interannual variability (_iav) columns
  expect_true(any(grepl("_iav$", colnames(hist$combined))))
  expect_equal(nrow(hist$combined), nrow(hist$typical))
})

test_that("network path requires sites, window and years", {
  expect_error(historical_envirotype(years = 2020),
               "Supply `sites`")
  expect_error(
    historical_envirotype(sites = data.frame(environment = "E1",
                                             latitude = 9, longitude = 7),
                          years = 2020, window = "06-01"),
    "start_md, end_md")
})

test_that(".resolve_windows supports global, per-site data frame and columns", {
  sites <- data.frame(environment = c("A", "B"), latitude = 1:2, longitude = 3:4)
  g <- .resolve_windows(sites, c("06-01", "09-30"))
  expect_equal(g$start_md, c("06-01", "06-01"))
  # per-site data frame (B has a different, year-crossing window)
  wdf <- data.frame(environment = c("A", "B"),
                    start_md = c("06-01", "11-01"), end_md = c("09-30", "02-28"))
  ps <- .resolve_windows(sites, wdf)
  expect_equal(ps$end_md[ps$environment == "B"], "02-28")
  expect_true(ps$end_md[2] < ps$start_md[2])            # -> flagged as cross-year
  # columns on sites
  s2 <- cbind(sites, start_md = c("05-01", "05-01"), end_md = c("08-01", "08-01"))
  expect_equal(.resolve_windows(s2, NULL)$start_md, c("05-01", "05-01"))
})

test_that(".merge_daily_series lets station data override and fill", {
  dl <- data.frame(environment = "E1", day = 0:2,
                   T2M = c(20, 21, 22), RAD = c(200, 210, 205))   # downloaded
  st <- data.frame(environment = "E1", day = 0:2,
                   T2M = c(19.5, NA, 22.5))                        # station (a gap)
  m <- .merge_daily_series(dl, st)
  expect_equal(m$T2M, c(19.5, 21, 22.5))       # station wins, downloaded fills the NA
  expect_true("RAD" %in% names(m))             # station-missing variable retained
})
