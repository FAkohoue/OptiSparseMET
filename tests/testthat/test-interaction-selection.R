# Tests for interaction selection: candidate enumeration, G x E variance-based
# screening, and the agronomic hypothesis catalogue.

# ---- list_candidate_interactions -------------------------------------------

test_that("enumeration labels type and column cost, drops environment", {
  df <- data.frame(environment = c("E1", "E2", "E3"),
                   temp = 1:3, rain = 4:6,
                   planting = c("a", "b", "a"),          # 2 levels
                   irrigation = c("x", "y", "x"))        # 2 levels
  cand <- list_candidate_interactions(df, order = 2)
  expect_false(any(grepl("environment", cand$term)))
  # numeric x numeric -> 1 column; cat x cat (2x2) -> 4; num x cat -> 2
  expect_equal(cand$est_columns[cand$term == "temp:rain"], 1)
  expect_equal(cand$est_columns[cand$term == "planting:irrigation"], 4)
  expect_equal(cand$type[cand$term == "temp:planting"], "mixed")
  # sorted cheapest first
  expect_true(!is.unsorted(cand$est_columns))
})

test_that("higher-order interactions and the two-covariate minimum", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  cand3 <- list_candidate_interactions(df, order = 3)
  expect_true("a:b:c" %in% cand3$term)
  expect_error(list_candidate_interactions(data.frame(environment = 1:3, a = 1:3)),
               "at least two")
})

# ---- screen_interactions ----------------------------------------------------

test_that("screening ranks the covariate that drives G x E on top", {
  set.seed(1)
  envs <- paste0("E", 1:6)
  x <- c(-2, -1, 0, 1, 2, 3)                  # the covariate driving G x E
  slopes <- rnorm(20)                          # genotype sensitivities to x
  gxe <- outer(slopes, x) + matrix(rnorm(120, 0, 0.05), 20)
  dimnames(gxe) <- list(paste0("L", 1:20), envs)
  cov <- cbind(x = x, noise = rnorm(6)); rownames(cov) <- envs

  sc <- screen_interactions(gxe, cov, method = "anova")
  expect_equal(sc$covariate[1], "x")           # x explains G x E best
  expect_gt(sc$importance[sc$covariate == "x"],
            sc$importance[sc$covariate == "noise"])
  expect_true(all(sc$importance >= 0 & sc$importance <= 1 + 1e-8))
})

test_that("screen_interactions needs enough shared environments", {
  gxe <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("E1", "E2")))
  cov <- matrix(rnorm(2), 2, 1, dimnames = list(c("E1", "E2"), "x"))
  expect_error(screen_interactions(gxe, cov), "at least 3 shared")
})

# ---- agronomic_interaction_hypotheses --------------------------------------

test_that("the agronomic catalogue has the expected structure", {
  h <- agronomic_interaction_hypotheses()
  expect_true(all(c("domain", "factor_1", "factor_2", "rationale") %in% names(h)))
  expect_gte(nrow(h), 8L)
  expect_true(any(grepl("water", h$factor_2, ignore.case = TRUE)))
})
