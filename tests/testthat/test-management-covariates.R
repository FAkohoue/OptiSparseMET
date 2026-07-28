# Tests for management-covariate refinements: not scaling dummies, and building
# management interaction covariates.

# ---- scale_dummies ----------------------------------------------------------

test_that("scale_dummies = FALSE leaves 0/1 columns untouched but scales the rest", {
  sites <- data.frame(environment = c("E1", "E2", "E3", "E4"))
  wx <- data.frame(environment = sites$environment,
                   mean_temp = c(18, 22, 26, 30))          # continuous
  mg <- data.frame(environment = sites$environment,
                   planting_mode = c("transplanting", "direct",
                                     "transplanting", "direct"),
                   fert_dose = c(120, 80, 100, 60))         # continuous

  X <- build_enviromic_covariates(sites, weather = wx, management = mg,
                                  standardize = TRUE, scale_dummies = FALSE)
  dummies <- grep("planting_mode_", colnames(X), value = TRUE)
  expect_true(length(dummies) == 2L)
  # dummy columns remain exactly 0/1
  expect_true(all(as.vector(X[, dummies]) %in% c(0, 1)))
  # continuous columns are standardised (mean ~ 0)
  cont <- setdiff(colnames(X), dummies)
  expect_equal(unname(colMeans(X[, cont])), rep(0, length(cont)), tolerance = 1e-8)
})

test_that("scale_dummies = TRUE (default) scales everything", {
  sites <- data.frame(environment = c("E1", "E2", "E3", "E4"))
  mg <- data.frame(environment = sites$environment,
                   planting_mode = c("a", "b", "a", "b"),
                   fert_dose = c(120, 80, 100, 60))
  X <- build_enviromic_covariates(sites, management = mg, standardize = TRUE)
  # a scaled dummy is centred, so not all values are 0/1 any more
  dummies <- grep("planting_mode_", colnames(X), value = TRUE)
  expect_false(all(as.vector(X[, dummies]) %in% c(0, 1)))
  expect_equal(unname(colMeans(X)), rep(0, ncol(X)), tolerance = 1e-8)
})

# ---- add_management_interactions --------------------------------------------

test_that("categorical x numeric gives dose-within-type columns", {
  mg <- data.frame(environment = c("E1", "E2", "E3"),
                   fert_type = c("NPK", "urea", "NPK"),
                   fert_dose = c(120, 80, 100))
  out <- add_management_interactions(mg, list(c("fert_type", "fert_dose")))
  expect_true(all(c("fert_type_NPK:fert_dose", "fert_type_urea:fert_dose") %in%
                    names(out)))
  # dose appears only within its own type, 0 otherwise
  expect_equal(out[["fert_type_NPK:fert_dose"]], c(120, 0, 100))
  expect_equal(out[["fert_type_urea:fert_dose"]], c(0, 80, 0))
  # original columns retained
  expect_true(all(c("fert_type", "fert_dose") %in% names(out)))
})

test_that("numeric x numeric gives a product and errors on missing vars", {
  mg <- data.frame(environment = c("E1", "E2"),
                   dose = c(2, 3), density = c(10, 5))
  out <- add_management_interactions(mg, list(c("dose", "density")))
  expect_equal(out[["dose:density"]], c(20, 15))
  expect_error(add_management_interactions(mg, list(c("dose", "nope"))),
               "not found")
  expect_error(add_management_interactions(mg, list("dose")),
               "at least two")
})

test_that("interactions flow through into the covariate matrix", {
  sites <- data.frame(environment = c("E1", "E2", "E3"))
  mg <- add_management_interactions(
    data.frame(environment = c("E1", "E2", "E3"),
               fert_type = c("NPK", "urea", "NPK"),
               fert_dose = c(120, 80, 100)),
    list(c("fert_type", "fert_dose")))
  X <- build_enviromic_covariates(sites, management = mg, standardize = FALSE)
  expect_true(any(grepl("fert_type_NPK:fert_dose", colnames(X), fixed = TRUE)))
})

# ---- nest_dose_within_type --------------------------------------------------

test_that("nest_dose_within_type nests dose in type and drops the raw dose", {
  mg <- data.frame(environment = c("E1", "E2", "E3"),
                   fert_type = c("NPK", "urea", "NPK"),
                   fert_dose = c(120, 80, 100))
  out <- nest_dose_within_type(mg, list(fert_dose = "fert_type"))
  expect_equal(out[["fert_type_NPK:fert_dose"]], c(120, 0, 100))
  expect_equal(out[["fert_type_urea:fert_dose"]], c(0, 80, 0))
  expect_false("fert_dose" %in% names(out))        # raw ungrouped dose removed
  expect_true("fert_type" %in% names(out))         # type kept for one-hot
})

test_that("nesting handles multiple domains and missing type as zero", {
  mg <- data.frame(environment = c("E1", "E2", "E3"),
                   fert_type = c("NPK", NA, "NPK"),
                   fert_dose = c(120, 80, 100),
                   irrigation_method = c("drip", "spray", "drip"),
                   water_dose = c(10, 20, 10))
  out <- nest_dose_within_type(mg, list(fert_dose = "fert_type",
                                        water_dose = "irrigation_method"))
  # E2 has an unknown fertiliser type -> its dose is not attributed to any type
  expect_equal(out[["fert_type_NPK:fert_dose"]], c(120, 0, 100))
  expect_equal(out[["irrigation_method_drip:water_dose"]], c(10, 0, 10))
  expect_equal(out[["irrigation_method_spray:water_dose"]], c(0, 20, 0))
  expect_false(any(c("fert_dose", "water_dose") %in% names(out)))
})

test_that("keep_type = FALSE and drop_raw_dose = FALSE are honoured; bad cols error", {
  mg <- data.frame(environment = c("E1", "E2"),
                   fert_type = c("NPK", "urea"), fert_dose = c(120, 80))
  out <- nest_dose_within_type(mg, list(fert_dose = "fert_type"),
                               drop_raw_dose = FALSE, keep_type = FALSE)
  expect_true("fert_dose" %in% names(out))
  expect_false("fert_type" %in% names(out))
  expect_error(nest_dose_within_type(mg, list(bad = "fert_type")), "Dose column")
  expect_error(nest_dose_within_type(mg, list(fert_dose = "nope")), "Type column")
})
