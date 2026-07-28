# Offline unit tests for the SoilGrids retrieval layer: conversion factors,
# request validation, and URL construction for the WCS and WebDAV backends.
# The live network paths are exercised (opt-in) in test-network-fetch.R.

# ---- conversion factors -----------------------------------------------------

test_that("conversion factors match the SoilGrids FAQ table", {
  expect_equal(.soilgrids_conversion_factor("bdod"), 100)
  expect_equal(.soilgrids_conversion_factor("nitrogen"), 100)
  for (p in c("cec", "cfvo", "clay", "sand", "silt", "soc", "ocd", "ocs",
              "phh2o", "wv0010", "wv0033", "wv1500"))
    expect_equal(.soilgrids_conversion_factor(p), 10, info = p)
  # unknown property -> no scaling
  expect_equal(.soilgrids_conversion_factor("not_a_prop"), 1)
})

# ---- request validation -----------------------------------------------------

test_that("invalid depth, quantile or property are rejected", {
  expect_error(.validate_soil_request("7cm", "mean", "clay"), "depth")
  expect_error(.validate_soil_request("0-5cm", "median", "clay"), "quantile")
  expect_error(.validate_soil_request("0-5cm", "mean", "unobtainium"),
               "Unknown SoilGrids")
  expect_true(.validate_soil_request("5-15cm", "Q0.5", c("clay", "sand")))
})

test_that("fetch_soilgrids validates the request before any network work", {
  sites <- data.frame(environment = "E1", latitude = 9, longitude = 7)
  expect_error(fetch_soilgrids(sites, depth = "bad-depth"), "depth")
  expect_error(fetch_soilgrids(sites, quantile = "Q0.99"), "quantile")
  expect_error(fetch_soilgrids(sites, properties = "gold"), "Unknown SoilGrids")
  expect_error(fetch_soilgrids(sites, backend = "nope"))   # match.arg
})

# ---- WCS GetCoverage URL ----------------------------------------------------

test_that("WCS GetCoverage URL is well-formed", {
  u <- .wcs_getcoverage_url("clay", "0-5cm", "mean",
                            bbox = c(-1000, 1000, 2000, 4000))
  expect_match(u, "map=/map/clay.map", fixed = TRUE)
  expect_match(u, "SERVICE=WCS", fixed = TRUE)
  expect_match(u, "VERSION=2.0.1", fixed = TRUE)
  expect_match(u, "REQUEST=GetCoverage", fixed = TRUE)
  expect_match(u, "COVERAGEID=clay_0-5cm_mean", fixed = TRUE)
  expect_match(u, "SUBSET=X(-1000", fixed = TRUE)
  expect_match(u, "SUBSET=Y(2000", fixed = TRUE)
  expect_match(u, "EPSG/0/152160", fixed = TRUE)     # Homolosine
})

test_that("WCS coverage id reflects property, depth and quantile", {
  u <- .wcs_getcoverage_url("nitrogen", "5-15cm", "Q0.5", c(0, 1, 0, 1))
  expect_match(u, "COVERAGEID=nitrogen_5-15cm_Q0.5", fixed = TRUE)
  expect_match(u, "map=/map/nitrogen.map", fixed = TRUE)
})

# ---- WebDAV VRT path --------------------------------------------------------

test_that("WebDAV VRT path uses the /vsicurl root with retries", {
  v <- .webdav_vrt_url("clay", "0-5cm", "mean")
  expect_match(v, "^/vsicurl\\?")
  expect_match(v, "max_retry=3", fixed = TRUE)
  expect_match(v, "list_dir=no", fixed = TRUE)
  expect_match(v, "url=https://files.isric.org/soilgrids/latest/data/",
               fixed = TRUE)
  expect_match(v, "clay/clay_0-5cm_mean.vrt$")
  v2 <- .webdav_vrt_url("soc", "30-60cm", "Q0.95")
  expect_match(v2, "soc/soc_30-60cm_Q0.95.vrt$")
})

# ---- backend dispatch / fallbacks -------------------------------------------

test_that("local backend requires a path for every property (when terra present)", {
  skip_if_not_installed("terra")
  sites <- data.frame(environment = "E1", latitude = 9, longitude = 7)
  expect_error(
    fetch_soilgrids(sites, backend = "local",
                    properties = c("clay", "sand"),
                    local_paths = list(clay = "clay.tif")),
    "local_paths")
})

test_that("build_enviromic_covariates rejects an invalid soil backend", {
  sites <- data.frame(environment = "E1", latitude = 9, longitude = 7)
  expect_error(
    build_enviromic_covariates(sites, fetch_soil = TRUE, soil_backend = "ftp"))
})

# ---- CRS regression (the WCS "output crs is not valid" defect) --------------

test_that("site coordinates project onto the SoilGrids Homolosine grid", {
  skip_if_not_installed("terra")
  sites <- data.frame(environment = c("E1", "E2"),
                      latitude = c(9.06, -1.29), longitude = c(7.49, 36.82))
  pts <- terra::vect(sites, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  crds <- .project_sites_igh(pts)
  expect_equal(nrow(crds), 2L)
  expect_true(all(c("x", "y") %in% names(crds)))
  # projection must yield finite metres, never NA/Inf (the bug produced an
  # invalid output CRS and hence unusable coordinates)
  expect_true(all(is.finite(as.matrix(crds[c("x", "y")]))))
  # a point built in the Homolosine CRS is valid and non-empty
  pt <- .soil_point_igh(crds$x[1], crds$y[1])
  expect_equal(terra::geomtype(pt), "points")
})
