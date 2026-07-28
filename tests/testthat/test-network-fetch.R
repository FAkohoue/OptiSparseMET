# Live integration tests for the weather and soil fetch paths.
#
# Design (per review): a live integration test must not turn network/parser
# failures into skips -- that hides broken integrations behind a green run.
# Instead these tests:
#   1. Never run on CRAN.
#   2. Run only when RUN_LIVE_API_TESTS=true is set explicitly.
#   3. When enabled, FAIL with the original API/network error (and any
#      underlying warnings) rather than skipping.
#
# Run locally with:
#   Sys.setenv(RUN_LIVE_API_TESTS = "true")
#   devtools::test(filter = "network-fetch")
# Disable again with:
#   Sys.unsetenv("RUN_LIVE_API_TESTS")

skip_unless_live_tests_enabled <- function() {
  enabled <- identical(
    tolower(Sys.getenv("RUN_LIVE_API_TESTS", unset = "false")), "true")
  testthat::skip_if_not(
    enabled, "Set RUN_LIVE_API_TESTS=true to run live API integration tests.")
}

# Evaluate a fetch expression, preserving the original error AND any warnings
# (the fetchers report the network cause -- e.g. DNS failure -- as a warning
# before returning NULL, so the top-level error is otherwise generic).
run_fetch <- function(expr) {
  warns <- character(0)
  res <- withCallingHandlers(
    tryCatch(list(value = expr, error = NULL),
             error = function(e) list(value = NULL, error = e)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  res$warnings <- warns
  res
}

format_live_api_error <- function(service, condition, warnings = character(0)) {
  msg <- paste0(
    service, " integration failed.\n",
    "Error class: ", paste(class(condition), collapse = ", "), "\n",
    "Error message: ", conditionMessage(condition))
  if (length(warnings)) {
    msg <- paste0(msg, "\nUnderlying warnings:\n  - ",
                  paste(warnings, collapse = "\n  - "))
  }
  msg
}


test_that("live NASA POWER weather fetch assembles finite covariates", {
  skip_on_cran()
  skip_unless_live_tests_enabled()
  skip_if_not_installed("nasapower")

  sites <- data.frame(
    environment = "E1",
    latitude    = 9.06,
    longitude   = 7.49,
    start_date  = "2020-06-01",
    end_date    = "2020-06-30",
    stringsAsFactors = FALSE)

  res <- run_fetch(build_enviromic_covariates(
    sites, fetch_weather = TRUE, standardize = FALSE))

  if (!is.null(res$error)) {
    fail(format_live_api_error("NASA POWER", res$error, res$warnings))
    return(invisible(NULL))
  }
  X <- res$value
  if (is.null(X)) {
    fail(paste0("NASA POWER fetch returned NULL without erroring. ",
                "Warnings: ", paste(res$warnings, collapse = "; ")))
    return(invisible(NULL))
  }

  expect_true(is.matrix(X) || is.data.frame(X))
  expect_equal(rownames(X), "E1")
  expect_true("mean_temp" %in% colnames(X),
    info = paste("Returned columns:", paste(colnames(X), collapse = ", ")))

  mean_temp <- as.numeric(X["E1", "mean_temp"])
  expect_true(is.finite(mean_temp),
    info = paste("Non-finite mean temperature:", mean_temp))
  expect_gt(mean_temp, 0)
  expect_lt(mean_temp, 50)
})


test_that("live SoilGrids soil fetch (WCS backend) assembles a covariate block", {
  skip_on_cran()
  skip_unless_live_tests_enabled()
  skip_if_not_installed("terra")   # WCS/WebDAV backends read rasters via terra

  sites <- data.frame(
    environment = "E1",
    latitude    = 9.06,
    longitude   = 7.49,
    stringsAsFactors = FALSE)

  # Default soil_backend is now "wcs" (robust OGC Web Coverage Service).
  res <- run_fetch(build_enviromic_covariates(
    sites, fetch_soil = TRUE, soil_backend = "wcs", standardize = FALSE))

  # Distinguish the failure modes. Only ONE is not our integration's fault: the
  # service was reached and read successfully but returned only null/NA values
  # for the point (a known, intermittent ISRIC outage, or a point off
  # coverage). That -- and only that -- is a legitimate skip. A network/read
  # error (DNS, connection, HTTP) surfaces as a "Soil fetch failed ..." warning
  # and must FAIL, even though it also ends in an all-NA block.
  fetch_error <- any(grepl("Soil fetch failed|resolve|cannot open|connection|HTTP|curl|status",
                           res$warnings, ignore.case = TRUE))
  only_null   <- any(grepl("only null", res$warnings))
  if (only_null && !fetch_error) {
    skip(paste0("SoilGrids was reachable but returned only null/NA values ",
                "(upstream ISRIC service degradation). The fetch/parse path ",
                "is exercised; only real data is absent."))
  }

  if (!is.null(res$error)) {
    fail(format_live_api_error("SoilGrids", res$error, res$warnings))
    return(invisible(NULL))
  }
  X <- res$value

  # A terrestrial point returning NULL is an integration failure, not a skip.
  # (SoilGrids intermittently answers valid points with all-null values; when
  # that happens the warnings above name it as a service degradation.)
  if (is.null(X)) {
    fail(paste0("SoilGrids returned NULL for lat = 9.06, lon = 7.49. ",
                "Warnings: ", paste(res$warnings, collapse = "; ")))
    return(invisible(NULL))
  }
  expect_true(is.matrix(X) || is.data.frame(X))
  expect_equal(rownames(X), "E1")
  expect_gte(ncol(X), 1L)

  soil_values <- as.numeric(X["E1", , drop = TRUE])
  expect_true(all(is.finite(soil_values) | is.na(soil_values)))
  expect_true(any(is.finite(soil_values)),
    info = "All returned SoilGrids covariates are missing.")

  # Clay is a fraction in conventional units (g/100g = %); a sane value is well
  # within [0, 100], confirming the conversion factor was applied.
  if ("clay" %in% colnames(X) && is.finite(X["E1", "clay"])) {
    expect_gte(X["E1", "clay"], 0)
    expect_lte(X["E1", "clay"], 100)
  }
})


test_that("live SoilGrids WebDAV backend returns finite clay data", {
  skip_on_cran()
  skip_unless_live_tests_enabled()
  skip_if_not_installed("terra")
  skip_if_not_installed("gdalUtilities")   # WebDAV crop uses gdal_translate

  # Several nearby rural coordinates: a single point can fall on an urban or
  # otherwise masked SoilGrids cell.
  sites <- data.frame(
    environment = c("E1", "E2", "E3"),
    latitude    = c(9.30, 9.50, 8.80),
    longitude   = c(7.10, 7.00, 7.20),
    stringsAsFactors = FALSE)

  res <- run_fetch(fetch_soilgrids(
    sites, backend = "webdav", properties = "clay",
    depth = "0-5cm", quantile = "mean"))

  fetch_error <- any(grepl(
    paste(c("Soil fetch failed", "crop failed", "could not create",
            "could not be read", "resolve", "cannot open", "connection",
            "HTTP", "curl", "status", "GDAL"), collapse = "|"),
    res$warnings, ignore.case = TRUE))
  only_null <- any(grepl("only null", res$warnings))
  if (only_null && !fetch_error) {
    skip(paste("SoilGrids WebDAV was reachable but returned only null values",
               "for all test locations."))
  }
  if (!is.null(res$error)) {
    fail(format_live_api_error("SoilGrids WebDAV", res$error, res$warnings))
    return(invisible(NULL))
  }
  X <- res$value
  if (is.null(X)) {
    fail(paste0("WebDAV fetch returned NULL. Warnings: ",
                paste(res$warnings, collapse = "; ")))
    return(invisible(NULL))
  }
  expect_true("clay" %in% colnames(X))
  expect_true(any(is.finite(X[, "clay"])),
    info = "No finite clay value returned for any WebDAV test location.")
})
