#' Assemble enviromic covariates for sites (weather, soil, management)
#'
#' @description
#' Turns site information into the environment-by-covariate matrix that
#' [build_environment_relationship()] needs. It merges three sources:
#' \describe{
#'   \item{Weather}{Optionally fetched from site GPS coordinates over a growing
#'     window (via the \pkg{nasapower} POWER API) and summarised to per-site
#'     features, or supplied directly.}
#'   \item{Soil}{Optionally fetched from coordinates via ISRIC SoilGrids using a
#'     robust backend (WCS or WebDAV; see [fetch_soilgrids()]), or supplied
#'     directly.}
#'   \item{Management}{User-supplied per-site management such as planting system,
#'     fertilisation, and irrigation. Numeric columns are kept; categorical
#'     columns are one-hot encoded.}
#' }
#'
#' Fetching requires internet access and the relevant optional packages; it is
#' never triggered by examples or tests. When a source cannot be fetched the
#' function warns and continues with whatever is available.
#'
#' @param sites Data frame with columns `environment`, `latitude`, `longitude`,
#'   and (for weather fetching) `start_date`, `end_date` (the growing window).
#' @param management Optional data frame with an `environment` column plus
#'   management covariates (numeric or categorical). Categorical columns are
#'   one-hot encoded. Management **doses** (fertiliser, pesticide, irrigation
#'   amount) should be nested within their type first with
#'   [nest_dose_within_type()] -- a raw dose is not comparable across products
#'   (100 kg NPK is not 100 kg urea). Use [add_management_interactions()] for
#'   other interactions.
#' @param weather,soil Optional pre-fetched matrices/data frames (one row per
#'   environment) to use instead of fetching.
#' @param fetch_weather,fetch_soil Logical; fetch from coordinates if `TRUE`.
#'   Fetched variables are **merged** with any `weather`/`soil` you supply -- your
#'   own data takes precedence and only the columns you are missing are fetched,
#'   so you can complete a partly-collected dataset automatically.
#' @param weather_pars Character vector of NASA POWER parameter codes to fetch
#'   (e.g. `c("T2M", "PRECTOTCORR", "WS2M", "GWETROOT", "EVPTRNS")`). Defaults to
#'   the six core variables. See [enviromic_variable_catalog()] for codes and
#'   meanings.
#' @param weather_stats Optional named vector giving how to aggregate each
#'   requested parameter over the window: `"mean"` (default), `"sum"` (e.g.
#'   precipitation, evapotranspiration), `"min"`, `"max"`, `"sd"`, or `"median"`.
#' @param weather_temporal NASA POWER temporal resolution to pull before
#'   aggregating: `"daily"` (default), `"monthly"`, or `"climatology"`
#'   (long-term). The covariate matrix always holds one aggregated value per
#'   site; for stage-specific covariates, set each site's `start_date`/`end_date`
#'   to the stage of interest and merge the results.
#' @param soil_backend Backend used when `fetch_soil = TRUE`: `"wcs"` (default),
#'   `"webdav"`, `"rest"`, or `"local"`. See [fetch_soilgrids()].
#' @param soil_properties,soil_depth,soil_quantile SoilGrids property names,
#'   depth interval, and prediction statistic to retrieve (see
#'   [fetch_soilgrids()]).
#' @param soil_local_paths Named list mapping properties to local raster files,
#'   required when `soil_backend = "local"`.
#' @param standardize Logical; scale the assembled numeric covariates. Default
#'   `TRUE`.
#' @param scale_dummies Logical; when `standardize = TRUE`, also centre/scale the
#'   0/1 dummy columns (one-hot management, indicators). Default `TRUE`; set
#'   `FALSE` to leave dummies as 0/1 while still scaling the continuous
#'   weather/soil/dose covariates.
#' @param min_coverage Minimum acceptable observed fraction per source-variable.
#' @param missing_action `"impute"` (default), `"warn"`, or `"error"`.
#' @param impute Explicit missing-data rule: `"median"` or `"none"`. Imputation
#'   and coverage are returned as matrix attributes rather than being silent.
#'
#' @return An environment-by-covariate numeric matrix (row names = environment),
#'   ready for [build_environment_relationship()] with `source = "enviromic"`.
#'   Attributes `environment_audit`, `environment_provenance`,
#'   `environment_imputation`, and `environment_blocks` preserve the source
#'   quality-control record.
#'
#' @seealso [build_environment_relationship()], [cluster_environments()].
#' @examples
#' mgmt <- data.frame(environment = c("E1", "E2"),
#'                    planting_system = c("conventional", "no_till"),
#'                    fertN = c(120, 80))
#' wx <- data.frame(environment = c("E1", "E2"),
#'                  mean_temp = c(21.5, 24.0), total_precip = c(410, 300))
#' X <- build_enviromic_covariates(
#'   sites = data.frame(environment = c("E1", "E2"),
#'                      latitude = c(-1.1, 9.0), longitude = c(37.0, 7.5)),
#'   management = mgmt, weather = wx)
#' \dontrun{
#' # Fetch weather and soil directly from the coordinates (needs internet):
#' X <- build_enviromic_covariates(sites, management = mgmt,
#'                                 fetch_weather = TRUE, fetch_soil = TRUE)
#' }
#' @export
build_enviromic_covariates <- function(sites, management = NULL,
                                       weather = NULL, soil = NULL,
                                       fetch_weather = FALSE, fetch_soil = FALSE,
                                       weather_pars = NULL, weather_stats = NULL,
                                       weather_temporal = "daily",
                                       soil_backend = c("wcs", "webdav",
                                                        "rest", "local"),
                                       soil_properties = c("bdod", "cec",
                                                           "cfvo", "clay",
                                                           "nitrogen", "ocd",
                                                           "phh2o", "sand",
                                                           "silt", "soc"),
                                       soil_depth = "0-5cm",
                                       soil_quantile = "mean",
                                       soil_local_paths = NULL,
                                       standardize = TRUE,
                                       scale_dummies = TRUE,
                                       min_coverage = 0.80,
                                       missing_action = c("impute", "warn",
                                                          "error"),
                                       impute = c("median", "none")) {
  soil_backend <- match.arg(soil_backend)
  missing_action <- match.arg(missing_action)
  impute <- match.arg(impute)
  sites <- as.data.frame(sites)
  if (!"environment" %in% names(sites))
    stop("`sites` must have an `environment` column.")
  envs <- as.character(sites$environment)
  if (anyNA(envs) || any(!nzchar(envs)) || anyDuplicated(envs))
    stop("`sites$environment` must be unique, non-missing, and non-empty.")

  by_env <- function(df) {
    df <- as.data.frame(df)
    if (!"environment" %in% names(df)) stop("source data must have an `environment` column.")
    rn <- as.character(df$environment)
    if (anyNA(rn) || any(!nzchar(rn)) || anyDuplicated(rn))
      stop("Every source must have unique, non-missing environment rows.")
    extra <- setdiff(rn, envs)
    if (length(extra))
      stop("Source contains environments not present in `sites`: ",
           paste(extra, collapse = ", "), ".")
    m  <- df[match(envs, rn), setdiff(names(df), "environment"), drop = FALSE]
    rownames(m) <- envs
    m
  }

  blocks <- list()

  ## weather -- merge your own data with fetched variables (yours takes
  ## precedence; only the columns you do not already have are added).
  w_sup <- if (!is.null(weather)) by_env(weather) else NULL
  w_fet <- if (fetch_weather) tryCatch(
    .fetch_weather_power(sites, pars = weather_pars, stats = weather_stats,
                         temporal = weather_temporal),
    error = function(e) {
      warning("Weather source failed independently: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }) else NULL
  w <- .merge_cov_blocks(w_sup, w_fet)
  if (!is.null(w)) blocks$weather <- w

  ## soil -- same merge logic
  s_sup <- if (!is.null(soil)) by_env(soil) else NULL
  s_fet <- if (fetch_soil) tryCatch(
    fetch_soilgrids(sites, backend = soil_backend, properties = soil_properties,
                    depth = soil_depth, quantile = soil_quantile,
                    local_paths = soil_local_paths),
    error = function(e) {
      warning("Soil source failed independently: ", conditionMessage(e),
              call. = FALSE)
      NULL
    }) else NULL
  s <- .merge_cov_blocks(s_sup, s_fet)
  if (!is.null(s)) blocks$soil <- s

  ## management (one-hot for categorical)
  if (!is.null(management)) blocks$management <- .onehot(by_env(management))

  if (!length(blocks))
    stop("No covariate sources available. Supply/fetch weather, soil, or management.")

  qc <- lapply(names(blocks), function(nm) {
    b <- as.data.frame(blocks[[nm]], check.names = FALSE,
                       stringsAsFactors = FALSE)
    b$environment <- rownames(blocks[[nm]])
    b <- b[, c("environment", setdiff(names(b), "environment")), drop = FALSE]
    qc_environmental_data(
      b, environments = envs, min_coverage = min_coverage,
      missing_action = missing_action, impute = impute,
      add_missing_indicators = TRUE, source = nm
    )
  })
  names(qc) <- names(blocks)
  clean_blocks <- lapply(qc, function(q)
    q$data[, setdiff(names(q$data), "environment"), drop = FALSE])
  X <- do.call(cbind, lapply(clean_blocks, data.matrix))
  rownames(X) <- envs
  if (standardize) {
    # Optionally leave 0/1 dummy columns (one-hot management, indicators)
    # untouched rather than centring/scaling them.
    scale_cols <- rep(TRUE, ncol(X))
    if (!scale_dummies) {
      is_binary <- apply(X, 2L, function(col) {
        u <- unique(col[!is.na(col)])
        length(u) > 0L && all(u %in% c(0, 1))
      })
      scale_cols <- !is_binary
    }
    if (any(scale_cols)) {
      for (j in which(scale_cols)) {
        observed <- is.finite(X[, j])
        mu <- if (any(observed)) mean(X[observed, j]) else NA_real_
        ss <- if (sum(observed) > 1L) stats::sd(X[observed, j]) else NA_real_
        if (is.finite(ss) && ss > 0) {
          X[observed, j] <- (X[observed, j] - mu) / ss
        } else {
          X[observed, j] <- 0
        }
      }
    }
    attr(X, "scaled:center") <- attr(X, "scaled:scale") <- NULL
  }
  X <- as.matrix(X)
  attr(X, "environment_audit") <- .bind_qc_component(qc, "audit")
  attr(X, "environment_provenance") <- .bind_qc_component(qc, "provenance")
  attr(X, "environment_imputation") <- .bind_qc_component(qc, "imputation")
  attr(X, "environment_blocks") <- stats::setNames(
    vapply(clean_blocks, ncol, integer(1)), names(clean_blocks)
  )
  X
}


# One-hot encode categorical columns; keep numeric columns.
.onehot <- function(df) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  out <- list()
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.numeric(col)) {
      out[[nm]] <- col
    } else {
      f <- factor(col)
      for (lv in levels(f))
        out[[paste0(nm, "_", lv)]] <- as.integer(f == lv)
    }
  }
  # check.names = FALSE keeps interaction column names such as "a:b" intact.
  m <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(m) <- rownames(df)
  m
}


# Friendly names for the default NASA POWER parameters (kept for continuity).
.POWER_FRIENDLY <- c(T2M = "mean_temp", T2M_MAX = "max_temp", T2M_MIN = "min_temp",
                     PRECTOTCORR = "total_precip",
                     ALLSKY_SFC_SW_DWN = "mean_radiation", RH2M = "mean_humidity")
.POWER_DEFAULT_STAT <- c(T2M = "mean", T2M_MAX = "mean", T2M_MIN = "mean",
                         PRECTOTCORR = "sum", ALLSKY_SFC_SW_DWN = "mean",
                         RH2M = "mean")

# Column name a parameter maps to (friendly for defaults, else the POWER code).
.power_colname <- function(p) if (p %in% names(.POWER_FRIENDLY))
  unname(.POWER_FRIENDLY[p]) else p

# Aggregate a daily/monthly POWER data frame to one value per parameter.
.aggregate_weather <- function(d, stat_map) {
  vals <- lapply(names(stat_map), function(p) {
    if (!p %in% names(d)) return(NA_real_)
    x <- d[[p]]; x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    switch(stat_map[[p]], mean = mean(x), sum = sum(x), min = min(x),
           max = max(x), sd = stats::sd(x), median = stats::median(x), mean(x))
  })
  out <- as.data.frame(stats::setNames(vals, vapply(names(stat_map), .power_colname,
                                                    character(1))),
                       stringsAsFactors = FALSE, check.names = FALSE)
  out
}

# Fetch and summarise weather over each site's window from NASA POWER. `pars`
# lets any POWER parameter be requested; `stats` sets how each is aggregated
# (mean/sum/min/max/sd/median); `temporal` is "daily", "monthly" or
# "climatology".
.fetch_weather_power <- function(sites, pars = NULL, stats = NULL,
                                 temporal = "daily") {
  if (!requireNamespace("nasapower", quietly = TRUE)) {
    warning("Package 'nasapower' not installed; skipping weather fetch.")
    return(NULL)
  }
  need <- if (temporal == "climatology") c("latitude", "longitude")
          else c("latitude", "longitude", "start_date", "end_date")
  if (!all(need %in% names(sites))) {
    warning("`sites` needs ", paste(need, collapse = ", "), " for weather fetch.")
    return(NULL)
  }

  if (is.null(pars)) {
    pars_vec <- names(.POWER_DEFAULT_STAT); stat_map <- as.list(.POWER_DEFAULT_STAT)
  } else {
    pars_vec <- as.character(pars)
    stat_map <- stats::setNames(rep(list("mean"), length(pars_vec)), pars_vec)
    if (!is.null(stats)) for (p in names(stats))
      if (p %in% pars_vec) stat_map[[p]] <- stats[[p]]
  }
  out_names <- vapply(pars_vec, .power_colname, character(1))
  na_row <- as.data.frame(stats::setNames(as.list(rep(NA_real_, length(out_names))),
                                          out_names), check.names = FALSE)

  rows <- lapply(seq_len(nrow(sites)), function(i) {
    dates <- if (temporal == "climatology") NULL else
      c(as.character(sites$start_date[i]), as.character(sites$end_date[i]))
    d <- .power_get_retry(
      community = "ag", pars = pars_vec, temporal_api = temporal,
      lonlat = c(sites$longitude[i], sites$latitude[i]), dates = dates,
      label = sites$environment[i])
    if (is.null(d)) return(NULL)
    .aggregate_weather(as.data.frame(d), stat_map)
  })
  if (all(vapply(rows, is.null, logical(1)))) return(NULL)
  out <- do.call(rbind, lapply(rows, function(r) if (is.null(r)) na_row else r))
  rownames(out) <- as.character(sites$environment)
  out
}

# Merge two per-environment covariate blocks cell by cell. Supplied non-missing
# values win; fetched values fill both missing columns and missing cells.
.merge_cov_blocks <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  a <- as.data.frame(a, check.names = FALSE)
  b <- as.data.frame(b, check.names = FALSE)
  envs <- union(rownames(a), rownames(b))
  a <- a[match(envs, rownames(a)), , drop = FALSE]
  b <- b[match(envs, rownames(b)), , drop = FALSE]
  rownames(a) <- rownames(b) <- envs
  shared <- intersect(names(a), names(b))
  for (nm in shared) {
    fill <- is.na(a[[nm]]) & !is.na(b[[nm]])
    a[[nm]][fill] <- b[[nm]][fill]
  }
  extra <- setdiff(colnames(b), colnames(a))
  if (length(extra))
    a <- cbind(a, b[, extra, drop = FALSE])
  a
}


# Fetch topsoil properties per site from the ISRIC SoilGrids REST API.
.fetch_soil_soilgrids <- function(sites) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("Package 'jsonlite' not installed; skipping soil fetch.")
    return(NULL)
  }
  props <- c("clay", "sand", "silt", "phh2o", "soc", "cec")
  q <- paste0("&property=", props, collapse = "")
  rows <- lapply(seq_len(nrow(sites)), function(i) {
    url <- sprintf(
      "https://rest.isric.org/soilgrids/v2.0/properties/query?lon=%f&lat=%f%s&depth=0-5cm&value=mean",
      sites$longitude[i], sites$latitude[i], q)
    js <- tryCatch(jsonlite::fromJSON(url), error = function(e) {
      warning(sprintf("SoilGrids request failed for '%s': %s",
                      sites$environment[i], conditionMessage(e)), call. = FALSE)
      NULL
    })
    if (is.null(js)) return(NULL)
    vals <- tryCatch({
      layers <- js$properties$layers
      raw <- vapply(seq_len(nrow(layers)), function(r)
        layers$depths[[r]]$values$mean[1], numeric(1))
      # SoilGrids reports integer "mapped units"; divide by the per-property
      # d_factor to convert to physical target units (e.g. clay g/kg -> %).
      dfac <- tryCatch(layers$unit_measure$d_factor, error = function(e) NULL)
      if (is.null(dfac) || length(dfac) != length(raw)) dfac <- rep(1, length(raw))
      dfac[!is.finite(dfac) | dfac == 0] <- 1
      stats::setNames(raw / dfac, layers$name)
    }, error = function(e) NULL)
    if (is.null(vals)) NULL else as.data.frame(as.list(vals))
  })
  if (all(vapply(rows, is.null, logical(1)))) {
    warning("SoilGrids fetch returned nothing; skipping soil.", call. = FALSE)
    return(NULL)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- as.character(sites$environment)[seq_len(nrow(out))]
  # SoilGrids intermittently answers with all-null values even on valid points;
  # an all-NA block is useless as a covariate, so treat it as no data.
  if (all(vapply(out, function(col) all(is.na(col)), logical(1)))) {
    warning("SoilGrids returned only null values (service may be degraded); ",
            "skipping soil.", call. = FALSE)
    return(NULL)
  }
  out
}
