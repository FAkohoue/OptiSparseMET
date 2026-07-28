#' Temporal envirotyping features from a daily weather series
#'
#' @description
#' Collapsing a season to one mean per site discards the within-season variation,
#' extremes, and *timing* that often drive genotype-by-environment interaction.
#' `envirotype_features()` keeps that information: from a daily (or sub-daily)
#' weather series per environment it builds a rich, fixed-length feature vector
#' per environment -- by time window / crop stage, with variability and stress
#' statistics, growing degree days, and an optional functional (spline)
#' representation of the whole trajectory. The resulting environment-by-feature
#' matrix plugs straight into [build_environment_relationship()] (or
#' [build_enviromic_covariates()] as `weather =`), so the environmental kinship
#' reflects how sites differ in their temporal dynamics, not just their means.
#'
#' @details
#' All features are computed on *phenological time* -- each series is indexed by
#' day-after-start, so sites with different season lengths are compared on a
#' common footing (this is why dynamic time warping is not used: warping would
#' erase the timing signal we want to keep).
#' \describe{
#'   \item{Windows / stages}{`windows` may be a single integer `K` (split each
#'     season into `K` equal intervals -- no crop knowledge needed), a numeric
#'     vector of day-after-start cut points (e.g. `c(30, 60, 90)`), or a named
#'     list of day ranges from crop physiology
#'     (e.g. `list(vegetative = c(0, 40), flowering = c(40, 70),
#'     grain_fill = c(70, 120))`). Expert stages give the most interpretable
#'     features; equal intervals or thermal-time cut points are the
#'     knowledge-free alternatives.}
#'   \item{Statistics}{Each `stats` (e.g. `"mean"`, `"sd"`, `"cv"`, `"sum"`,
#'     `"min"`, `"max"`, `"median"`) is computed per window per parameter; `"sd"`
#'     and `"cv"` capture the day-to-day variability the mean hides.}
#'   \item{Stress days}{`thresholds` counts days meeting a condition per window,
#'     e.g. `list(heat_days = list(par = "T2M_MAX", above = 35))`.}
#'   \item{Growing degree days}{`gdd` accumulates degree days per window, e.g.
#'     `list(gdd = list(par = "T2M", base = 10, upper = 30))`.}
#'   \item{Functional}{`functional` fits a B-spline basis over normalised season
#'     time and returns the coefficients as shape features, e.g.
#'     `list(temp_shape = list(par = "T2M", df = 4))`.}
#' }
#'
#' @param daily A data frame with an environment column, a day column (numeric
#'   day-after-start; if absent it is derived from row order within each
#'   environment), and one column per weather parameter. Long format (many rows
#'   per environment).
#' @param windows Integer `K`, numeric cut points, or a named list of day ranges
#'   (see Details). Default 1 (whole season).
#' @param stats Character vector of per-window statistics. Default `"mean"`.
#' @param thresholds Named list of stress-day definitions (see Details).
#' @param gdd Named list of growing-degree-day definitions (see Details).
#' @param functional Named list of functional (spline) specifications (see
#'   Details).
#' @param stress Named list of cumulative stress-severity definitions. Each uses
#'   `par` and either `above` or `below`; severity is the accumulated distance
#'   beyond the threshold rather than only a count of stress days.
#' @param spells Named list of consecutive-event definitions using `par` and
#'   `above` or `below`. Both the longest run and number of runs are returned.
#' @param phenology Optional data frame with one row per environment-stage and
#'   columns `environment`, `stage`, plus either `start_day`/`end_day` or
#'   `start_date`/`end_date`. These observed or breeder-specified stages replace
#'   the common `windows` for the corresponding environment.
#' @param env_col,day_col,date_col Column names for environment id, day index,
#'   and calendar date. When day is absent but date is available, the day index
#'   is computed from actual dates after sorting. When both are unavailable,
#'   row order is used within environment and one consolidated warning names
#'   every affected environment.
#' @return An environment-by-feature numeric matrix (row names = environments).
#' @references
#' Costa-Neto, G., Fritsche-Neto, R., & Crossa, J. (2021). Nonlinear kernels,
#' dominance, and envirotyping data increase prediction accuracy... *Heredity*.
#' Jarquin, D. et al. (2014). A reaction norm model for genomic selection using
#' high-dimensional genomic and environmental data. *TAG* 127, 595-607.
#' @seealso [fetch_weather_series()], [build_environment_relationship()],
#'   [build_enviromic_covariates()].
#' @examples
#' set.seed(1)
#' mk <- function(e, len, base) data.frame(
#'   environment = e, day = 0:(len - 1),
#'   T2M = base + sin((0:(len - 1)) / 10) * 3 + rnorm(len, 0, 0.5),
#'   T2M_MAX = base + 6 + rnorm(len, 0, 0.5),
#'   PRECTOTCORR = pmax(0, rnorm(len, 3, 3)))
#' daily <- rbind(mk("E1", 90, 24), mk("E2", 80, 28))
#' W <- envirotype_features(
#'   daily, windows = list(veg = c(0, 40), flower = c(40, 70), fill = c(70, 120)),
#'   stats = c("mean", "sd"),
#'   thresholds = list(heat_days = list(par = "T2M_MAX", above = 33)),
#'   gdd = list(gdd = list(par = "T2M", base = 10)),
#'   functional = list(temp_shape = list(par = "T2M", df = 4)))
#' dim(W)
#' @export
envirotype_features <- function(daily, windows = 1, stats = "mean",
                                thresholds = NULL, gdd = NULL,
                                functional = NULL, stress = NULL, spells = NULL,
                                phenology = NULL,
                                env_col = "environment", day_col = "day",
                                date_col = "date") {
  df <- as.data.frame(daily, stringsAsFactors = FALSE)
  if (!env_col %in% names(df)) stop("`daily` needs an `", env_col, "` column.")
  if (!day_col %in% names(df)) df[[day_col]] <- NA_real_
  if (date_col %in% names(df)) {
    df[[date_col]] <- as.Date(df[[date_col]])
    if (anyNA(df[[date_col]]))
      stop("`", date_col, "` contains invalid or missing dates.")
  }

  envs <- unique(as.character(df[[env_col]]))
  row_order_envs <- if (date_col %in% names(df)) {
    character(0)
  } else {
    envs[vapply(envs, function(e) {
      day <- suppressWarnings(as.numeric(
        df[df[[env_col]] == e, day_col, drop = TRUE]
      ))
      all(is.na(day))
    }, logical(1))]
  }
  if (length(row_order_envs)) {
    warning(
      "No day or date information for environment",
      if (length(row_order_envs) == 1L) " " else "s ",
      paste0("'", row_order_envs, "'", collapse = ", "),
      "; deriving day from row order within each environment.",
      call. = FALSE
    )
  }
  metadata <- c(env_col, day_col, date_col, "year", "YEAR", "MO", "DY",
                "month", "day_of_month", "cache_hit", "weather_source")
  par_cols <- setdiff(names(df), metadata)
  par_cols <- par_cols[vapply(df[par_cols], is.numeric, logical(1))]
  if (!length(par_cols)) stop("No numeric weather parameters found in `daily`.")
  specs <- c(thresholds, gdd, functional, stress, spells)
  requested <- unique(vapply(specs, function(z)
    if (is.null(z$par)) NA_character_ else as.character(z$par), character(1)))
  requested <- requested[!is.na(requested)]
  absent <- setdiff(requested, par_cols)
  if (length(absent))
    stop("Feature definitions reference missing weather parameters: ",
         paste(absent, collapse = ", "), ".")

  rows <- lapply(envs, function(e) {
    de <- df[df[[env_col]] == e, , drop = FALSE]
    if (date_col %in% names(de))
      de <- de[order(de[[date_col]]), , drop = FALSE]
    day <- suppressWarnings(as.numeric(de[[day_col]]))
    if (all(is.na(day)) && date_col %in% names(de))
      day <- as.integer(de[[date_col]] - min(de[[date_col]]))
    if (all(is.na(day)))
      day <- seq_len(nrow(de)) - 1L
    if (anyNA(day) || anyDuplicated(day))
      stop("Day indices must be complete and unique within environment '", e, "'.")
    de[[day_col]] <- day
    wins <- .environment_feature_windows(
      phenology, e, windows, day,
      if (date_col %in% names(de)) de[[date_col]] else NULL
    )

    r <- list()
    for (w in wins) {
      sel <- if (w$last) day >= w$lo & day <= w$hi else day >= w$lo & day < w$hi
      sub <- de[sel, , drop = FALSE]
      for (p in par_cols) for (st in stats)
        r[[paste0(p, "_", st, "_", w$name)]] <- .agg_stat(sub[[p]], st)
      for (tn in names(thresholds))
        r[[paste0(tn, "_", w$name)]] <-
          .count_threshold(sub[[thresholds[[tn]]$par]], thresholds[[tn]])
      for (gn in names(gdd))
        r[[paste0(gn, "_", w$name)]] <- .gdd_sum(sub[[gdd[[gn]]$par]], gdd[[gn]])
      for (sn in names(stress))
        r[[paste0(sn, "_severity_", w$name)]] <-
          .stress_severity(sub[[stress[[sn]]$par]], stress[[sn]])
      for (rn in names(spells)) {
        sp <- .spell_summary(sub[[spells[[rn]]$par]], spells[[rn]])
        r[[paste0(rn, "_longest_", w$name)]] <- sp[["longest"]]
        r[[paste0(rn, "_runs_", w$name)]] <- sp[["runs"]]
      }
    }
    for (fn in names(functional)) {
      sp <- functional[[fn]]
      k  <- if (is.null(sp$df)) 4L else as.integer(sp$df)
      coefs <- .functional_coef(day, de[[sp$par]], k)
      for (j in seq_along(coefs)) r[[paste0(fn, "_f", j)]] <- coefs[j]
    }
    r
  })

  allnames <- unique(unlist(lapply(rows, names)))
  M <- matrix(NA_real_, length(envs), length(allnames),
              dimnames = list(envs, allnames))
  for (i in seq_along(envs)) {
    ri <- rows[[i]]
    if (length(ri)) M[i, names(ri)] <- unlist(ri, use.names = FALSE)
  }
  M
}


#' Rice-focused weather and crop-stress envirotyping
#'
#' @description
#' Adds physiologically interpretable rice weather variables before calling
#' [envirotype_features()]: vapour-pressure deficit, diurnal temperature range,
#' precipitation-minus-evapotranspiration water balance, heat/cold severity,
#' hot nights, dry/hot/wet spells, root-zone moisture deficit, solar radiation,
#' and stage-specific summaries. Stages may come from observed phenology,
#' breeder-specified calendar windows, or thermal-time windows.
#'
#' @param daily Long daily weather data with NASA POWER-style names.
#' @param windows Stage windows used when `phenology` is absent.
#' @param phenology Optional observed environment-stage table; see
#'   [envirotype_features()].
#' @param time_scale `"calendar"` or `"thermal"`. Thermal time accumulates GDD
#'   from `base_temperature` and interprets `windows` on that scale.
#' @param base_temperature Base temperature for rice thermal time.
#' @param stats Per-stage summary statistics.
#' @param heat_threshold,hot_night_threshold,cold_threshold Temperature stress
#'   thresholds in degrees C.
#' @param dry_day_threshold,wet_day_threshold Rain thresholds in mm/day.
#' @param root_moisture_threshold Threshold for `GWETROOT` when present.
#' @param env_col,day_col,date_col Column names.
#'
#' @return An environment-by-feature numeric matrix. Attributes include the
#'   derived-variable definitions.
#' @export
rice_weather_features <- function(
    daily,
    windows = list(vegetative = c(0, 40), reproductive = c(40, 75),
                   grain_fill = c(75, 120)),
    phenology = NULL, time_scale = c("calendar", "thermal"),
    base_temperature = 10, stats = c("mean", "sd", "min", "max", "sum"),
    heat_threshold = 35, hot_night_threshold = 24, cold_threshold = 15,
    dry_day_threshold = 1, wet_day_threshold = 20,
    root_moisture_threshold = 0.30,
    env_col = "environment", day_col = "day", date_col = "date") {
  time_scale <- match.arg(time_scale)
  d <- as.data.frame(daily, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c(env_col, "T2M") %in% names(d)))
    stop("`daily` needs environment and T2M columns.")
  if (all(c("T2M", "RH2M") %in% names(d))) {
    es <- 0.6108 * exp(17.27 * d$T2M / (d$T2M + 237.3))
    d$VPD <- pmax(0, es * (1 - pmin(100, pmax(0, d$RH2M)) / 100))
  }
  if (all(c("T2M_MAX", "T2M_MIN") %in% names(d)))
    d$DIURNAL_TEMP_RANGE <- d$T2M_MAX - d$T2M_MIN
  if ("PRECTOTCORR" %in% names(d)) {
    et_name <- intersect(c("EVPTRNS", "ET0"), names(d))
    if (length(et_name))
      d$WATER_BALANCE <- d$PRECTOTCORR - d[[et_name[1L]]]
  }

  use_day <- day_col
  if (time_scale == "thermal") {
    d$thermal_time <- NA_real_
    for (e in unique(as.character(d[[env_col]]))) {
      ii <- which(d[[env_col]] == e)
      if (date_col %in% names(d))
        ii <- ii[order(as.Date(d[[date_col]][ii]))]
      else if (day_col %in% names(d))
        ii <- ii[order(as.numeric(d[[day_col]][ii]))]
      d$thermal_time[ii] <- cumsum(pmax(0, d$T2M[ii] - base_temperature))
    }
    use_day <- "thermal_time"
    if (!is.null(phenology))
      warning("`phenology` takes precedence over thermal-time windows.",
              call. = FALSE)
  }

  thresholds <- list()
  stress <- list()
  spells <- list()
  if ("T2M_MAX" %in% names(d)) {
    thresholds$heat_days <- list(par = "T2M_MAX", above = heat_threshold)
    stress$heat <- list(par = "T2M_MAX", above = heat_threshold)
    spells$hot_spell <- list(par = "T2M_MAX", above = heat_threshold)
  }
  if ("T2M_MIN" %in% names(d)) {
    thresholds$hot_nights <- list(par = "T2M_MIN",
                                  above = hot_night_threshold)
    thresholds$cold_nights <- list(par = "T2M_MIN", below = cold_threshold)
    stress$hot_night <- list(par = "T2M_MIN",
                             above = hot_night_threshold)
    stress$cold <- list(par = "T2M_MIN", below = cold_threshold)
  }
  if ("PRECTOTCORR" %in% names(d)) {
    thresholds$dry_days <- list(par = "PRECTOTCORR",
                                below = dry_day_threshold)
    thresholds$heavy_rain_days <- list(par = "PRECTOTCORR",
                                       above = wet_day_threshold)
    spells$dry_spell <- list(par = "PRECTOTCORR",
                             below = dry_day_threshold)
    spells$wet_spell <- list(par = "PRECTOTCORR",
                             above = wet_day_threshold)
  }
  if ("GWETROOT" %in% names(d)) {
    thresholds$root_water_deficit_days <- list(
      par = "GWETROOT", below = root_moisture_threshold
    )
    stress$root_water_deficit <- list(
      par = "GWETROOT", below = root_moisture_threshold
    )
  }
  out <- envirotype_features(
    d, windows = windows, stats = stats, thresholds = thresholds,
    gdd = list(gdd = list(par = "T2M", base = base_temperature,
                          upper = 42)),
    stress = stress, spells = spells, phenology = phenology,
    env_col = env_col, day_col = use_day, date_col = date_col
  )
  attr(out, "derived_variables") <- c(
    VPD = "saturation vapour pressure x (1 - RH/100), kPa",
    DIURNAL_TEMP_RANGE = "T2M_MAX - T2M_MIN, degrees C",
    WATER_BALANCE = "PRECTOTCORR - evapotranspiration, mm/day"
  )[c("VPD", "DIURNAL_TEMP_RANGE", "WATER_BALANCE") %in% names(d)]
  out
}


#' Fetch a daily weather series from NASA POWER
#'
#' @description
#' Returns the per-site **daily** weather series (long format) rather than a
#' single aggregated value, so it can be turned into rich temporal features with
#' [envirotype_features()]. Requires internet and the \pkg{nasapower} package.
#'
#' @param sites Data frame with `environment`, `latitude`, `longitude`,
#'   `start_date`, `end_date`.
#' @param pars NASA POWER parameter codes. Defaults to the six core variables;
#'   see [enviromic_variable_catalog()].
#' @param temporal `"daily"` (default) or `"hourly"` where available.
#' @param max_tries Number of attempts per site before giving up (the POWER API
#'   sometimes times out). Default 3.
#' @param pause Seconds to wait between retries. Default 3.
#' @param cache_dir Optional directory to cache each site's raw POWER response
#'   (keyed by coordinates/dates/parameters). Re-runs read the cache instead of
#'   re-downloading -- a large speed-up for planning. Default `NULL` (no cache).
#' @param workers Number of parallel workers for the fetch. `> 1` uses
#'   \pkg{future.apply} (install it) and is much faster for many site-years.
#'   Default 1 (sequential).
#' @return A long data frame with `environment`, the actual `date`, `day`
#'   (0-based day-after-start computed from that date), and one column per
#'   parameter; or `NULL` if nothing was retrieved. A `provenance` attribute
#'   records request dates, returned dates, row counts, and cache use by site.
#' @seealso [envirotype_features()], [build_enviromic_covariates()].
#' @export
fetch_weather_series <- function(sites, pars = NULL, temporal = "daily",
                                 max_tries = 3L, pause = 3,
                                 cache_dir = NULL, workers = 1L) {
  if (!requireNamespace("nasapower", quietly = TRUE)) {
    warning("Package 'nasapower' not installed; skipping weather fetch.")
    return(NULL)
  }
  need <- c("latitude", "longitude", "start_date", "end_date")
  if (!all(need %in% names(sites))) {
    warning("`sites` needs ", paste(need, collapse = ", "), ".")
    return(NULL)
  }
  if (is.null(pars))
    pars <- c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR",
              "ALLSKY_SFC_SW_DWN", "RH2M")
  if (length(pars) > 20L)
    warning("NASA POWER allows at most 20 daily parameters per request; ",
            length(pars), " requested -- some may be dropped.", call. = FALSE)
  if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  one <- function(i) {
    cache_hit <- FALSE
    cfile <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(
      .power_cache_key(sites$longitude[i], sites$latitude[i], sites$start_date[i],
                       sites$end_date[i], temporal, pars), ".rds")) else NULL
    if (!is.null(cfile) && file.exists(cfile)) {
      d <- readRDS(cfile)
      cache_hit <- TRUE
    } else {
      d <- .power_get_retry(
        community = "ag", pars = pars, temporal_api = temporal,
        lonlat = c(sites$longitude[i], sites$latitude[i]),
        dates = c(as.character(sites$start_date[i]), as.character(sites$end_date[i])),
        max_tries = max_tries, pause = pause, label = sites$environment[i])
      if (!is.null(d) && !is.null(cfile)) saveRDS(d, cfile)
    }
    if (is.null(d)) return(NULL)
    d <- as.data.frame(d)
    actual_date <- .power_response_dates(d)
    date_source <- "response"
    if (is.null(actual_date)) {
      actual_date <- as.Date(sites$start_date[i]) + seq_len(nrow(d)) - 1L
      date_source <- "request_sequence"
      warning("POWER response for '", sites$environment[i],
              "' had no recognised date fields; dates reconstructed from the ",
              "request start.", call. = FALSE)
    }
    if (anyNA(actual_date) || anyDuplicated(actual_date))
      stop("POWER response for '", sites$environment[i],
           "' has invalid or duplicate dates.")
    oo <- order(actual_date)
    actual_date <- actual_date[oo]
    d <- d[oo, , drop = FALSE]
    keep <- intersect(pars, names(d))
    dat <- data.frame(
      environment = as.character(sites$environment[i]),
      date = actual_date,
      day = as.integer(actual_date - as.Date(sites$start_date[i])),
      d[, keep, drop = FALSE],
      check.names = FALSE, stringsAsFactors = FALSE
    )
    expected <- seq(as.Date(sites$start_date[i]), as.Date(sites$end_date[i]),
                    by = "day")
    provenance <- data.frame(
      environment = as.character(sites$environment[i]),
      source = "NASA POWER",
      longitude = as.numeric(sites$longitude[i]),
      latitude = as.numeric(sites$latitude[i]),
      requested_start = as.Date(sites$start_date[i]),
      requested_end = as.Date(sites$end_date[i]),
      returned_start = min(actual_date),
      returned_end = max(actual_date),
      expected_days = length(expected),
      returned_days = length(unique(actual_date)),
      coverage = length(intersect(actual_date, expected)) / length(expected),
      cache_hit = cache_hit,
      date_source = date_source,
      temporal = temporal,
      parameters = paste(pars, collapse = ";"),
      api_endpoint = "NASA POWER daily point API via nasapower",
      nasapower_version = as.character(utils::packageVersion("nasapower")),
      accessed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      cache_file = if (is.null(cfile)) NA_character_ else
        normalizePath(cfile, winslash = "/", mustWork = FALSE),
      cache_md5 = if (!is.null(cfile) && file.exists(cfile))
        unname(tools::md5sum(cfile)) else NA_character_,
      stringsAsFactors = FALSE
    )
    list(data = dat, provenance = provenance)
  }

  idx <- seq_len(nrow(sites))
  out <- if (workers > 1L && requireNamespace("future.apply", quietly = TRUE) &&
             requireNamespace("future", quietly = TRUE)) {
    oplan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(oplan), add = TRUE)
    future.apply::future_lapply(idx, one, future.seed = TRUE)
  } else {
    if (workers > 1L)
      warning("Install 'future' + 'future.apply' for parallel fetch; ",
              "running sequentially.", call. = FALSE)
    lapply(idx, one)
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(NULL)
  ans <- do.call(rbind, lapply(out, `[[`, "data"))
  attr(ans, "provenance") <- do.call(rbind, lapply(out, `[[`, "provenance"))
  ans
}

# Deterministic, filesystem-safe cache key for one POWER request.
.power_cache_key <- function(lon, lat, start, end, temporal, pars) {
  version <- if (requireNamespace("nasapower", quietly = TRUE))
    as.character(utils::packageVersion("nasapower")) else "unknown"
  substr(gsub("[^0-9A-Za-z]", "",
              paste0("power", round(lon, 6), round(lat, 6), start, end, temporal,
                     paste(pars, collapse = ""), version)), 1, 180)
}


# Call nasapower::get_power with retries and a longer timeout (the POWER service
# intermittently times out on the default). Returns the data frame or NULL.
.power_get_retry <- function(..., max_tries = 3L, pause = 3, label = "") {
  old <- getOption("timeout"); on.exit(options(timeout = old), add = TRUE)
  options(timeout = max(120, old))
  d <- NULL
  for (k in seq_len(max_tries)) {
    d <- tryCatch(nasapower::get_power(...), error = function(e) e)
    if (!inherits(d, "error")) return(d)
    if (k < max_tries) Sys.sleep(pause)
  }
  warning(sprintf("NASA POWER fetch failed for '%s' after %d tries: %s",
                  label, max_tries, conditionMessage(d)), call. = FALSE)
  NULL
}

# Extract actual dates from known nasapower response schemas.
.power_response_dates <- function(d) {
  nm <- names(d)
  date_name <- intersect(c("DATE", "Date", "date", "YYYYMMDD"), nm)
  if (length(date_name)) {
    x <- as.character(d[[date_name[1L]]])
    out <- suppressWarnings(as.Date(x))
    if (all(is.na(out)) && all(grepl("^[0-9]{8}$", x)))
      out <- as.Date(x, format = "%Y%m%d")
    if (!all(is.na(out))) return(out)
  }
  yy <- intersect(c("YEAR", "year", "Year"), nm)
  mm <- intersect(c("MO", "MONTH", "month", "Month"), nm)
  dd <- intersect(c("DY", "DAY", "day_of_month", "Day"), nm)
  if (length(yy) && length(mm) && length(dd)) {
    return(as.Date(sprintf(
      "%04d-%02d-%02d",
      as.integer(d[[yy[1L]]]), as.integer(d[[mm[1L]]]),
      as.integer(d[[dd[1L]]])
    )))
  }
  rn <- rownames(d)
  if (!is.null(rn) && all(grepl("^[0-9]{8}$", rn)))
    return(as.Date(rn, format = "%Y%m%d"))
  NULL
}


# ---- helpers ----------------------------------------------------------------

.make_windows <- function(windows, day) {
  lo0 <- min(day, na.rm = TRUE); hi0 <- max(day, na.rm = TRUE)
  if (is.list(windows) && !is.null(names(windows))) {
    n <- length(windows)
    return(lapply(seq_len(n), function(i) list(
      name = names(windows)[i], lo = windows[[i]][1], hi = windows[[i]][2],
      last = i == n)))
  }
  if (length(windows) == 1L) {
    K <- max(1L, as.integer(windows))
    br <- seq(lo0, hi0, length.out = K + 1L)
  } else {
    br <- unique(sort(c(lo0, as.numeric(windows), hi0)))
  }
  n <- length(br) - 1L
  lapply(seq_len(n), function(i) list(
    name = paste0("w", i), lo = br[i], hi = br[i + 1L], last = i == n))
}

.environment_feature_windows <- function(phenology, environment, windows, day,
                                         date = NULL) {
  if (is.null(phenology)) return(.make_windows(windows, day))
  p <- as.data.frame(phenology, stringsAsFactors = FALSE)
  need <- c("environment", "stage")
  if (!all(need %in% names(p)))
    stop("`phenology` needs environment and stage columns.")
  p <- p[as.character(p$environment) == environment, , drop = FALSE]
  if (!nrow(p)) return(.make_windows(windows, day))
  if (all(c("start_day", "end_day") %in% names(p))) {
    lo <- as.numeric(p$start_day); hi <- as.numeric(p$end_day)
  } else if (all(c("start_date", "end_date") %in% names(p))) {
    if (is.null(date))
      stop("Date-based `phenology` requires a date column in `daily`.")
    origin <- min(as.Date(date))
    lo <- as.integer(as.Date(p$start_date) - origin)
    hi <- as.integer(as.Date(p$end_date) - origin)
  } else {
    stop("`phenology` needs start_day/end_day or start_date/end_date.")
  }
  if (any(!is.finite(lo)) || any(!is.finite(hi)) || any(hi < lo))
    stop("Invalid phenology limits for environment '", environment, "'.")
  oo <- order(lo, hi)
  lapply(seq_along(oo), function(k) {
    i <- oo[k]
    list(name = as.character(p$stage[i]), lo = lo[i], hi = hi[i],
         last = TRUE)
  })
}

.agg_stat <- function(x, st) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  switch(st,
         mean = mean(x), sum = sum(x), min = min(x), max = max(x),
         sd = stats::sd(x), median = stats::median(x),
         cv = if (mean(x) == 0) NA_real_ else stats::sd(x) / mean(x),
         mean(x))
}

.count_threshold <- function(x, sp) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  if (!is.null(sp$above)) return(sum(x > sp$above))
  if (!is.null(sp$below)) return(sum(x < sp$below))
  if (!is.null(sp$between)) return(sum(x >= sp$between[1] & x <= sp$between[2]))
  NA_real_
}

.stress_severity <- function(x, sp) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  if (!is.null(sp$above)) return(sum(pmax(0, x - sp$above)))
  if (!is.null(sp$below)) return(sum(pmax(0, sp$below - x)))
  stop("Stress definition needs `above` or `below`.")
}

.spell_summary <- function(x, sp) {
  event <- rep(FALSE, length(x))
  ok <- is.finite(x)
  if (!is.null(sp$above)) event[ok] <- x[ok] > sp$above
  else if (!is.null(sp$below)) event[ok] <- x[ok] < sp$below
  else stop("Spell definition needs `above` or `below`.")
  event[!ok] <- FALSE
  rr <- rle(event)
  lengths <- rr$lengths[rr$values]
  c(longest = if (length(lengths)) max(lengths) else 0,
    runs = length(lengths))
}

.gdd_sum <- function(x, sp) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  up <- if (is.null(sp$upper)) Inf else sp$upper
  sum(pmax(0, pmin(x, up) - sp$base))
}

# B-spline coefficients over normalised season time (0-1), so different season
# lengths are compared on a common phenological basis.
.functional_coef <- function(day, y, df_basis = 4L) {
  ok <- is.finite(day) & is.finite(y); day <- day[ok]; y <- y[ok]
  if (length(unique(day)) < df_basis + 1L) return(rep(NA_real_, df_basis))
  u <- (day - min(day)) / (max(day) - min(day) + 1e-9)
  B <- splines::bs(u, df = df_basis, intercept = TRUE)
  tryCatch(as.numeric(solve(crossprod(B), crossprod(B, y))),
           error = function(e) rep(NA_real_, ncol(B)))
}


#' Bias-correct gridded weather with on-station observations
#'
#' @description
#' Calibrates downloaded daily weather against overlapping station records,
#' separately by environment and variable. Additive mean-bias correction is
#' used for temperature-like variables; a non-negative ratio correction is used
#' for precipitation-like variables. Alternatively, a linear calibration is
#' fitted. The downloaded series remains complete and the station observations
#' replace overlapping values.
#'
#' @param downloaded,station Long daily data frames.
#' @param method `"mean_bias"` or `"linear"`.
#' @param variables Optional variables to calibrate; defaults to shared numeric
#'   weather columns.
#' @param min_overlap Minimum paired observations required per fit.
#' @param env_col,date_col,day_col Key-column names. Date is preferred.
#'
#' @return A list with corrected `data` and a fit `diagnostics` table.
#' @export
calibrate_weather_series <- function(downloaded, station,
                                     method = c("mean_bias", "linear"),
                                     variables = NULL, min_overlap = 20L,
                                     env_col = "environment",
                                     date_col = "date", day_col = "day") {
  method <- match.arg(method)
  dl <- as.data.frame(downloaded, check.names = FALSE,
                      stringsAsFactors = FALSE)
  st <- as.data.frame(station, check.names = FALSE,
                      stringsAsFactors = FALSE)
  if (!env_col %in% names(dl) || !env_col %in% names(st))
    stop("Both data sets need an environment column.")
  key_time <- if (date_col %in% names(dl) && date_col %in% names(st))
    date_col else day_col
  if (!key_time %in% names(dl) || !key_time %in% names(st))
    stop("Both data sets need a shared date or day key.")
  if (key_time == date_col) {
    dl[[date_col]] <- as.Date(dl[[date_col]])
    st[[date_col]] <- as.Date(st[[date_col]])
    if (anyNA(dl[[date_col]]) || anyNA(st[[date_col]]))
      stop("Downloaded/station date keys contain invalid dates.")
  }
  keys <- c(env_col, key_time)
  if (anyDuplicated(dl[keys]) || anyDuplicated(st[keys]))
    stop("Downloaded and station data must have unique environment-time keys.")
  shared <- intersect(names(dl), names(st))
  shared <- setdiff(shared, keys)
  shared <- shared[vapply(shared, function(nm)
    is.numeric(dl[[nm]]) && is.numeric(st[[nm]]), logical(1))]
  if (is.null(variables)) variables <- shared
  absent <- setdiff(variables, shared)
  if (length(absent))
    stop("Calibration variables are not shared numeric columns: ",
         paste(absent, collapse = ", "), ".")
  if (length(min_overlap) != 1L || !is.finite(min_overlap) ||
      min_overlap < 3)
    stop("`min_overlap` must be at least 3.")

  corrected <- dl
  diagnostics <- list()
  for (e in intersect(unique(as.character(dl[[env_col]])),
                      unique(as.character(st[[env_col]])))) {
    id <- which(as.character(dl[[env_col]]) == e)
    is <- which(as.character(st[[env_col]]) == e)
    mt <- match(dl[[key_time]][id], st[[key_time]][is])
    for (v in variables) {
      x <- dl[[v]][id]
      y <- st[[v]][is[mt]]
      ok <- is.finite(x) & is.finite(y)
      status <- "insufficient_overlap"
      intercept <- 0; slope <- 1
      if (sum(ok) >= min_overlap) {
        precip_like <- grepl("PREC|RAIN", v, ignore.case = TRUE)
        if (method == "linear") {
          fit <- stats::lm(y[ok] ~ x[ok])
          cf <- stats::coef(fit)
          if (all(is.finite(cf))) {
            intercept <- unname(cf[1L]); slope <- unname(cf[2L])
            status <- "linear"
          }
        } else if (precip_like) {
          den <- mean(x[ok])
          slope <- if (is.finite(den) && den > 0)
            mean(y[ok]) / den else 1
          intercept <- 0
          status <- "multiplicative_mean_bias"
        } else {
          intercept <- mean(y[ok] - x[ok])
          slope <- 1
          status <- "additive_mean_bias"
        }
      }
      pred <- intercept + slope * x
      if (grepl("PREC|RAIN", v, ignore.case = TRUE)) pred <- pmax(0, pred)
      pred[!is.na(mt) & is.finite(y)] <- y[!is.na(mt) & is.finite(y)]
      corrected[[v]][id] <- pred
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        environment = e, variable = v, method = status,
        n_overlap = sum(ok), intercept = intercept, slope = slope,
        rmse_before = if (sum(ok)) sqrt(mean((y[ok] - x[ok])^2)) else NA_real_,
        rmse_after = if (sum(ok))
          sqrt(mean((y[ok] - (intercept + slope * x[ok]))^2)) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    data = corrected,
    diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else
      data.frame()
  )
}
