#' Characterise candidate sites from past seasons (planning-stage envirotyping)
#'
#' @description
#' At the design stage the upcoming season has not happened, so its weather
#' cannot be fetched. `historical_envirotype()` characterises each candidate site
#' from **several past seasons** over the intended calendar window: it builds
#' [envirotype_features()] for every year and summarises them into (i) the
#' **typical** environmental profile (the across-year mean, used to predict the
#' environment relationship and select sites) and (ii) the **interannual
#' variability** of each feature (the across-year standard deviation or CV -- the
#' site's environmental *risk*, which feeds robust design). This is the
#' planning-appropriate substitute for this-season data.
#'
#' @details
#' For the network path, supply `sites` (with `latitude`/`longitude`) and a
#' calendar `window` of month-day strings; the growing window is applied to each
#' year in `years`, the daily series fetched with [fetch_weather_series()], and
#' features computed with [envirotype_features()]. For offline use (or your own
#' records) pass `daily_by_year`, a named list of long daily data frames (one per
#' year). The typical profile plugs into [build_environment_relationship()]; the
#' variability profile flags sites whose environment swings year to year, which
#' informs descriptive mega-environment stability and which separate weather
#' structures should enter [robust_scenarios()]. It does not identify genetic
#' `Sigma_E` without historical MET values.
#'
#' @param sites Data frame with `environment`, `latitude`, `longitude` (network
#'   path). May be `NULL` when `daily_by_year` is given.
#' @param years Integer vector of past years to characterise the site with.
#' @param window The growing window, as **month-day** strings. Either a length-2
#'   vector `c(start_md, end_md)` applied to every site, OR a data frame with
#'   columns `environment`, `start_md`, `end_md` for **site-specific** windows
#'   (e.g. Peru's season differs from Colombia's), OR supply `start_md`/`end_md`
#'   columns directly on `sites`. A window whose `end_md` precedes its `start_md`
#'   is treated as crossing into the next calendar year (e.g. `"11-01"`->`"02-28"`).
#' @param station_daily Optional on-station (manually collected) daily weather to
#'   merge with the downloaded series: a long data frame with `environment`,
#'   `day` (0-based day-after-start), optionally `year`, and weather columns.
#'   Station values take precedence; downloaded data fills the variables and days
#'   the station did not measure.
#' @param pars NASA POWER parameter codes (see [enviromic_variable_catalog()]).
#' @param envirotype A list of arguments passed to [envirotype_features()]
#'   (`windows`, `stats`, `thresholds`, `gdd`, `functional`).
#' @param variability `"sd"` (default) or `"cv"` for the interannual-variability
#'   summary.
#' @param daily_by_year Optional named list (by year) of long daily data frames
#'   to use instead of fetching -- enables offline use and reproducibility.
#' @param max_tries Retries per site-year for the NASA POWER fetch (the service
#'   intermittently times out). Default 3.
#' @param cache_dir,workers Passed to [fetch_weather_series()] to cache responses
#'   (fast re-runs) and fetch in parallel.
#' @param station_correction Optional correction of downloaded weather against
#'   overlapping station records: `"none"`, `"mean_bias"`, or `"linear"`.
#' @param summary_components Historical summaries included in `combined`.
#'   Available components are `"mean"`, `"sd"` (or `"cv"` according to
#'   `variability`), `"q10"`, `"q50"`, `"q90"`, `"trend"`, and
#'   `"n_observed"`.
#' @param min_years Minimum observed seasons required per environment-feature.
#'   Summaries with fewer seasons are set to missing and recorded in
#'   `year_coverage`.
#' @param feature_function `"generic"` uses [envirotype_features()];
#'   `"rice"` uses [rice_weather_features()] and therefore adds rice-specific
#'   stress physiology before across-year summarisation.
#' @param min_daily_coverage Minimum within-season coverage required for every
#'   weather variable at every environment.
#' @param daily_missing_action `"impute"` (linear within site, then median),
#'   `"warn"`, or `"error"`.
#' @param weather_ranges Named physical plausibility ranges for daily weather.
#'   Defaults cover common NASA POWER variables.
#' @return A list with `typical` (environment-by-feature matrix, across-year
#'   mean), `variability` (across-year SD or CV), `combined` (typical +
#'   variability columns, suffixed `_iav`, so the between-year variation can
#'   enter the relationship matrix), `n_years`, and `per_year` (the list of
#'   per-year feature matrices). Build `D` from `typical` (mean conditions),
#'   `combined` (mean + interannual risk), or the consensus of the per-year
#'   matrices from [assess_envirotype_stability()].
#' @seealso [envirotype_features()], [fetch_weather_series()],
#'   [build_environment_relationship()], [robust_scenarios()].
#' @examples
#' # Offline: three past seasons of daily data supplied directly.
#' mk <- function(base) rbind(
#'   data.frame(environment = "E1", day = 0:59,
#'              T2M = base + rnorm(60), T2M_MAX = base + 6 + rnorm(60)),
#'   data.frame(environment = "E2", day = 0:59,
#'              T2M = base + 4 + rnorm(60), T2M_MAX = base + 10 + rnorm(60)))
#' hist <- historical_envirotype(
#'   years = 2020:2022,
#'   daily_by_year = list(`2020` = mk(24), `2021` = mk(26), `2022` = mk(22)),
#'   envirotype = list(windows = 2, stats = "mean"))
#' hist$typical
#' hist$variability
#' @export
historical_envirotype <- function(sites = NULL, years = NULL, window = NULL,
                                  pars = NULL, envirotype = list(),
                                  variability = c("sd", "cv"),
                                  daily_by_year = NULL, max_tries = 3L,
                                  cache_dir = NULL, workers = 1L,
                                  station_daily = NULL,
                                  station_correction = c("none", "mean_bias",
                                                         "linear"),
                                  summary_components = c("mean", "sd", "q10",
                                                         "q50", "q90", "trend"),
                                  min_years = 3L,
                                  feature_function = c("generic", "rice"),
                                  min_daily_coverage = 0.80,
                                  daily_missing_action = c("impute", "warn",
                                                           "error"),
                                  weather_ranges = .default_weather_ranges()) {
  variability <- match.arg(variability)
  station_correction <- match.arg(station_correction)
  allowed_components <- c("mean", "sd", "cv", "q10", "q50", "q90",
                          "trend", "n_observed")
  if (any(!summary_components %in% allowed_components))
    stop("Unknown `summary_components`: ",
         paste(setdiff(summary_components, allowed_components), collapse = ", "),
         ".")
  if (!is.numeric(min_years) || length(min_years) != 1L ||
      !is.finite(min_years) || min_years < 1)
    stop("`min_years` must be a positive integer.")
  min_years <- as.integer(min_years)
  feature_function <- match.arg(feature_function)
  daily_missing_action <- match.arg(daily_missing_action)

  if (!is.null(daily_by_year)) {
    dailies <- daily_by_year
  } else {
    if (is.null(sites) || is.null(years))
      stop("Supply `sites` + `years` (+ `window`), or `daily_by_year`.")
    wr <- .resolve_windows(sites, window)              # per-site growing window
    dailies <- lapply(years, function(y) {
      s <- sites
      s$start_date <- paste0(y, "-", wr$start_md)
      end_year     <- ifelse(wr$end_md < wr$start_md, y + 1L, y)  # cross-year window
      s$end_date   <- paste0(end_year, "-", wr$end_md)
      dl <- fetch_weather_series(s, pars = pars, max_tries = max_tries,
                                 cache_dir = cache_dir, workers = workers)
      power_provenance <- attr(dl, "provenance")
      ## merge in on-station (manually collected) daily data for this year, if any
      st <- if (!is.null(station_daily) && "year" %in% names(station_daily))
        station_daily[station_daily$year == y,
                      setdiff(names(station_daily), "year"), drop = FALSE] else
        station_daily
      correction <- NULL
      if (!is.null(dl) && !is.null(st) && station_correction != "none") {
        correction <- calibrate_weather_series(
          dl, st, method = station_correction
        )
        dl <- correction$data
      }
      merged <- .merge_daily_series(dl, st)
      attr(merged, "station_calibration") <- if (is.null(correction))
        NULL else correction$diagnostics
      attr(merged, "provenance") <- power_provenance
      merged
    })
    names(dailies) <- as.character(years)
  }

  station_diagnostics <- lapply(dailies, attr, "station_calibration")
  weather_provenance <- lapply(dailies, attr, "provenance")
  daily_qc <- lapply(seq_along(dailies), function(k) {
    dy <- dailies[[k]]
    if (is.null(dy)) return(NULL)
    .qc_daily_weather(
      dy, source = paste0("weather_", names(dailies)[k]),
      min_coverage = min_daily_coverage,
      missing_action = daily_missing_action,
      ranges = weather_ranges
    )
  })
  names(daily_qc) <- names(dailies)
  dailies <- lapply(daily_qc, function(z) if (is.null(z)) NULL else z$data)

  per_year <- lapply(dailies, function(dy) {
    if (is.null(dy)) return(NULL)
    fun <- if (feature_function == "rice") rice_weather_features else
      envirotype_features
    do.call(fun, c(list(daily = dy), envirotype))
  })
  per_year <- per_year[!vapply(per_year, is.null, logical(1))]
  if (!length(per_year)) {
    warning("No usable seasons; nothing to summarise.", call. = FALSE)
    return(NULL)
  }

  agg <- .aggregate_across_years(per_year, variability, min_years = min_years)
  components <- list(
    mean = agg$typical,
    sd = if (variability == "sd") agg$variability else agg$sd,
    cv = if (variability == "cv") agg$variability else agg$cv,
    q10 = agg$q10, q50 = agg$q50, q90 = agg$q90,
    trend = agg$trend, n_observed = agg$n_observed
  )
  selected <- components[summary_components]
  suffix <- c(mean = "", sd = "_iav", cv = "_iav", q10 = "_q10",
              q50 = "_q50", q90 = "_q90", trend = "_trend",
              n_observed = "_n_years")
  selected <- Map(function(M, s) {
    M <- as.matrix(M)
    colnames(M) <- paste0(colnames(M), s)
    M
  }, selected, suffix[names(selected)])
  combined <- do.call(cbind, selected)
  c(agg, list(combined = combined, n_years = length(per_year),
              per_year = per_year,
              summary_components = summary_components,
              station_calibration = station_diagnostics,
              weather_provenance = weather_provenance,
              daily_audit = .bind_qc_component(daily_qc, "audit"),
              daily_imputation = .bind_qc_component(daily_qc, "imputation"),
              daily_qc_provenance = .bind_qc_component(daily_qc,
                                                       "provenance")))
}


#' Assess whether the environment relationship is stable across years
#'
#' @description
#' Is the environmental kinship a reliable planning input, or does it swing from
#' year to year? `assess_envirotype_stability()` builds the environment
#' relationship matrix `D` **separately for each past year** (from
#' `historical_envirotype()$per_year`) and measures how consistent those yearly
#' matrices are -- a Mantel-type correlation between the off-diagonal similarities
#' of each pair of years. High mean stability means the site relationships (and
#' the mega-environment structure) are dependable; low stability warns that the
#' descriptive environmental structure is unstable and that distinct yearly
#' weather kernels should be retained as robust-design stress tests. It does not
#' by itself estimate genetic `Sigma_E`.
#'
#' @param hist A [historical_envirotype()] result (uses its `per_year`).
#' @param kernel,variables,weights Passed to [build_environment_relationship()]
#'   when building each year's `D`.
#' @param consensus_method How to combine the per-year matrices into
#'   `consensus_D`; passed to [consensus_relationship()]. Default
#'   `"rv_weighted"` (STATIS -- down-weights anomalous years), rather than a plain
#'   mean.
#' @param bandwidth_multipliers Gaussian bandwidth ensemble around the median.
#' @param n_boot Number of environment-pair bootstrap replicates for the
#'   confidence interval of mean stability.
#' @param conf_level Bootstrap confidence level.
#' @param seed Optional bootstrap seed.
#' @return A list with `per_year_D` (a `D` per year), `stability` (year-by-year
#'   correlation matrix of the off-diagonal similarities), `mean_stability` (the
#'   mean across year pairs), `consensus_D` (a robust consensus relationship to
#'   design on), and `consensus_weights` (per-year weights when applicable).
#' @seealso [historical_envirotype()], [consensus_relationship()],
#'   [build_environment_relationship()], [robust_scenarios()].
#' @export
assess_envirotype_stability <- function(hist, kernel = "gaussian",
                                        variables = NULL, weights = NULL,
                                        consensus_method = "rv_weighted",
                                        bandwidth_multipliers = c(0.5, 1, 2),
                                        n_boot = 500L, conf_level = 0.95,
                                        seed = NULL) {
  per <- hist$per_year
  if (is.null(per) || length(per) < 2L)
    stop("Need at least two years in `hist$per_year` to assess stability.")
  common <- Reduce(intersect, lapply(per, rownames))
  if (length(common) < 3L)
    stop("Fewer than 3 environments are shared across all years.")
  Ds <- lapply(per, function(W) {
    X <- W[common, , drop = FALSE]
    if (kernel == "gaussian" && length(bandwidth_multipliers) > 1L) {
      base <- build_environment_relationship(
        X, source = "enviromic", kernel = kernel, variables = variables,
        weights = weights
      )
      h0 <- attr(base, "bandwidth")
      set <- lapply(h0 * bandwidth_multipliers, function(h)
        build_environment_relationship(
          X, source = "enviromic", kernel = kernel, variables = variables,
          weights = weights, bandwidth = h
        ))
      Reduce(`+`, set) / length(set)
    } else {
      build_environment_relationship(
        X, source = "enviromic", kernel = kernel, variables = variables,
        weights = weights
      )
    }
  })
  ny <- length(Ds); yrs <- names(per)
  offd <- function(D) D[upper.tri(D)]
  cormat <- matrix(1, ny, ny, dimnames = list(yrs, yrs))
  for (i in seq_len(ny)) for (j in seq_len(ny)) if (i < j) {
    cc <- suppressWarnings(stats::cor(offd(Ds[[i]]), offd(Ds[[j]])))
    cormat[i, j] <- cormat[j, i] <- if (is.finite(cc)) cc else NA_real_
  }
  # Robust consensus relationship (STATIS by default: down-weights anomalous
  # years), rather than a plain arithmetic mean. A more dependable design input
  # than a D built from year-averaged features -- so the analysis is actionable.
  names(Ds) <- names(per)
  consensus <- consensus_relationship(
    Ds, method = consensus_method, use_off_diagonal = TRUE
  )
  mean_stability <- mean(cormat[upper.tri(cormat)], na.rm = TRUE)
  ci <- c(lower = NA_real_, upper = NA_real_)
  boot <- numeric()
  if (n_boot > 0L) {
    if (!is.numeric(n_boot) || length(n_boot) != 1L ||
        !is.finite(n_boot) || n_boot < 1)
      stop("`n_boot` must be a non-negative integer.")
    if (!is.numeric(conf_level) || length(conf_level) != 1L ||
        !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
      stop("`conf_level` must lie strictly between 0 and 1.")
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }, add = TRUE)
      set.seed(seed)
    }
    pair_vectors <- lapply(Ds, offd)
    np <- length(pair_vectors[[1L]])
    boot <- replicate(as.integer(n_boot), {
      ii <- sample.int(np, np, replace = TRUE)
      vals <- numeric()
      for (a in seq_len(ny - 1L)) for (b in (a + 1L):ny) {
        cc <- suppressWarnings(stats::cor(
          pair_vectors[[a]][ii], pair_vectors[[b]][ii]
        ))
        if (is.finite(cc)) vals <- c(vals, cc)
      }
      if (length(vals)) mean(vals) else NA_real_
    })
    alpha <- (1 - conf_level) / 2
    ci <- stats::quantile(boot, c(alpha, 1 - alpha), na.rm = TRUE,
                          names = FALSE)
    names(ci) <- c("lower", "upper")
  }
  list(per_year_D = Ds, stability = cormat,
       mean_stability = mean_stability,
       mean_stability_ci = ci, bootstrap_stability = boot,
       consensus_D = consensus,
       consensus_weights = attr(consensus, "weights"))
}


# Resolve a per-site growing window. `window` may be a length-2 c(start_md,
# end_md) applied to all sites, a data frame with columns environment/start_md/
# end_md, or start_md/end_md columns already on `sites`.
.resolve_windows <- function(sites, window) {
  envs <- as.character(sites$environment)
  if (is.data.frame(window)) {
    if (!all(c("environment", "start_md", "end_md") %in% names(window)))
      stop("`window` data frame needs columns environment, start_md, end_md.")
    idx <- match(envs, as.character(window$environment))
    if (anyNA(idx)) stop("`window` is missing rows for some environments.")
    return(data.frame(environment = envs,
                      start_md = as.character(window$start_md)[idx],
                      end_md = as.character(window$end_md)[idx],
                      stringsAsFactors = FALSE))
  }
  if (all(c("start_md", "end_md") %in% names(sites)))
    return(data.frame(environment = envs,
                      start_md = as.character(sites$start_md),
                      end_md = as.character(sites$end_md),
                      stringsAsFactors = FALSE))
  if (is.null(window) || length(window) != 2L)
    stop("Supply `window` = c(start_md, end_md), a per-site data frame ",
         "(environment/start_md/end_md), or start_md/end_md columns in `sites`.")
  data.frame(environment = envs, start_md = window[1], end_md = window[2],
             stringsAsFactors = FALSE)
}

# Merge downloaded and on-station daily series by environment + date (preferred)
# or day. On-station
# values take precedence for shared variables; downloaded fills variables (and
# rows) the station lacks -- so measured data and gap-filling online data combine.
.merge_daily_series <- function(dl, st) {
  if (is.null(st)) return(dl)
  if (is.null(dl)) return(st)
  st <- as.data.frame(st, check.names = FALSE, stringsAsFactors = FALSE)
  time_key <- if ("date" %in% names(dl) && "date" %in% names(st))
    "date" else "day"
  key <- c("environment", time_key)
  if (!all(key %in% names(st))) return(dl)
  if (time_key == "date") {
    dl$date <- as.Date(dl$date)
    st$date <- as.Date(st$date)
    if (anyNA(dl$date) || anyNA(st$date))
      stop("Downloaded/station date keys contain invalid dates.")
  }
  m <- merge(dl, st, by = key, all = TRUE, suffixes = c(".dl", ".st"))
  shared <- intersect(setdiff(names(dl), key), setdiff(names(st), key))
  for (v in shared) {
    vs <- paste0(v, ".st"); vd <- paste0(v, ".dl")
    m[[v]] <- ifelse(!is.na(m[[vs]]), m[[vs]], m[[vd]])   # station priority
    m[[vs]] <- NULL; m[[vd]] <- NULL
  }
  m[order(m$environment, m[[time_key]]), , drop = FALSE]
}


# Across-year mean (typical) and SD/CV (interannual variability) of aligned
# per-year env-by-feature matrices.
.aggregate_across_years <- function(mats, variability = "sd", min_years = 1L) {
  envs  <- Reduce(union, lapply(mats, rownames))
  feats <- Reduce(union, lapply(mats, colnames))
  arr <- array(NA_real_, dim = c(length(envs), length(feats), length(mats)),
               dimnames = list(envs, feats, NULL))
  for (k in seq_along(mats)) {
    m <- mats[[k]]
    arr[rownames(m), colnames(m), k] <- m
  }
  typ <- apply(arr, c(1, 2), function(z) {
    z <- z[is.finite(z)]; if (!length(z)) NA_real_ else mean(z)
  })
  sdv <- apply(arr, c(1, 2), function(z) {
    z <- z[is.finite(z)]; if (length(z) < 2L) NA_real_ else stats::sd(z)
  })
  var_mat <- if (variability == "cv") {
    cv <- sdv / abs(typ); cv[!is.finite(cv)] <- NA_real_; cv
  } else sdv
  cv <- sdv / abs(typ); cv[!is.finite(cv)] <- NA_real_
  qfun <- function(prob) apply(arr, c(1, 2), function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) NA_real_ else
      as.numeric(stats::quantile(z, prob, names = FALSE, type = 8))
  })
  nobs <- apply(arr, c(1, 2), function(z) sum(is.finite(z)))
  yr <- suppressWarnings(as.numeric(names(mats)))
  if (length(yr) != length(mats) || any(!is.finite(yr)))
    yr <- seq_along(mats)
  trend <- apply(arr, c(1, 2), function(z) {
    ok <- is.finite(z)
    if (sum(ok) < max(2L, min_years)) return(NA_real_)
    unname(stats::coef(stats::lm(z[ok] ~ yr[ok]))[2L])
  })
  low <- nobs < min_years
  mats_out <- list(
    typical = typ, variability = var_mat, sd = sdv, cv = cv,
    q10 = qfun(0.10), q50 = qfun(0.50), q90 = qfun(0.90),
    trend = trend, n_observed = nobs
  )
  for (nm in setdiff(names(mats_out), "n_observed"))
    mats_out[[nm]][low] <- NA_real_
  for (nm in names(mats_out))
    dimnames(mats_out[[nm]]) <- list(envs, feats)
  mats_out$year_coverage <- data.frame(
    environment = rep(envs, each = length(feats)),
    feature = rep(feats, times = length(envs)),
    n_observed = as.vector(t(nobs)),
    sufficient = as.vector(t(nobs >= min_years)),
    stringsAsFactors = FALSE
  )
  mats_out
}

.default_weather_ranges <- function() {
  list(
    T2M = c(-30, 60), T2M_MAX = c(-25, 65), T2M_MIN = c(-40, 55),
    T2MDEW = c(-50, 45), RH2M = c(0, 100),
    PRECTOTCORR = c(0, 1000), ALLSKY_SFC_SW_DWN = c(0, 50),
    WS2M = c(0, 100), GWETTOP = c(0, 1), GWETROOT = c(0, 1),
    GWETPROF = c(0, 1), EVPTRNS = c(-20, 100)
  )
}

.qc_daily_weather <- function(daily, source, min_coverage,
                              missing_action, ranges) {
  source_provenance <- attr(daily, "provenance")
  d <- as.data.frame(daily, check.names = FALSE, stringsAsFactors = FALSE)
  synthetic_date <- !"date" %in% names(d)
  if (synthetic_date) {
    if (!"day" %in% names(d))
      stop("Daily weather requires a date or day column before QC.")
    d$date <- as.Date("1970-01-01") + as.integer(d$day)
  } else {
    d$date <- as.Date(d$date)
  }
  expected <- do.call(rbind, lapply(unique(as.character(d$environment)),
                                    function(e) {
    z <- d$date[as.character(d$environment) == e]
    if (!is.null(source_provenance) &&
        all(c("environment", "requested_start", "requested_end") %in%
            names(source_provenance)) &&
        e %in% source_provenance$environment) {
      pp <- source_provenance[
        match(e, as.character(source_provenance$environment)), ,
        drop = FALSE]
      dates <- seq(as.Date(pp$requested_start), as.Date(pp$requested_end),
                   by = "day")
    } else {
      dates <- seq(min(z), max(z), by = "day")
    }
    data.frame(environment = e, date = dates,
               stringsAsFactors = FALSE)
  }))
  d <- merge(expected, d, by = c("environment", "date"), all.x = TRUE,
             sort = FALSE)
  if ("day" %in% names(d)) {
    for (e in unique(as.character(d$environment))) {
      ii <- which(as.character(d$environment) == e)
      origin <- min(d$date[ii])
      d$day[ii] <- as.integer(d$date[ii] - origin)
    }
  }
  q <- qc_environmental_data(
    d, environments = unique(as.character(d$environment)),
    date_col = "date", expected_dates = expected,
    ranges = ranges[intersect(names(ranges), names(d))],
    min_coverage = min_coverage,
    duplicate_action = "error", missing_action = missing_action,
    impute = if (missing_action == "impute") "linear" else "none",
    add_missing_indicators = FALSE, source = source
  )
  if (synthetic_date) q$data$date <- NULL
  q
}
