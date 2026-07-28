#' Quality-control environmental source data
#'
#' @description
#' Applies explicit, auditable quality control before weather, soil, management,
#' or other environmental data are converted to kernels. The function validates
#' environment coverage and duplicate keys, converts out-of-range observations
#' to missing, applies a declared imputation rule, and records every affected
#' environment-variable combination. It never silently replaces missing values.
#'
#' @param data A data frame containing environmental observations.
#' @param environments Optional complete character vector of required
#'   environments.
#' @param environment_col Name of the environment column.
#' @param date_col Optional date column. When supplied, duplicate and coverage
#'   checks use environment-date keys and the output is ordered by them.
#' @param expected_dates Optional vector of expected dates, or a data frame with
#'   `environment` and the column named by `date_col`. If omitted, coverage is
#'   relative to the observed union of dates.
#' @param ranges Optional named list. Each element is a length-two valid range
#'   for a numeric variable; values outside it are set to `NA`.
#' @param min_coverage Minimum acceptable observed fraction for each
#'   environment-variable combination.
#' @param duplicate_action `"error"` (default) or `"aggregate"`. Aggregation
#'   averages numeric values and takes the first non-missing non-numeric value.
#' @param missing_action `"impute"`, `"warn"`, or `"error"`.
#' @param impute Imputation used when `missing_action = "impute"`:
#'   `"linear"` interpolates within environment and then uses the cross-
#'   environment median, `"median"` uses the variable median, and `"none"`
#'   leaves values missing.
#' @param add_missing_indicators Add `<variable>__was_missing` columns.
#' @param source Short source label stored in the audit and provenance.
#'
#' @return A list with `data`, `audit`, `provenance`, and `imputation`.
#' @export
qc_environmental_data <- function(data, environments = NULL,
                                  environment_col = "environment",
                                  date_col = NULL, expected_dates = NULL,
                                  ranges = NULL, min_coverage = 0.80,
                                  duplicate_action = c("error", "aggregate"),
                                  missing_action = c("impute", "warn", "error"),
                                  impute = c("linear", "median", "none"),
                                  add_missing_indicators = TRUE,
                                  source = "environmental") {
  duplicate_action <- match.arg(duplicate_action)
  missing_action <- match.arg(missing_action)
  impute <- match.arg(impute)
  d <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  if (!environment_col %in% names(d))
    stop("`data` has no `", environment_col, "` column.")
  if (!is.null(date_col) && !date_col %in% names(d))
    stop("`data` has no `", date_col, "` column.")
  if (!is.numeric(min_coverage) || length(min_coverage) != 1L ||
      !is.finite(min_coverage) || min_coverage < 0 || min_coverage > 1)
    stop("`min_coverage` must be one finite number in [0, 1].")

  d[[environment_col]] <- as.character(d[[environment_col]])
  if (anyNA(d[[environment_col]]) || any(!nzchar(d[[environment_col]])))
    stop("Environment identifiers must be non-missing and non-empty.")
  if (is.null(environments)) environments <- unique(d[[environment_col]])
  environments <- as.character(environments)
  if (anyDuplicated(environments))
    stop("`environments` must contain unique identifiers.")
  unknown <- setdiff(unique(d[[environment_col]]), environments)
  if (length(unknown))
    stop("Source contains environments outside `environments`: ",
         paste(unknown, collapse = ", "), ".")

  if (!is.null(date_col)) {
    d[[date_col]] <- as.Date(d[[date_col]])
    if (anyNA(d[[date_col]]))
      stop("`", date_col, "` contains missing or invalid dates.")
  }
  key_cols <- c(environment_col, date_col)
  duplicated_key <- duplicated(d[key_cols]) | duplicated(d[key_cols], fromLast = TRUE)
  n_duplicate_rows <- sum(duplicated_key)
  if (n_duplicate_rows && duplicate_action == "error")
    stop("Source '", source, "' has ", n_duplicate_rows,
         " rows involved in duplicate ", paste(key_cols, collapse = "-"),
         " keys.")
  if (n_duplicate_rows)
    d <- .aggregate_environment_duplicates(d, key_cols)

  missing_env <- setdiff(environments, unique(d[[environment_col]]))
  if (length(missing_env) && missing_action == "error")
    stop("Source '", source, "' is missing environments: ",
         paste(missing_env, collapse = ", "), ".")

  value_cols <- setdiff(names(d), key_cols)
  numeric_cols <- value_cols[vapply(d[value_cols], is.numeric, logical(1))]
  invalid <- matrix(FALSE, nrow(d), length(numeric_cols),
                    dimnames = list(NULL, numeric_cols))
  if (length(ranges)) {
    bad_names <- setdiff(names(ranges), numeric_cols)
    if (length(bad_names))
      warning("Ignoring ranges for absent/non-numeric variables: ",
              paste(bad_names, collapse = ", "), ".", call. = FALSE)
    for (nm in intersect(names(ranges), numeric_cols)) {
      rr <- as.numeric(ranges[[nm]])
      if (length(rr) != 2L || any(!is.finite(rr)) || rr[1] > rr[2])
        stop("Range for `", nm, "` must be finite c(lower, upper).")
      bad <- !is.na(d[[nm]]) & (d[[nm]] < rr[1] | d[[nm]] > rr[2])
      invalid[, nm] <- bad
      d[[nm]][bad] <- NA_real_
    }
  }

  originally_missing <- lapply(numeric_cols, function(nm) is.na(d[[nm]]))
  names(originally_missing) <- numeric_cols
  n_imputed <- stats::setNames(integer(length(numeric_cols)), numeric_cols)
  imputation_rows <- list()
  if (missing_action == "impute" && impute != "none") {
    for (nm in numeric_cols) {
      miss0 <- is.na(d[[nm]])
      if (any(miss0) && impute == "linear" && !is.null(date_col)) {
        for (ee in unique(d[[environment_col]])) {
          ii <- which(d[[environment_col]] == ee)
          oo <- order(d[[date_col]][ii])
          jj <- ii[oo]
          ok <- is.finite(d[[nm]][jj])
          if (sum(ok) >= 2L) {
            xp <- as.numeric(d[[date_col]][jj][ok])
            d[[nm]][jj] <- stats::approx(
              x = xp, y = d[[nm]][jj][ok],
              xout = as.numeric(d[[date_col]][jj]),
              rule = 1, ties = mean
            )$y
          }
        }
      }
      fill <- stats::median(d[[nm]], na.rm = TRUE)
      if (!is.finite(fill)) fill <- 0
      still <- is.na(d[[nm]])
      d[[nm]][still] <- fill
      changed <- miss0 & is.finite(d[[nm]])
      n_imputed[nm] <- sum(changed)
      if (any(changed)) {
        imputation_rows[[length(imputation_rows) + 1L]] <- data.frame(
          source = source,
          environment = d[[environment_col]][changed],
          date = if (is.null(date_col)) as.Date(NA) else d[[date_col]][changed],
          variable = nm,
          method = if (impute == "linear" && !is.null(date_col))
            "within-environment linear, then global median" else "global median",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (add_missing_indicators) {
    for (nm in numeric_cols) {
      flag <- originally_missing[[nm]]
      if (any(flag))
        d[[paste0(nm, "__was_missing")]] <- as.integer(flag)
    }
  }

  expected <- .environment_expected_counts(
    d, environments, environment_col, date_col, expected_dates
  )
  audit <- if (!length(numeric_cols)) {
    data.frame(
      source = character(), environment = character(), variable = character(),
      n_expected = integer(), n_observed = integer(), coverage = numeric(),
      n_invalid = integer(), n_imputed = integer(), status = character(),
      stringsAsFactors = FALSE
    )
  } else do.call(rbind, lapply(environments, function(ee) {
      ii <- d[[environment_col]] == ee
      do.call(rbind, lapply(numeric_cols, function(nm) {
        n_obs <- sum(!originally_missing[[nm]][ii])
        n_exp <- expected[[ee]]
        coverage <- if (n_exp > 0) n_obs / n_exp else 0
        data.frame(
          source = source, environment = ee, variable = nm,
          n_expected = n_exp, n_observed = n_obs,
          coverage = min(1, coverage),
          n_invalid = sum(invalid[ii, nm]),
          n_imputed = sum(originally_missing[[nm]][ii] & !is.na(d[[nm]][ii])),
          status = if (coverage < min_coverage) "below_threshold" else "ok",
          stringsAsFactors = FALSE
        )
      }))
    }))
  low <- audit$status == "below_threshold"
  if (any(low)) {
    msg <- paste0(
      "Source '", source, "' has ", sum(low),
      " environment-variable combinations below min_coverage = ",
      format(min_coverage), "."
    )
    if (missing_action == "error") stop(msg) else warning(msg, call. = FALSE)
  }
  if (missing_action == "warn" && anyNA(d[numeric_cols]))
    warning("Source '", source, "' retains missing numeric values.",
            call. = FALSE)

  if (!is.null(date_col))
    d <- d[order(match(d[[environment_col]], environments), d[[date_col]]),
           , drop = FALSE]
  provenance <- data.frame(
    source = source,
    rows_output = nrow(d),
    environments_expected = length(environments),
    environments_observed = length(unique(d[[environment_col]])),
    variables_numeric = length(numeric_cols),
    duplicate_rows = n_duplicate_rows,
    invalid_values = sum(invalid),
    imputed_values = sum(n_imputed),
    minimum_coverage = min_coverage,
    missing_action = missing_action,
    imputation = impute,
    stringsAsFactors = FALSE
  )
  imp <- if (length(imputation_rows)) do.call(rbind, imputation_rows) else
    data.frame(source = character(), environment = character(),
               date = as.Date(character()), variable = character(),
               method = character(), stringsAsFactors = FALSE)
  rownames(d) <- NULL
  list(data = d, audit = audit, provenance = provenance, imputation = imp)
}

.aggregate_environment_duplicates <- function(d, keys) {
  id <- interaction(d[keys], drop = TRUE, lex.order = TRUE)
  spl <- split(seq_len(nrow(d)), id)
  rows <- lapply(spl, function(ii) {
    z <- d[ii, , drop = FALSE]
    out <- z[1L, , drop = FALSE]
    for (nm in setdiff(names(d), keys)) {
      x <- z[[nm]]
      if (is.numeric(x)) {
        out[[nm]] <- if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      } else {
        ok <- which(!is.na(x) & nzchar(as.character(x)))
        out[[nm]] <- if (length(ok)) x[ok[1L]] else x[1L]
      }
    }
    out
  })
  do.call(rbind, rows)
}

.environment_expected_counts <- function(d, environments, environment_col,
                                         date_col, expected_dates) {
  if (is.null(date_col))
    return(stats::setNames(rep(1L, length(environments)), environments))
  if (is.null(expected_dates)) {
    all_dates <- sort(unique(d[[date_col]]))
    return(stats::setNames(rep(length(all_dates), length(environments)),
                           environments))
  }
  if (is.data.frame(expected_dates)) {
    if (!all(c(environment_col, date_col) %in% names(expected_dates)))
      stop("`expected_dates` data frame needs `", environment_col,
           "` and `", date_col, "`.")
    tab <- table(factor(expected_dates[[environment_col]], levels = environments))
    return(stats::setNames(as.integer(tab), environments))
  }
  stats::setNames(rep(length(unique(as.Date(expected_dates))),
                      length(environments)), environments)
}
