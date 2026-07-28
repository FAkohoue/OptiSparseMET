#' Build one environment relationship matrix for exploratory selection
#'
#' @description
#' [select_environments()] needs an environment-by-environment similarity matrix
#' `D`. `build_environment_relationship()` constructs one such matrix from
#' historical multi-environment trial responses or one matrix of environmental
#' covariates. The environmental route can include unseen sites for which
#' covariates, but no phenotypes, are available. It describes exposure
#' similarity; it is not automatically a trait-specific genetic covariance.
#' Use [build_environment_kernels()] when weather, soil, management, geography,
#' and their interactions must remain separate.
#'
#' @param x Input data.
#'   \describe{
#'     \item{`source = "genetic_correlation"`}{A genotype-by-environment matrix
#'       of genetic values or best linear unbiased estimates (BLUEs; genotypes
#'       in rows, environments in
#'       columns). `D` is the correlation between environment columns.}
#'     \item{`source = "enviromic"`}{An environment-by-covariate matrix
#'       (environments in rows, standardised or raw covariates in columns).}
#'   }
#' @param source Either `"genetic_correlation"` or `"enviromic"`.
#' @param kernel For `source = "enviromic"`: `"correlation"` (correlation of
#'   environments across standardised covariates) or `"gaussian"` (Gaussian
#'   kernel on covariate distance).
#'
#'   **Interactions and the choice of kernel.** The Gaussian radial basis
#'   function kernel is non-linear, but that does not identify a biologically
#'   interpretable weather-by-soil or weather-by-management contribution.
#'   Explicit scientific hypotheses should be represented by separate
#'   cross-modality kernels from [build_environment_kernels()] or
#'   [add_environment_kernel_interactions()]. The `"correlation"` kernel is
#'   linear and contains no implicit interaction feature space.
#' @param bandwidth Optional Gaussian-kernel bandwidth; defaults to the median
#'   pairwise distance.
#' @param variables For `source = "enviromic"`: optional character vector (column
#'   names) or integer indices selecting which covariates to use. This is how a
#'   breeder restricts the kinship to the *limiting factors* they judge relevant
#'   (e.g. only `mean_temp`, `total_precip`, `clay`) instead of every assembled
#'   covariate. The valid names are exactly `colnames(x)`; call
#'   [enviromic_variable_catalog()] to see every fetchable variable with its
#'   meaning and units. Defaults to all columns.
#' @param weights For `source = "enviromic"`: optional numeric weights for the
#'   selected covariates (a named vector is matched by name, otherwise taken in
#'   column order). Higher weights make a variable count more toward the
#'   similarity. Weighting is
#'   applied as \eqn{\sqrt{w}} on the standardised columns, so both the
#'   correlation and Gaussian kernels honour it. The default gives equal weight
#'   to the selected columns. Non-equal weights require external agronomic or
#'   empirical justification. They must not be used to manufacture
#'   trait-specific modality weights when historical MET responses are absent.
#' @param geo Logical. If `TRUE`, `x` is treated as site coordinates and
#'   great-circle (haversine) distance is used instead of Euclidean -- correct
#'   for GPS latitude/longitude at continental scale. The columns are taken as
#'   `latitude`, `longitude` (by those names if present, else the first two).
#'   Only used with `source = "enviromic", kernel = "gaussian"`.
#' @param missing Missing-covariate policy for `source = "enviromic"`:
#'   `"error"` (default) or explicit column-median imputation. Prefer
#'   [qc_environmental_data()] when an audit ledger is required.
#'
#' @return A square environment relationship matrix with row/column names, ready
#'   to pass to [select_environments()].
#'
#' @seealso [select_environments()].
#' @export
build_environment_relationship <- function(x,
                                           source = c("genetic_correlation",
                                                      "enviromic"),
                                           kernel = c("correlation", "gaussian"),
                                           bandwidth = NULL, variables = NULL,
                                           weights = NULL, geo = FALSE,
                                           missing = c("error", "median")) {
  source <- match.arg(source)
  kernel <- match.arg(kernel)
  missing <- match.arg(missing)
  x <- as.matrix(x)

  if (source == "genetic_correlation") {
    if (is.null(colnames(x)))
      colnames(x) <- paste0("E", seq_len(ncol(x)))
    D <- stats::cor(x, use = "pairwise.complete.obs")
    if (anyNA(diag(D)))
      stop("At least one environment has no estimable variance for a genetic ",
           "correlation.")
    if (anyNA(D))
      warning("Some pairwise genetic correlations are not estimable and remain ",
              "missing.", call. = FALSE)
    return(D)
  }

  # enviromic: environments in rows, covariates in columns
  if (is.null(rownames(x)))
    rownames(x) <- paste0("E", seq_len(nrow(x)))
  if (is.null(colnames(x)))
    colnames(x) <- paste0("V", seq_len(ncol(x)))

  # Restrict to the user-chosen limiting variables (by name or index).
  if (!is.null(variables) && !geo) {
    if (is.character(variables)) {
      miss <- setdiff(variables, colnames(x))
      if (length(miss))
        stop("`variables` not found in covariate columns: ",
             paste(miss, collapse = ", "), ".\nAvailable columns are: ",
             paste(colnames(x), collapse = ", "),
             ".\nSee enviromic_variable_catalog() for what each name means.")
      x <- x[, variables, drop = FALSE]
    } else {
      x <- x[, variables, drop = FALSE]
    }
    if (ncol(x) < 1L) stop("`variables` selected no covariates.")
  }

  if (geo) {
    cn <- tolower(colnames(x))
    lat_i <- if ("latitude" %in% cn) which(cn == "latitude")[1] else 1L
    lon_i <- if ("longitude" %in% cn) which(cn == "longitude")[1] else 2L
    xy <- x[, c(lat_i, lon_i), drop = FALSE]
    if (any(!is.finite(xy)))
      stop("Geographic coordinates must be finite; they are not imputed.")
    dd <- .haversine_matrix(xy[, 1L], xy[, 2L])
    h  <- if (is.null(bandwidth)) stats::median(dd[upper.tri(dd)]) else bandwidth
    if (length(h) != 1L)
      stop("`bandwidth` must be one positive number.")
    if (!is.finite(h) || h <= 0) h <- 1
    D  <- exp(-(dd^2) / (2 * h^2))
    dimnames(D) <- list(rownames(x), rownames(x))
    attr(D, "bandwidth") <- h
    return(D)
  }

  if (any(is.infinite(x)))
    stop("Enviromic covariates contain infinite values.")
  if (anyNA(x)) {
    if (missing == "error")
      stop("Enviromic covariates contain missing values. Use ",
           "`qc_environmental_data()` or set `missing = 'median'` explicitly.")
    for (j in seq_len(ncol(x))) {
      fill <- stats::median(x[, j], na.rm = TRUE)
      if (!is.finite(fill))
        stop("Covariate '", colnames(x)[j],
             "' has no finite value for median imputation.")
      x[is.na(x[, j]), j] <- fill
    }
  }

  # Validate weights before removing constant variables so unnamed weights
  # remain tied to the user's original selected column order.
  w <- NULL
  if (!is.null(weights)) {
    w <- weights
    if (!is.null(names(w))) {
      miss <- setdiff(colnames(x), names(w))
      if (length(miss))
        stop("`weights` is missing entries for covariates: ",
             paste(miss, collapse = ", "), ".")
      w <- w[colnames(x)]
    }
    if (length(w) != ncol(x))
      stop("`weights` must have one value per selected covariate (length ",
           ncol(x), ").")
    if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0)
      stop("`weights` must be finite, non-negative, and not all zero.")
  }
  informative <- vapply(seq_len(ncol(x)), function(j) {
    s <- stats::sd(x[, j])
    is.finite(s) && s > 0
  }, logical(1))
  dropped <- colnames(x)[!informative]
  if (!any(informative))
    stop("No non-constant enviromic covariates are available.")
  x <- x[, informative, drop = FALSE]
  if (!is.null(w)) w <- w[informative]
  Z <- scale(x)

  # Emphasise the most limiting covariates: sqrt(w) column scaling makes the
  # weighted Euclidean/covariance geometry honour the weights for both kernels.
  if (!is.null(w)) {
    Z <- sweep(Z, 2L, sqrt(w), `*`)
  }

  if (kernel == "correlation") {
    row_sd <- apply(Z, 1L, stats::sd)
    if (any(!is.finite(row_sd) | row_sd <= 0))
      stop("The correlation kernel is undefined for an environment with no ",
           "within-profile variation; use the Gaussian kernel or revise the ",
           "covariate block.")
    D <- stats::cor(t(Z))
    if (anyNA(D))
      stop("The correlation kernel is undefined for an environment with no ",
           "within-profile variation; use the Gaussian kernel or revise the ",
           "covariate block.")
  } else {
    dd <- as.matrix(stats::dist(Z))
    h  <- if (is.null(bandwidth)) stats::median(dd[upper.tri(dd)]) else bandwidth
    if (length(h) != 1L)
      stop("`bandwidth` must be one positive number.")
    if (!is.finite(h) || h <= 0) h <- 1
    D  <- exp(-(dd^2) / (2 * h^2))
    attr(D, "bandwidth") <- h
  }
  dimnames(D) <- list(rownames(x), rownames(x))
  attr(D, "dropped_constant_variables") <- dropped
  D
}


# Great-circle (haversine) distance matrix in km from latitude/longitude vectors.
.haversine_matrix <- function(lat, lon, R = 6371) {
  lat <- as.numeric(lat); lon <- as.numeric(lon); n <- length(lat)
  rad <- pi / 180
  phi <- lat * rad; lam <- lon * rad
  D <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) if (i < j) {
    dphi <- phi[j] - phi[i]; dl <- lam[j] - lam[i]
    a <- sin(dphi / 2)^2 + cos(phi[i]) * cos(phi[j]) * sin(dl / 2)^2
    D[i, j] <- D[j, i] <- 2 * R * asin(min(1, sqrt(a)))
  }
  D
}
