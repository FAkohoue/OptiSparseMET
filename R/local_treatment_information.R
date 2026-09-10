#' Treatment information from a realised field layout
#'
#' @description
#' Builds the treatment-information matrix contributed by one realised field
#' layout after absorbing the intercept, checks, fixed blocking factors, and
#' optional random nuisance effects. Unlike a single efficiency multiplier,
#' the returned matrix retains which treatment contrasts the layout measures
#' well or poorly and can be passed directly to [met_information()].
#'
#' @param field_book A data frame containing `Treatment`, `Row`, and `Column`.
#' @param treatments Target treatments whose information is required. Defaults
#'   to all non-check treatments in the field book.
#' @param check_treatments Treatments modelled as fixed checks and excluded from
#'   the target information matrix.
#' @param fixed_effects Field-book columns absorbed as fixed nuisance factors.
#'   Missing names are ignored. Default `"Block"`.
#' @param random_effect_variances Named list or numeric vector mapping
#'   field-book factor columns to positive random-effect variances. Defaults to
#'   empty (no random nuisance effects).
#' @param residual_structure,rho_row,rho_col,nugget,kernel_range,matern_nu
#'   Spatial residual model, as in [met_evaluate_famoptg_efficiency()].
#' @param sigma_e2 Positive residual variance.
#' @return A symmetric positive-semidefinite treatment-information matrix with
#'   rows and columns in `treatments` order.
#' @export
local_treatment_information <- function(
    field_book, treatments = NULL, check_treatments = character(0),
    fixed_effects = "Block", random_effect_variances = list(),
    residual_structure = c("IID", "AR1", "AR1xAR1", "AR1xAR1_nugget",
                           "exponential", "gaussian", "matern"),
    rho_row = 0, rho_col = 0, nugget = 0, kernel_range = NULL,
    matern_nu = 1.5, sigma_e2 = 1) {
  residual_structure <- match.arg(residual_structure)
  if (!is.data.frame(field_book) ||
      !all(c("Treatment", "Row", "Column") %in% names(field_book)))
    stop("`field_book` must contain Treatment, Row, and Column columns.")
  keep <- !is.na(field_book$Treatment) & nzchar(as.character(field_book$Treatment))
  fb <- field_book[keep, , drop = FALSE]
  if (!nrow(fb)) stop("`field_book` contains no occupied plots.")
  trt_obs <- as.character(fb$Treatment)
  check_treatments <- unique(as.character(check_treatments))
  if (is.null(treatments))
    treatments <- setdiff(unique(trt_obs), check_treatments)
  treatments <- unique(as.character(treatments))
  if (!length(treatments) || anyNA(treatments) || any(!nzchar(treatments)))
    stop("`treatments` must contain at least one non-empty treatment ID.")
  if (length(intersect(treatments, check_treatments)))
    stop("`treatments` and `check_treatments` must not overlap.")
  if (!is.numeric(sigma_e2) || length(sigma_e2) != 1L ||
      !is.finite(sigma_e2) || sigma_e2 <= 0)
    stop("`sigma_e2` must be one finite positive number.")

  n <- nrow(fb); p <- length(treatments)
  zcol <- match(trt_obs, treatments)
  zrows <- which(!is.na(zcol))
  Z <- Matrix::sparseMatrix(i = zrows, j = zcol[zrows], x = 1,
                            dims = c(n, p),
                            dimnames = list(NULL, treatments))

  if (residual_structure == "IID") {
    Q <- Matrix::Diagonal(n, 1 / sigma_e2)
  } else {
    Q <- .build_residual_precision(
      fb$Row, fb$Column, max(as.integer(fb$Row)), max(as.integer(fb$Column)),
      residual_structure, rho_row, rho_col, sigma_e2,
      nugget = nugget, kernel_range = kernel_range, matern_nu = matern_nu)
  }

  X <- Matrix::Matrix(rep(1, n), ncol = 1, sparse = TRUE)
  colnames(X) <- "(Intercept)"
  add_fixed_factor <- function(X, values, prefix) {
    f <- factor(values)
    if (nlevels(f) <= 1L) return(X)
    mm <- stats::model.matrix(~ f)[, -1L, drop = FALSE]
    colnames(mm) <- paste0(prefix, "_", levels(f)[-1L])
    cbind(X, Matrix::Matrix(mm, sparse = TRUE))
  }
  for (nm in intersect(as.character(fixed_effects), names(fb)))
    X <- add_fixed_factor(X, fb[[nm]], nm)
  if (length(check_treatments)) {
    for (chk in check_treatments) {
      v <- as.numeric(trt_obs == chk)
      if (any(v)) X <- cbind(X, Matrix::Matrix(v, ncol = 1, sparse = TRUE))
    }
  }

  if (is.numeric(random_effect_variances) && !is.list(random_effect_variances))
    random_effect_variances <- as.list(random_effect_variances)
  if (!is.list(random_effect_variances) ||
      (length(random_effect_variances) &&
       (is.null(names(random_effect_variances)) ||
        any(names(random_effect_variances) == ""))))
    stop("`random_effect_variances` must be a named list or vector.")
  W <- NULL; random_prec <- numeric(0)
  for (nm in names(random_effect_variances)) {
    vv <- random_effect_variances[[nm]]
    if (!is.numeric(vv) || length(vv) != 1L || !is.finite(vv) || vv <= 0)
      stop("Every random-effect variance must be one finite positive number.")
    if (!nm %in% names(fb))
      stop("Random-effect column `", nm, "` is absent from `field_book`.")
    inc <- stats::model.matrix(~ factor(fb[[nm]]) - 1)
    W <- if (is.null(W)) Matrix::Matrix(inc, sparse = TRUE) else
      cbind(W, Matrix::Matrix(inc, sparse = TRUE))
    random_prec <- c(random_prec, rep(1 / vv, ncol(inc)))
  }

  N <- if (is.null(W)) X else cbind(X, W)
  penalty <- c(rep(0, ncol(X)), random_prec)
  NtQN <- Matrix::crossprod(N, Q %*% N) + Matrix::Diagonal(x = penalty)
  NtQZ <- Matrix::crossprod(N, Q %*% Z)
  absorbed <- .solve_C(NtQN, NtQZ)
  info <- as.matrix(Matrix::crossprod(Z, Q %*% Z) -
                      t(NtQZ) %*% absorbed)
  # Numerical symmetrisation and removal of negligible negative eigenvalues.
  info <- (info + t(info)) / 2
  ee <- eigen(info, symmetric = TRUE)
  tol <- 1e-10 * max(abs(ee$values), 1)
  if (min(ee$values) < -100 * tol)
    stop("The absorbed treatment-information matrix is materially indefinite.")
  info <- ee$vectors %*% (t(ee$vectors) * pmax(ee$values, 0))
  dimnames(info) <- list(treatments, treatments)
  attr(info, "model") <- list(
    residual_structure = residual_structure, sigma_e2 = sigma_e2,
    fixed_effects = intersect(as.character(fixed_effects), names(fb)),
    random_effect_variances = random_effect_variances)
  info
}
