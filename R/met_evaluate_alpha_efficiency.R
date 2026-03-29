# ==============================================================================
# met_evaluate_alpha_efficiency.R
#
# OptiSparseMET version of evaluate_alpha_efficiency() from OptiDesign.
# The met_ prefix prevents namespace conflicts when both packages are loaded.
# Function body is identical to evaluate_alpha_efficiency(). Only the
# name changes, and varcomp uses sigma_rep2 + sigma_ib2 (Rep/IBlock nesting).
# ==============================================================================

#' Evaluate the statistical efficiency of an alpha-lattice design
#'
#' @description
#' `met_evaluate_alpha_efficiency()` is the OptiSparseMET version of
#' `evaluate_alpha_efficiency()` from the OptiDesign package. The `met_`
#' prefix avoids namespace conflicts when both packages are loaded
#' simultaneously. All arguments, return values, and internal logic are
#' identical to `evaluate_alpha_efficiency()` in OptiDesign.
#'
#' `met_evaluate_alpha_efficiency()` takes a `field_book` produced by
#' [met_alpha_rc_stream()] and computes optimality criteria for the design
#' under a user-specified mixed model. It is fully decoupled from design
#' construction, so a single design can be evaluated multiple times under
#' different model assumptions without rebuilding the layout.
#'
#' @details
#' ## Mixed model
#'
#' The model is:
#' \deqn{y = X\beta + Zu + e}
#'
#' **Fixed part** \eqn{X}: intercept, optional check fixed effects
#' (`check_as_fixed = TRUE`), entry treatment fixed effects when
#' `treatment_effect = "fixed"`.
#'
#' **Random part** \eqn{Z}: replicate, incomplete block within replicate, row,
#' column, and - when `treatment_effect = "random"` - entry treatment effects.
#'
#' **Residual structure** \eqn{R = \sigma_e^2 \Sigma}:
#' - `"IID"`: \eqn{\Sigma = I}
#' - `"AR1"`: row-only AR1,
#'   \eqn{Q_\text{AR1}(\rho_\text{row}) \otimes I_\text{cols}}
#' - `"AR1xAR1"`: separable row x column AR1,
#'   \eqn{Q_\text{AR1}(\rho_\text{col}) \otimes Q_\text{AR1}(\rho_\text{row})}
#'
#' AR1 precision matrices are tridiagonal sparse matrices. For an
#' \eqn{n}-dimensional AR1 process with parameter \eqn{\rho}, let
#' \eqn{a = 1/(1-\rho^2)}:
#' - Interior diagonal: \eqn{(1+\rho^2) \times a}
#' - Edge diagonal (positions 1 and \eqn{n}): \eqn{a}
#' - Off-diagonal: \eqn{-\rho \times a}
#'
#' ## Mixed model coefficient matrix
#'
#' \deqn{C = \begin{pmatrix} X^\top Q X & X^\top Q Z \\ Z^\top Q X &
#'   Z^\top Q Z + G^{-1} \end{pmatrix}}
#'
#' where \eqn{Q = R^{-1}} is the residual precision matrix.
#'
#' ## Random-effect structure G
#'
#' - `"IID"`: \eqn{G^{-1}_\text{entry} = \sigma_g^{-2} I}
#' - `"GBLUP"` / `"PBLUP"`: \eqn{G^{-1}_\text{entry} = \sigma_g^{-2} K^{-1}},
#'   computed via sparse Cholesky; falls back to Moore-Penrose pseudoinverse
#'   for near-singular `K`.
#'
#' ## Optimality criteria
#'
#' ### Fixed treatment effects
#'
#' **A-criterion** (lower = better):
#' \deqn{A_\text{criterion} = \bar{v}_\text{diff} =
#'   \frac{2}{p(p-1)} \sum_{i<j} \text{Var}(\hat{\tau}_i - \hat{\tau}_j)}
#'
#' **D-criterion** (lower = better):
#' \deqn{D_\text{criterion} = \exp\!\left(\frac{\log\det(HVH)}{p-1}\right)}
#' where \eqn{H = I_p - p^{-1}J_p} is the centering matrix and \eqn{V} is
#' the treatment variance-covariance submatrix.
#'
#' **Efficiency forms** (higher = better):
#' \deqn{A_\text{efficiency} = 1 / A_\text{criterion}, \quad
#'       D_\text{efficiency} = 1 / D_\text{criterion}}
#'
#' ### Random treatment effects (genomic prediction)
#'
#' **PEV-criterion** (lower = better): mean prediction error variance across
#' all entry lines.
#' \deqn{\text{PEV}_\text{criterion} = \frac{1}{p}\sum_{i=1}^{p}
#'   \text{Var}(\hat{u}_i - u_i)}
#'
#' **CDmean** (higher = better): mean coefficient of determination for genomic
#' breeding value prediction (Rincent et al. 2012, *Genetics* 192:715-728).
#' \deqn{\text{CDmean} = 1 - \frac{\text{PEV}_\text{criterion}}{\sigma_g^2}}
#'
#' CDmean measures the proportion of genetic variance explained by GEBV
#' predictions on average across lines. Values close to 1 indicate highly
#' reliable predictions; values close to 0 indicate near-uninformative
#' predictions. CDmean is only defined when `treatment_effect = "random"`.
#'
#' **CD per line**: when the number of treatments does not exceed
#' `eff_full_max`, a per-line coefficient of determination vector is also
#' returned:
#' \deqn{CD_i = 1 - \text{PEV}_i / \sigma_g^2}
#'
#' ## Large-design approximation
#'
#' When the number of treatments exceeds `eff_full_max`, a Hutchinson
#' stochastic trace estimator with `eff_trace_samples` Rademacher vectors
#' replaces exact submatrix extraction. The `mode` field carries the suffix
#' `_APPROX` and `D_criterion`, `D_efficiency`, and `CD_per_line` are `NA`.
#'
#' @param field_book Data frame produced by [met_alpha_rc_stream()]. Must
#'   contain columns: `Plot`, `Row`, `Column`, `Rep`, `IBlock`, `BlockInRep`,
#'   `Treatment`, `Check`.
#'
#' @param n_rows Positive integer. Number of field rows (must match the
#'   original design).
#'
#' @param n_cols Positive integer. Number of field columns (must match the
#'   original design).
#'
#' @param check_treatments Character vector of check treatment identifiers
#'   (must match those used in [met_alpha_rc_stream()]).
#'
#' @param treatment_effect Character. Whether entry treatments are modelled as
#'   `"fixed"` (BLUE-based A and D criteria) or `"random"` (PEV-based
#'   criterion and CDmean). Note: `"fixed"` does not produce CDmean; `"random"`
#'   does not produce D-criterion.
#'
#' @param prediction_type Character. Random-effect prediction model for entries.
#'   `"IID"` assumes independent entries (\eqn{G^{-1} = \sigma_g^{-2} I}).
#'   `"GBLUP"` and `"PBLUP"` use the supplied `K` matrix
#'   (\eqn{G^{-1} = \sigma_g^{-2} K^{-1}}). `"none"` is invalid when
#'   `treatment_effect = "random"` and raises an error.
#'
#' @param check_as_fixed Logical. If `TRUE`, checks are included as fixed
#'   effects in the model.
#'
#' @param residual_structure Character. Residual covariance structure.
#'   `"IID"` (independence), `"AR1"` (row AR1 only), or `"AR1xAR1"`
#'   (separable row x column AR1).
#'
#' @param rho_row Numeric in \eqn{(-1, 1)}. AR1 autocorrelation parameter
#'   along rows. Used when `residual_structure %in% c("AR1", "AR1xAR1")`.
#'
#' @param rho_col Numeric in \eqn{(-1, 1)}. AR1 autocorrelation parameter
#'   along columns. Used only when `residual_structure = "AR1xAR1"`.
#'
#' @param varcomp Named list of variance components. Required components:
#'   `sigma_e2` (residual), `sigma_g2` (entry genetic), `sigma_rep2`
#'   (replicate), `sigma_ib2` (incomplete block), `sigma_r2` (row), `sigma_c2`
#'   (column). All components must be present even if the corresponding effect
#'   is absent from the model. `sigma_g2` is used as the denominator of CDmean
#'   when `treatment_effect = "random"`. Note: this function uses `sigma_rep2`
#'   and `sigma_ib2` for the Rep and IBlock(Rep) random effects, which differs
#'   from [met_evaluate_famoptg_efficiency()] which uses `sigma_b2`.
#'
#' @param K Optional square numeric matrix with rownames and colnames. Required
#'   when `prediction_type %in% c("GBLUP", "PBLUP")`.
#'
#' @param line_id_map Optional data frame with columns `Treatment` and
#'   `LineID`. Required when treatment labels do not match `rownames(K)`.
#'
#' @param spatial_engine Character. Computational engine for matrix operations.
#'   `"auto"` uses `"dense"` when used plots \eqn{\leq} `dense_max_n`,
#'   otherwise `"sparse"`.
#'
#' @param dense_max_n Integer. Threshold for `spatial_engine = "auto"`.
#'
#' @param eff_trace_samples Positive integer. Number of Rademacher vectors for
#'   the Hutchinson stochastic trace estimator (large designs only).
#'
#' @param eff_full_max Positive integer. Maximum number of treatments for exact
#'   \eqn{C^{-1}} submatrix extraction. Designs with more treatments use the
#'   stochastic estimator (`mode` suffix `_APPROX`).
#'
#' @return A named list containing model metadata and criterion values:
#' \describe{
#'   \item{`model`}{Character. Model string.}
#'   \item{`treatment_effect`}{Character. As supplied.}
#'   \item{`prediction_type`}{Character or `NA`.}
#'   \item{`residual_structure_requested`}{Character. As supplied.}
#'   \item{`residual_structure_used`}{Character. As resolved.}
#'   \item{`spatial_engine_used`}{Character. `"dense"` or `"sparse"`.}
#'   \item{`mode`}{Character. Computation mode: `"FIXED_TREATMENT_BLUE_CONTRAST"`,
#'     `"FIXED_TREATMENT_BLUE_APPROX"`, `"RANDOM_TREATMENT_PEV"`, or
#'     `"RANDOM_TREATMENT_PEV_APPROX"`.}
#'   \item{`n_trt`}{Integer. Number of treatment columns evaluated.}
#'   \item{`n_contrasts`}{Integer. `n_trt - 1`. Fixed mode only.}
#'   \item{`A_criterion`}{Numeric. Mean pairwise contrast variance (fixed) or
#'     mean PEV (random). Lower is better.}
#'   \item{`D_criterion`}{Numeric or `NA`. Geometric mean of contrast
#'     covariance eigenvalues. Fixed full mode only. Lower is better.}
#'   \item{`A_efficiency`}{Numeric. `1 / A_criterion`. Higher is better.}
#'   \item{`D_efficiency`}{Numeric or `NA`. `1 / D_criterion`. Higher is
#'     better.}
#'   \item{`A`}{Numeric. Alias for `A_efficiency` (backward compatibility).}
#'   \item{`D`}{Numeric or `NA`. Alias for `D_efficiency` (backward
#'     compatibility).}
#'   \item{`mean_VarDiff`}{Numeric. Mean pairwise contrast variance. Fixed
#'     full mode only.}
#'   \item{`PEV_criterion`}{Numeric. Mean PEV. Random mode only.}
#'   \item{`mean_PEV`}{Numeric. Alias for `PEV_criterion`. Random mode only.}
#'   \item{`CDmean`}{Numeric. Mean coefficient of determination for GEBV
#'     prediction: \eqn{1 - \text{mean\_PEV} / \sigma_g^2}. Range \eqn{[0,1]}.
#'     Higher is better. Random mode only.}
#'   \item{`CD_per_line`}{Numeric vector or `NA`. Per-line coefficient of
#'     determination: \eqn{1 - \text{PEV}_i / \sigma_g^2}. Available in full
#'     mode only (`NA` in `_APPROX` mode). Random mode only.}
#' }
#'
#' @references
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P.,
#' ... & Moreau, L. (2012). Maximizing the reliability of genomic selection by
#' optimizing the calibration set of reference individuals: comparison of
#' methods in two diverse groups of maize inbreds (*Zea mays* L.).
#' *Genetics*, 192(2), 715-728.
#'
#' @seealso [met_alpha_rc_stream()] to construct the design.
#'   [met_optimize_alpha_rc()] to search for a criterion-optimal design.
#'   [met_evaluate_famoptg_efficiency()] for the equivalent function for
#'   [met_prep_famoptg()] designs (uses `sigma_b2` instead of `sigma_rep2`
#'   + `sigma_ib2`).
#'
#' @examples
#' design <- met_alpha_rc_stream(
#'   check_treatments = c("CHK1", "CHK2", "CHK3"),
#'   check_families   = c("CHECK", "CHECK", "CHECK"),
#'   entry_treatments = paste0("G", 1:167),
#'   entry_families   = rep(paste0("F", 1:7), length.out = 167),
#'   n_reps = 3, n_rows = 30, n_cols = 20,
#'   min_block_size = 19, max_block_size = 20
#' )
#'
#' ## Fixed-effect evaluation (BLUE contrasts, AR1xAR1 spatial)
#' eff_fixed <- met_evaluate_alpha_efficiency(
#'   field_book         = design$field_book,
#'   n_rows = 30, n_cols = 20,
#'   check_treatments   = c("CHK1", "CHK2", "CHK3"),
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row = 0.10, rho_col = 0.10
#' )
#' eff_fixed$A_criterion   # lower is better
#' eff_fixed$D_criterion   # lower is better
#' eff_fixed$A_efficiency  # higher is better
#'
#' ## Random-effect evaluation (GBLUP with CDmean)
#' \dontrun{
#' eff_gblup <- met_evaluate_alpha_efficiency(
#'   field_book         = design$field_book,
#'   n_rows = 30, n_cols = 20,
#'   check_treatments   = c("CHK1", "CHK2", "CHK3"),
#'   treatment_effect   = "random",
#'   prediction_type    = "GBLUP",
#'   K                  = my_kinship_matrix,
#'   varcomp            = list(sigma_g2 = 0.4, sigma_e2 = 0.6,
#'                             sigma_rep2 = 0.1, sigma_ib2 = 0.05,
#'                             sigma_r2 = 0.02, sigma_c2 = 0.02),
#'   residual_structure = "AR1xAR1",
#'   rho_row = 0.10, rho_col = 0.10
#' )
#' eff_gblup$CDmean       # mean prediction reliability [0, 1]
#' eff_gblup$CD_per_line  # per-line reliability vector
#' eff_gblup$mean_PEV     # mean prediction error variance
#' }
#'
#' @importFrom Matrix Diagonal sparseMatrix crossprod solve t Cholesky
#' @export
met_evaluate_alpha_efficiency <- function(
    field_book,
    n_rows,
    n_cols,
    check_treatments,

    # -- Model specification --------------------------------------------------
    treatment_effect   = c("random", "fixed"),
    prediction_type    = c("IID", "GBLUP", "PBLUP", "none"),
    check_as_fixed     = TRUE,
    residual_structure = c("IID", "AR1", "AR1xAR1"),
    rho_row            = 0,
    rho_col            = 0,

    # -- Variance components --------------------------------------------------
    varcomp = list(
      sigma_e2   = 1,
      sigma_g2   = 1,
      sigma_rep2 = 1,
      sigma_ib2  = 1,
      sigma_r2   = 1,
      sigma_c2   = 1
    ),

    # -- Genomic / pedigree relationship matrix -------------------------------
    K            = NULL,
    line_id_map  = NULL,

    # -- Computation controls -------------------------------------------------
    spatial_engine    = c("auto", "sparse", "dense"),
    dense_max_n       = 5000,
    eff_trace_samples = 80,
    eff_full_max      = 400
) {

  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' is required.")

  # -- Argument matching -------------------------------------------------------
  treatment_effect   <- match.arg(treatment_effect)
  prediction_type    <- match.arg(prediction_type)
  residual_structure <- match.arg(residual_structure)
  spatial_engine     <- match.arg(spatial_engine)

  # -- Pre-flight guard --------------------------------------------------------
  # Efficiency is only well-defined when there is something to estimate.
  if (treatment_effect == "random" && prediction_type == "none") {
    stop(paste0(
      "Cannot compute efficiency: treatment_effect = 'random' and prediction_type = 'none'.\n",
      "Either set treatment_effect = 'fixed', or choose prediction_type in ",
      "c('IID', 'GBLUP', 'PBLUP')."
    ))
  }

  if (prediction_type %in% c("GBLUP", "PBLUP") && is.null(K)) {
    stop("K must be provided when prediction_type is 'GBLUP' or 'PBLUP'.")
  }

  if (!is.list(varcomp) ||
      !all(c("sigma_e2", "sigma_g2", "sigma_rep2",
             "sigma_ib2", "sigma_r2", "sigma_c2") %in% names(varcomp)))
    stop("varcomp must contain: sigma_e2, sigma_g2, sigma_rep2, sigma_ib2, sigma_r2, sigma_c2.")

  if (!is.numeric(rho_row) || length(rho_row) != 1 || abs(rho_row) >= 1)
    stop("rho_row must satisfy |rho_row| < 1.")
  if (!is.numeric(rho_col) || length(rho_col) != 1 || abs(rho_col) >= 1)
    stop("rho_col must satisfy |rho_col| < 1.")

  required_cols <- c("Plot", "Row", "Column", "Rep", "IBlock",
                     "BlockInRep", "Treatment", "Check")
  missing_cols  <- setdiff(required_cols, names(field_book))
  if (length(missing_cols) > 0)
    stop(paste0("field_book is missing columns: ", paste(missing_cols, collapse = ", ")))

  # -- Unpack varcomp ----------------------------------------------------------
  sigma_e2   <- varcomp$sigma_e2
  sigma_g2   <- varcomp$sigma_g2
  sigma_rep2 <- varcomp$sigma_rep2
  sigma_ib2  <- varcomp$sigma_ib2
  sigma_r2   <- varcomp$sigma_r2
  sigma_c2   <- varcomp$sigma_c2

  # -- Working subset: rows with assigned treatments ---------------------------
  fb  <- field_book[!is.na(field_book$Treatment), , drop = FALSE]
  nn  <- nrow(fb)
  if (nn == 0) stop("field_book contains no assigned treatments.")

  trt      <- as.character(fb$Treatment)
  is_check <- as.logical(fb$Check)
  blk      <- paste0("IBlock_", fb$IBlock)
  repfac   <- paste0("Rep_",    fb$Rep)
  rw       <- fb$Row
  cl       <- fb$Column

  # -- Spatial engine ----------------------------------------------------------
  spatial_engine_use <- if (spatial_engine == "auto") {
    if (nn <= dense_max_n) "dense" else "sparse"
  } else {
    spatial_engine
  }

  # -- Residual precision matrix Q ---------------------------------------------
  residual_use <- residual_structure

  if (residual_use == "IID") {
    Q <- Matrix::Diagonal(nn, 1 / sigma_e2)

  } else if (residual_use == "AR1") {
    Qrow       <- .ar1_precision_sparse(n_rows, rho_row)
    grid_index <- (as.integer(rw) - 1L) * n_cols + as.integer(cl)
    Qfull      <- Matrix::kronecker(Matrix::Diagonal(n_cols, 1), Qrow) * (1 / sigma_e2)
    Q          <- Qfull[grid_index, grid_index, drop = FALSE]

  } else {
    # AR1xAR1
    Qrow       <- .ar1_precision_sparse(n_rows, rho_row)
    Qcol       <- .ar1_precision_sparse(n_cols, rho_col)
    Qgrid      <- Matrix::kronecker(Qcol, Qrow) * (1 / sigma_e2)
    grid_index <- (as.integer(rw) - 1L) * n_cols + as.integer(cl)
    Q          <- Qgrid[grid_index, grid_index, drop = FALSE]
  }

  # -- Fixed effects matrix X --------------------------------------------------
  X_int         <- Matrix::Matrix(rep(1, nn), ncol = 1, sparse = TRUE)
  colnames(X_int) <- "(Intercept)"
  X             <- X_int

  if (isTRUE(check_as_fixed)) {
    chk_vec <- ifelse(is_check, trt, NA)
    chk_inc <- .make_sparse_incidence(chk_vec)
    if (ncol(chk_inc$M) > 0) {
      keep <- intersect(check_treatments, colnames(chk_inc$M))
      chkM <- chk_inc$M[, keep, drop = FALSE]
      chkM <- chkM[, check_treatments[check_treatments %in% keep], drop = FALSE]
      colnames(chkM) <- paste0("Check_", colnames(chkM))
      X <- cbind(X, chkM)
    }
  }

  if (treatment_effect == "fixed") {
    trt_fix_vec <- ifelse(!is_check, trt, NA)
    trt_inc     <- .make_sparse_incidence(trt_fix_vec)
    if (ncol(trt_inc$M) < 2)
      stop("Not enough fixed non-check treatments to compute efficiency.")
    colnames(trt_inc$M) <- paste0("Line_", colnames(trt_inc$M))
    X <- cbind(X, trt_inc$M)
    if (isTRUE(check_as_fixed))
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  # -- Random effects matrices Z -----------------------------------------------
  Zrep <- .make_sparse_incidence(repfac)$M
  Zb   <- .make_sparse_incidence(blk)$M
  Zr   <- .make_sparse_incidence(as.character(rw))$M
  Zc   <- .make_sparse_incidence(as.character(cl))$M

  colnames(Zrep) <- paste0("Rep_",   colnames(Zrep))
  colnames(Zb)   <- paste0("Block_", colnames(Zb))
  colnames(Zr)   <- paste0("Row_",   colnames(Zr))
  colnames(Zc)   <- paste0("Col_",   colnames(Zc))

  Z_list  <- list(Zrep = Zrep, Zb = Zb, Zr = Zr, Zc = Zc)
  Zg      <- NULL
  trt_idx <- integer(0)

  if (treatment_effect == "random") {
    g_vec <- ifelse(!is_check, trt, NA)
    g_inc <- .make_sparse_incidence(g_vec)
    Zg    <- g_inc$M
    if (ncol(Zg) < 2)
      stop("Not enough random non-check treatments to compute efficiency.")
    colnames(Zg) <- paste0("Line_", colnames(Zg))
    Z_list$Zg    <- Zg
  }

  Z <- do.call(cbind, Z_list)

  # -- Mixed model equations ---------------------------------------------------
  QZ    <- Q %*% Z
  XtQX  <- Matrix::Matrix(Matrix::crossprod(X, Q %*% X), sparse = TRUE)
  XtQZ  <- Matrix::Matrix(Matrix::crossprod(X, QZ),      sparse = TRUE)
  ZtQX  <- Matrix::t(XtQZ)
  ZtQZ  <- Matrix::Matrix(Matrix::crossprod(Z, QZ),      sparse = TRUE)

  # -- G inverse --------------------------------------------------------------
  Ginv <- Matrix::Diagonal(ncol(Z), 0)
  idx0 <- 0L

  pRep <- ncol(Zrep); pB <- ncol(Zb); pR <- ncol(Zr); pC <- ncol(Zc)

  if (pRep > 0) {
    ii <- (idx0 + 1L):(idx0 + pRep)
    Ginv[ii, ii] <- Matrix::Diagonal(pRep, 1 / sigma_rep2)
    idx0 <- idx0 + pRep
  }
  if (pB > 0) {
    ii <- (idx0 + 1L):(idx0 + pB)
    Ginv[ii, ii] <- Matrix::Diagonal(pB, 1 / sigma_ib2)
    idx0 <- idx0 + pB
  }
  if (pR > 0) {
    ii <- (idx0 + 1L):(idx0 + pR)
    Ginv[ii, ii] <- Matrix::Diagonal(pR, 1 / sigma_r2)
    idx0 <- idx0 + pR
  }
  if (pC > 0) {
    ii <- (idx0 + 1L):(idx0 + pC)
    Ginv[ii, ii] <- Matrix::Diagonal(pC, 1 / sigma_c2)
    idx0 <- idx0 + pC
  }

  if (treatment_effect == "random") {
    pG      <- ncol(Zg)
    trt_idx <- (idx0 + 1L):(idx0 + pG)

    if (prediction_type == "IID") {
      Ginv[trt_idx, trt_idx] <- Matrix::Diagonal(pG, 1 / sigma_g2)

    } else if (prediction_type %in% c("GBLUP", "PBLUP")) {
      if (is.null(rownames(K)) || is.null(colnames(K)))
        stop("K must have rownames and colnames.")

      line_trt_order <- sub("^Line_", "", colnames(Zg))
      if (is.null(line_id_map)) {
        line_ids <- setNames(line_trt_order, line_trt_order)
      } else {
        if (!is.data.frame(line_id_map) ||
            !all(c("Treatment", "LineID") %in% names(line_id_map)))
          stop("line_id_map must be a data.frame with columns: Treatment, LineID")
        line_ids <- setNames(line_id_map$LineID, line_id_map$Treatment)
      }

      ids  <- unname(line_ids[line_trt_order])
      miss <- setdiff(ids, rownames(K))
      if (length(miss) > 0) stop("Some line IDs are missing in K rownames/colnames.")

      Ksub2    <- K[ids, ids, drop = FALSE]
      Kinv_try <- try({
        Km   <- Matrix::Matrix(Ksub2, sparse = FALSE)
        facK <- Matrix::Cholesky(Km, LDL = FALSE, Imult = 0)
        Matrix::solve(facK, Matrix::Diagonal(nrow(Km), 1))
      }, silent = TRUE)

      Kinv <- if (!inherits(Kinv_try, "try-error")) Kinv_try else
        Matrix::Matrix(.pinv_sym_dense(Ksub2), sparse = FALSE)

      Ginv[trt_idx, trt_idx] <- (1 / sigma_g2) * Kinv
    }
  }

  # -- Assemble and optionally densify C --------------------------------------
  Cmat <- rbind(cbind(XtQX, XtQZ), cbind(ZtQX, ZtQZ + Ginv))
  Cmat <- if (spatial_engine_use == "dense") as.matrix(Cmat) else
    Matrix::Matrix(Cmat, sparse = TRUE)

  # -- Compute criterion ------------------------------------------------------
  eff <- list(
    model                      = "mu + Check + Entry + Rep + IBlock(Rep) + Row + Column + e",
    treatment_effect           = treatment_effect,
    prediction_type            = if (treatment_effect == "random") prediction_type else NA_character_,
    residual_structure_requested = residual_structure,
    residual_structure_used    = residual_use,
    spatial_engine_used        = spatial_engine_use
  )

  # -- FIXED treatment branch -------------------------------------------------
  if (treatment_effect == "fixed") {
    xnames   <- colnames(X)
    trt_cols <- grep("^(Check_|Line_)", xnames)
    if (length(trt_cols) < 2)
      stop("Not enough fixed treatment columns to compute efficiency.")
    p <- length(trt_cols)
    q <- p - 1L

    if (p <= eff_full_max) {
      B    <- Matrix::sparseMatrix(i = trt_cols, j = seq_along(trt_cols),
                                   x = 1, dims = c(nrow(Cmat), p))
      Xsol <- .solve_C(Cmat, B)
      Vsub <- as.matrix(Xsol[trt_cols, , drop = FALSE])

      mean_var_diff <- .pairwise_diff_mean_var(Vsub)
      H             <- diag(p) - matrix(1 / p, p, p)
      Vctr          <- H %*% Vsub %*% H
      logdet        <- .safe_logdet_psd_dense(Vctr)

      A_criterion  <- mean_var_diff
      D_criterion  <- if (is.finite(logdet)) exp(logdet / q) else NA_real_
      A_efficiency <- 1 / A_criterion
      D_efficiency <- if (!is.na(D_criterion) && D_criterion > 0) 1 / D_criterion else NA_real_

      eff$mode         <- "FIXED_TREATMENT_BLUE_CONTRAST"
      eff$A_criterion  <- A_criterion
      eff$D_criterion  <- D_criterion
      eff$A_efficiency <- A_efficiency
      eff$D_efficiency <- D_efficiency
      eff$A            <- A_efficiency          # backward-compatible alias
      eff$D            <- D_efficiency          # backward-compatible alias
      eff$mean_VarDiff <- mean_var_diff
      eff$n_trt        <- p
      eff$n_contrasts  <- q

    } else {
      tr_est       <- .trace_subinv_est(Cmat, trt_cols, m = eff_trace_samples, seed_local = 1)
      mean_var_coef <- tr_est / p

      A_criterion  <- mean_var_coef
      A_efficiency <- 1 / A_criterion

      eff$mode         <- "FIXED_TREATMENT_BLUE_APPROX"
      eff$A_criterion  <- A_criterion
      eff$D_criterion  <- NA_real_
      eff$A_efficiency <- A_efficiency
      eff$D_efficiency <- NA_real_
      eff$A            <- A_efficiency
      eff$D            <- NA_real_
      eff$mean_VarCoef <- mean_var_coef
      eff$mean_Var     <- mean_var_coef
      eff$n_trt        <- p
      eff$n_contrasts  <- q
    }

  # -- RANDOM treatment branch ------------------------------------------------
  } else {
    p <- length(trt_idx)
    if (p < 2) stop("Not enough random non-check treatments to compute efficiency.")

    if (p <= eff_full_max) {
      B      <- Matrix::sparseMatrix(i = trt_idx, j = seq_along(trt_idx),
                                     x = 1, dims = c(nrow(Cmat), p))
      Xsol   <- .solve_C(Cmat, B)
      PEVsub <- as.matrix(Xsol[trt_idx, , drop = FALSE])
      mean_pev <- mean(diag(PEVsub))

      # -- CDmean (Rincent et al. 2012) ----------------------------------------
      # CDmean = mean coefficient of determination for GEBV prediction.
      # For each line i: CD_i = 1 - PEV_i / sigma_g2
      # CDmean = mean(CD_i) = 1 - mean_PEV / sigma_g2
      #
      # Interpretation: proportion of genetic variance explained by prediction.
      # Range [0, 1]; higher is better.
      #
      # When prediction_type = "GBLUP"/"PBLUP", sigma_g2 is the genomic
      # variance scalar used in G^{-1}. For IID, it is the iid variance.
      # CDmean is only meaningful for random treatment effects.
      CDmean <- 1 - mean_pev / sigma_g2

      # Per-line CD vector (diagonal of I - PEV/sigma_g2)
      CD_per_line <- 1 - diag(PEVsub) / sigma_g2

      eff$mode          <- "RANDOM_TREATMENT_PEV"
      eff$PEV_criterion <- mean_pev
      eff$A_criterion   <- mean_pev
      eff$D_criterion   <- NA_real_
      eff$A_efficiency  <- 1 / mean_pev
      eff$D_efficiency  <- NA_real_
      eff$A             <- eff$A_efficiency
      eff$D             <- NA_real_
      eff$mean_PEV      <- mean_pev
      eff$CDmean        <- CDmean
      eff$CD_per_line   <- CD_per_line
      eff$n_trt         <- p

    } else {
      tr_est   <- .trace_subinv_est(Cmat, trt_idx, m = eff_trace_samples, seed_local = 1)
      mean_pev <- tr_est / p

      # CDmean approximation from stochastic trace estimate
      CDmean <- 1 - mean_pev / sigma_g2

      eff$mode          <- "RANDOM_TREATMENT_PEV_APPROX"
      eff$PEV_criterion <- mean_pev
      eff$A_criterion   <- mean_pev
      eff$D_criterion   <- NA_real_
      eff$A_efficiency  <- 1 / mean_pev
      eff$D_efficiency  <- NA_real_
      eff$A             <- eff$A_efficiency
      eff$D             <- NA_real_
      eff$mean_PEV      <- mean_pev
      eff$CDmean        <- CDmean
      eff$CD_per_line   <- NA_real_   # not available in approximate mode
      eff$n_trt         <- p
    }
  }

  eff
}
