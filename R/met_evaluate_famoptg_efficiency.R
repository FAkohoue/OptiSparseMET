# ==============================================================================
# met_evaluate_famoptg_efficiency.R
#
# OptiSparseMET version of evaluate_famoptg_efficiency() from OptiDesign.
# The met_ prefix prevents namespace conflicts when both packages are loaded.
# Function body is identical to evaluate_famoptg_efficiency(). Only the
# name changes, and varcomp uses sigma_b2 (no Rep/IBlock nesting).
# ==============================================================================

#' Evaluate the statistical efficiency of a repeated-check block design
#'
#' @description
#' `met_evaluate_famoptg_efficiency()` is the OptiSparseMET version of
#' `evaluate_famoptg_efficiency()` from the OptiDesign package. The `met_`
#' prefix avoids namespace conflicts when both packages are loaded
#' simultaneously. All arguments, return values, and internal logic are
#' identical to `evaluate_famoptg_efficiency()` in OptiDesign.
#'
#' `met_evaluate_famoptg_efficiency()` takes a `field_book` produced by
#' [met_prep_famoptg()] and computes optimality criteria under a
#' user-specified mixed model. It is fully decoupled from design construction
#' so that a single design can be evaluated multiple times under different
#' model assumptions without rebuilding the layout.
#'
#' This function is the sibling of `met_evaluate_alpha_efficiency()` for
#' [met_alpha_rc_stream()] designs. The key structural difference is that
#' [met_prep_famoptg()] designs have a single flat blocking structure (one set
#' of `n_blocks` blocks with no replicate or incomplete-block nesting), so
#' the random effect model here contains **Block + Row + Column** rather than
#' the **Rep + IBlock(Rep) + Row + Column** structure of alpha-lattice designs.
#'
#' @details
#' ## Mixed model
#'
#' \deqn{y = X\beta + Zu + e}
#'
#' **Fixed part** \eqn{X}: intercept, optional check fixed effects
#' (`check_as_fixed = TRUE`), entry fixed effects when
#' `treatment_effect = "fixed"`.
#'
#' **Random part** \eqn{Z}: block, row, column, and - when
#' `treatment_effect = "random"` - entry treatment effects. There is no
#' replicate or incomplete-block term because [met_prep_famoptg()] has no
#' such nesting.
#'
#' **Residual structure** \eqn{R = \sigma_e^2 \Sigma}:
#' - `"IID"`: \eqn{\Sigma = I}
#' - `"AR1"`: row-only AR1,
#'   \eqn{Q_\text{AR1}(\rho_\text{row}) \otimes I_\text{cols}}
#' - `"AR1xAR1"`: separable row x column AR1,
#'   \eqn{Q_\text{AR1}(\rho_\text{col}) \otimes Q_\text{AR1}(\rho_\text{row})}
#'
#' ## Mixed model coefficient matrix
#'
#' \deqn{C = \begin{pmatrix} X^\top Q X & X^\top Q Z \\ Z^\top Q X &
#'   Z^\top Q Z + G^{-1} \end{pmatrix}}
#'
#' where \eqn{Q = R^{-1}} and \eqn{G^{-1}} is block-diagonal with:
#' - Block random effects: \eqn{\sigma_b^{-2} I}
#' - Row random effects: \eqn{\sigma_r^{-2} I}
#' - Column random effects: \eqn{\sigma_c^{-2} I}
#' - Entry random effects (when `treatment_effect = "random"`):
#'   \eqn{\sigma_g^{-2} I} (IID) or \eqn{\sigma_g^{-2} K^{-1}} (GBLUP/PBLUP)
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
#' \deqn{D_\text{criterion} =
#'   \exp\!\left(\frac{\log\det(HVH)}{p-1}\right)}
#' where \eqn{H = I_p - p^{-1}J_p} and \eqn{V} is the treatment
#' variance-covariance submatrix of \eqn{C^{-1}}.
#'
#' **Efficiency forms** (higher = better):
#' \deqn{A_\text{efficiency} = 1 / A_\text{criterion}, \quad
#'       D_\text{efficiency} = 1 / D_\text{criterion}}
#'
#' ### Random treatment effects (genomic prediction)
#'
#' **PEV-criterion** (lower = better): mean prediction error variance.
#' \deqn{\text{PEV}_\text{criterion} = \frac{1}{p}
#'   \sum_{i=1}^{p} \text{Var}(\hat{u}_i - u_i)}
#'
#' **CDmean** (higher = better): mean coefficient of determination for
#' GEBV prediction (Rincent et al. 2012, *Genetics* 192:715-728):
#' \deqn{\text{CDmean} = 1 - \frac{\text{PEV}_\text{criterion}}{\sigma_g^2}}
#'
#' **CD per line**: per-line coefficient of determination:
#' \deqn{CD_i = 1 - \text{PEV}_i / \sigma_g^2}
#' Available in full mode only; `NA` in `_APPROX` mode.
#'
#' ## Large-design approximation
#'
#' When the number of treatments exceeds `eff_full_max`, a Hutchinson
#' stochastic trace estimator with `eff_trace_samples` Rademacher vectors
#' replaces exact submatrix extraction. The `mode` field carries the suffix
#' `_APPROX` and `D_criterion`, `D_efficiency`, and `CD_per_line` are `NA`.
#'
#' @param field_book Data frame produced by [met_prep_famoptg()]. Must contain
#'   columns: `Treatment`, `Family`, `Block`, `Plot`, `Row`, `Column`.
#'
#' @param n_rows Positive integer. Number of field rows (must match the
#'   original design).
#'
#' @param n_cols Positive integer. Number of field columns (must match the
#'   original design).
#'
#' @param check_treatments Character vector of check treatment identifiers
#'   (must match those used in [met_prep_famoptg()]).
#'
#' @param treatment_effect Character. Whether entry treatments are modelled
#'   as `"fixed"` (BLUE-based A and D criteria) or `"random"` (PEV-based
#'   criterion and CDmean). Note: `"fixed"` does not produce CDmean;
#'   `"random"` does not produce D-criterion.
#'
#' @param prediction_type Character. Random-effect prediction model:
#' \describe{
#'   \item{`"IID"`}{\eqn{G^{-1}_\text{entry} = \sigma_g^{-2} I}. No
#'     relationship matrix required.}
#'   \item{`"GBLUP"` / `"PBLUP"`}{\eqn{G^{-1}_\text{entry} =
#'     \sigma_g^{-2} K^{-1}}. Requires `K`.}
#'   \item{`"none"`}{Invalid when `treatment_effect = "random"`.}
#' }
#'
#' @param check_as_fixed Logical. If `TRUE`, checks are included as fixed
#'   effects in the model. Default `TRUE`.
#'
#' @param residual_structure Character. Residual covariance structure:
#'   `"IID"`, `"AR1"` (row AR1 only), or `"AR1xAR1"` (separable row x
#'   column AR1).
#'
#' @param rho_row Numeric in \eqn{(-1, 1)}. AR1 autocorrelation along rows.
#'   Used when `residual_structure %in% c("AR1", "AR1xAR1")`.
#'
#' @param rho_col Numeric in \eqn{(-1, 1)}. AR1 autocorrelation along
#'   columns. Used only when `residual_structure = "AR1xAR1"`.
#'
#' @param varcomp Named list of variance components. Required components:
#' \describe{
#'   \item{`sigma_e2`}{Residual variance.}
#'   \item{`sigma_g2`}{Entry genetic variance. Used as denominator of
#'     CDmean when `treatment_effect = "random"`.}
#'   \item{`sigma_b2`}{Block variance. Corresponds to the single flat
#'     blocking level in [met_prep_famoptg()] designs (no replicate or
#'     incomplete-block nesting). This differs from
#'     `met_evaluate_alpha_efficiency()` which uses `sigma_rep2` and
#'     `sigma_ib2`.}
#'   \item{`sigma_r2`}{Row variance.}
#'   \item{`sigma_c2`}{Column variance.}
#' }
#'
#' @param K Optional square numeric matrix with rownames and colnames.
#'   Required when `prediction_type %in% c("GBLUP", "PBLUP")`.
#'
#' @param line_id_map Optional data frame with columns `Treatment` and
#'   `LineID`. Required when treatment labels do not match `rownames(K)`.
#'
#' @param spatial_engine Character. Computational engine: `"auto"` uses
#'   `"dense"` when used plots \eqn{\leq} `dense_max_n`, otherwise
#'   `"sparse"`.
#'
#' @param dense_max_n Integer. Threshold for `spatial_engine = "auto"`.
#'
#' @param eff_trace_samples Positive integer. Number of Rademacher vectors
#'   for the Hutchinson stochastic trace estimator (large designs only).
#'
#' @param eff_full_max Positive integer. Maximum number of treatments for
#'   exact \eqn{C^{-1}} submatrix extraction. Designs with more treatments
#'   use the stochastic estimator (`mode` suffix `_APPROX`).
#'
#' @return A named list containing model metadata and criterion values:
#' \describe{
#'   \item{`model`}{Character. Model string.}
#'   \item{`treatment_effect`}{Character. As supplied.}
#'   \item{`prediction_type`}{Character or `NA`.}
#'   \item{`residual_structure_requested`}{Character. As supplied.}
#'   \item{`residual_structure_used`}{Character. As resolved.}
#'   \item{`spatial_engine_used`}{Character. `"dense"` or `"sparse"`.}
#'   \item{`mode`}{Character. One of `"FIXED_TREATMENT_BLUE_CONTRAST"`,
#'     `"FIXED_TREATMENT_BLUE_APPROX"`, `"RANDOM_TREATMENT_PEV"`,
#'     `"RANDOM_TREATMENT_PEV_APPROX"`.}
#'   \item{`n_trt`}{Integer. Number of treatment columns evaluated.}
#'   \item{`n_contrasts`}{Integer. `n_trt - 1`. Fixed mode only.}
#'   \item{`A_criterion`}{Numeric. Mean pairwise contrast variance (fixed)
#'     or mean PEV (random). Lower is better.}
#'   \item{`D_criterion`}{Numeric or `NA`. Fixed full mode only.
#'     Lower is better.}
#'   \item{`A_efficiency`}{Numeric. `1 / A_criterion`. Higher is better.}
#'   \item{`D_efficiency`}{Numeric or `NA`. `1 / D_criterion`. Higher is
#'     better.}
#'   \item{`A`}{Numeric. Alias for `A_efficiency` (backward compatibility).}
#'   \item{`D`}{Numeric or `NA`. Alias for `D_efficiency` (backward
#'     compatibility).}
#'   \item{`mean_VarDiff`}{Numeric. Mean pairwise contrast variance.
#'     Fixed full mode only.}
#'   \item{`PEV_criterion`}{Numeric. Mean PEV. Random mode only.}
#'   \item{`mean_PEV`}{Numeric. Alias for `PEV_criterion`. Random mode only.}
#'   \item{`CDmean`}{Numeric. Mean coefficient of determination:
#'     \eqn{1 - \text{mean\_PEV} / \sigma_g^2}. Range \eqn{[0, 1]}.
#'     Higher is better. Random mode only.}
#'   \item{`CD_per_line`}{Numeric vector or `NA`. Per-line coefficient of
#'     determination. Available in full mode only. Random mode only.}
#' }
#'
#' @references
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P.,
#' ... & Moreau, L. (2012). Maximizing the reliability of genomic selection by
#' optimizing the calibration set of reference individuals.
#' *Genetics*, 192(2), 715-728.
#'
#' @seealso [met_prep_famoptg()] to construct the design.
#'   [met_optimize_famoptg()] to search for a criterion-optimal design.
#'   [met_evaluate_alpha_efficiency()] for the equivalent function for
#'   [met_alpha_rc_stream()] designs (different random effect structure:
#'   uses `sigma_rep2` and `sigma_ib2` instead of `sigma_b2`).
#'
#' @examples
#' design <- met_prep_famoptg(
#'   check_treatments        = c("CHK1", "CHK2"),
#'   check_families          = c("CHECK", "CHECK"),
#'   p_rep_treatments        = paste0("P", 1:20),
#'   p_rep_reps              = rep(2L, 20),
#'   p_rep_families          = rep(paste0("F", 1:4), 5),
#'   unreplicated_treatments = paste0("U", 1:60),
#'   unreplicated_families   = rep(paste0("F", 1:4), 15),
#'   n_blocks = 5, n_rows = 15, n_cols = 20
#' )
#'
#' ## Fixed-effect evaluation (BLUE contrasts, IID residuals)
#' eff_fixed <- met_evaluate_famoptg_efficiency(
#'   field_book       = design$field_book,
#'   n_rows = 15, n_cols = 20,
#'   check_treatments = c("CHK1", "CHK2"),
#'   treatment_effect = "fixed"
#' )
#' eff_fixed$A_criterion   # lower is better
#' eff_fixed$D_criterion   # lower is better
#' eff_fixed$A_efficiency  # higher is better
#'
#' ## AR1xAR1 spatial model
#' eff_spatial <- met_evaluate_famoptg_efficiency(
#'   field_book         = design$field_book,
#'   n_rows = 15, n_cols = 20,
#'   check_treatments   = c("CHK1", "CHK2"),
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row = 0.10, rho_col = 0.10
#' )
#' eff_spatial$A_criterion
#'
#' \dontrun{
#' ## Random-effect evaluation with GBLUP and CDmean
#' eff_gblup <- met_evaluate_famoptg_efficiency(
#'   field_book       = design$field_book,
#'   n_rows = 15, n_cols = 20,
#'   check_treatments = c("CHK1", "CHK2"),
#'   treatment_effect = "random",
#'   prediction_type  = "GBLUP",
#'   K                = my_kinship_matrix,
#'   varcomp          = list(
#'     sigma_g2 = 0.4, sigma_e2 = 0.6,
#'     sigma_b2 = 0.1, sigma_r2 = 0.02, sigma_c2 = 0.02
#'   )
#' )
#' eff_gblup$CDmean       # mean GEBV prediction reliability [0, 1]
#' eff_gblup$CD_per_line  # per-line reliability vector
#' eff_gblup$mean_PEV     # mean prediction error variance
#' }
#'
#' @importFrom Matrix Diagonal sparseMatrix crossprod solve t Cholesky
#' @export
met_evaluate_famoptg_efficiency <- function(
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
    # Note: sigma_b2 = block variance (replaces sigma_rep2 + sigma_ib2
    # from evaluate_alpha_efficiency; prep_famoptg has no rep/IBlock structure)
    varcomp = list(
      sigma_e2 = 1,
      sigma_g2 = 1,
      sigma_b2 = 1,
      sigma_r2 = 1,
      sigma_c2 = 1
    ),

    # -- Genomic / pedigree relationship matrix -------------------------------
    K           = NULL,
    line_id_map = NULL,

    # -- Computation controls -------------------------------------------------
    spatial_engine    = c("auto", "sparse", "dense"),
    dense_max_n       = 5000,
    eff_trace_samples = 80,
    eff_full_max      = 400
) {

  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required.")

  # -- Argument matching -------------------------------------------------------
  treatment_effect   <- match.arg(treatment_effect)
  prediction_type    <- match.arg(prediction_type)
  residual_structure <- match.arg(residual_structure)
  spatial_engine     <- match.arg(spatial_engine)

  # -- Pre-flight guards -------------------------------------------------------
  if (treatment_effect == "random" && prediction_type == "none") {
    stop(paste0(
      "Cannot compute efficiency: treatment_effect = 'random' and ",
      "prediction_type = 'none'.\n",
      "Set treatment_effect = 'fixed', or choose prediction_type in ",
      "c('IID', 'GBLUP', 'PBLUP')."
    ))
  }
  if (prediction_type %in% c("GBLUP", "PBLUP") && is.null(K))
    stop("K must be provided when prediction_type is 'GBLUP' or 'PBLUP'.")

  if (!is.list(varcomp) ||
      !all(c("sigma_e2", "sigma_g2", "sigma_b2", "sigma_r2", "sigma_c2") %in% names(varcomp)))
    stop("varcomp must contain: sigma_e2, sigma_g2, sigma_b2, sigma_r2, sigma_c2.")

  if (!is.numeric(rho_row) || length(rho_row) != 1 || abs(rho_row) >= 1)
    stop("rho_row must satisfy |rho_row| < 1.")
  if (!is.numeric(rho_col) || length(rho_col) != 1 || abs(rho_col) >= 1)
    stop("rho_col must satisfy |rho_col| < 1.")

  required_cols <- c("Treatment", "Family", "Block", "Plot", "Row", "Column")
  missing_cols  <- setdiff(required_cols, names(field_book))
  if (length(missing_cols) > 0)
    stop(paste0("field_book is missing columns: ", paste(missing_cols, collapse = ", ")))

  # -- Unpack varcomp ----------------------------------------------------------
  sigma_e2 <- varcomp$sigma_e2
  sigma_g2 <- varcomp$sigma_g2
  sigma_b2 <- varcomp$sigma_b2
  sigma_r2 <- varcomp$sigma_r2
  sigma_c2 <- varcomp$sigma_c2

  # -- Working data ------------------------------------------------------------
  fb       <- field_book
  nn       <- nrow(fb)
  if (nn == 0) stop("field_book is empty.")

  trt      <- as.character(fb$Treatment)
  is_check <- trt %in% check_treatments
  blk      <- paste0("Block_", fb$Block)
  rw       <- as.character(fb$Row)
  cl       <- as.character(fb$Column)

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
    grid_index <- (as.integer(fb$Row) - 1L) * n_cols + as.integer(fb$Column)
    Qfull      <- Matrix::kronecker(Matrix::Diagonal(n_cols, 1), Qrow) * (1 / sigma_e2)
    Q          <- Qfull[grid_index, grid_index, drop = FALSE]

  } else {
    # AR1xAR1
    Qrow       <- .ar1_precision_sparse(n_rows, rho_row)
    Qcol       <- .ar1_precision_sparse(n_cols, rho_col)
    Qgrid      <- Matrix::kronecker(Qcol, Qrow) * (1 / sigma_e2)
    grid_index <- (as.integer(fb$Row) - 1L) * n_cols + as.integer(fb$Column)
    Q          <- Qgrid[grid_index, grid_index, drop = FALSE]
  }

  # -- Fixed effects matrix X --------------------------------------------------
  X_int           <- Matrix::Matrix(rep(1, nn), ncol = 1, sparse = TRUE)
  colnames(X_int) <- "(Intercept)"
  X               <- X_int

  if (isTRUE(check_as_fixed)) {
    chk_vec <- ifelse(is_check, trt, NA_character_)
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
    trt_fix_vec <- ifelse(!is_check, trt, NA_character_)
    trt_inc     <- .make_sparse_incidence(trt_fix_vec)
    if (ncol(trt_inc$M) < 2)
      stop("Not enough fixed non-check treatments to compute efficiency.")
    colnames(trt_inc$M) <- paste0("Line_", colnames(trt_inc$M))
    X <- cbind(X, trt_inc$M)
    # Remove redundant intercept when full treatment dummy set is present
    if (isTRUE(check_as_fixed))
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  # -- Random effects matrices Z -----------------------------------------------
  # prep_famoptg model: Block + Row + Column (no Rep, no IBlock)
  Zb <- .make_sparse_incidence(blk)$M
  Zr <- .make_sparse_incidence(rw)$M
  Zc <- .make_sparse_incidence(cl)$M

  colnames(Zb) <- paste0("Block_", colnames(Zb))
  colnames(Zr) <- paste0("Row_",   colnames(Zr))
  colnames(Zc) <- paste0("Col_",   colnames(Zc))

  Z_list  <- list(Zb = Zb, Zr = Zr, Zc = Zc)
  Zg      <- NULL
  trt_idx <- integer(0)

  if (treatment_effect == "random") {
    g_vec <- ifelse(!is_check, trt, NA_character_)
    g_inc <- .make_sparse_incidence(g_vec)
    Zg    <- g_inc$M
    if (ncol(Zg) < 2)
      stop("Not enough random non-check treatments to compute efficiency.")
    colnames(Zg) <- paste0("Line_", colnames(Zg))
    Z_list$Zg    <- Zg
  }

  Z <- do.call(cbind, Z_list)

  # -- Mixed model equations ---------------------------------------------------
  QZ   <- Q %*% Z
  XtQX <- Matrix::Matrix(Matrix::crossprod(X, Q %*% X), sparse = TRUE)
  XtQZ <- Matrix::Matrix(Matrix::crossprod(X, QZ),      sparse = TRUE)
  ZtQX <- Matrix::t(XtQZ)
  ZtQZ <- Matrix::Matrix(Matrix::crossprod(Z, QZ),      sparse = TRUE)

  # -- G inverse --------------------------------------------------------------
  Ginv <- Matrix::Diagonal(ncol(Z), 0)
  idx0 <- 0L

  pB <- ncol(Zb); pR <- ncol(Zr); pC <- ncol(Zc)

  if (pB > 0) {
    ii <- (idx0 + 1L):(idx0 + pB)
    Ginv[ii, ii] <- Matrix::Diagonal(pB, 1 / sigma_b2)
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
          stop("line_id_map must be a data.frame with columns: Treatment, LineID.")
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

  # -- Assemble C matrix -------------------------------------------------------
  Cmat <- rbind(cbind(XtQX, XtQZ), cbind(ZtQX, ZtQZ + Ginv))
  Cmat <- if (spatial_engine_use == "dense") as.matrix(Cmat) else
    Matrix::Matrix(Cmat, sparse = TRUE)

  # -- Result list skeleton ----------------------------------------------------
  eff <- list(
    model                        = "mu + Check + Entry + Block + Row + Column + e",
    treatment_effect             = treatment_effect,
    prediction_type              = if (treatment_effect == "random") prediction_type else NA_character_,
    residual_structure_requested = residual_structure,
    residual_structure_used      = residual_use,
    spatial_engine_used          = spatial_engine_use
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
      eff$A            <- A_efficiency       # backward-compatible alias
      eff$D            <- D_efficiency       # backward-compatible alias
      eff$mean_VarDiff <- mean_var_diff
      eff$n_trt        <- p
      eff$n_contrasts  <- q

    } else {
      tr_est        <- .trace_subinv_est(Cmat, trt_cols, m = eff_trace_samples, seed_local = 1)
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

      # CDmean: mean coefficient of determination for GEBV prediction
      # CDmean = 1 - mean_PEV / sigma_g2
      # Range [0,1]; higher = more reliable genomic predictions
      CDmean      <- 1 - mean_pev / sigma_g2
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
