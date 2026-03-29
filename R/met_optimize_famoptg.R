# ==============================================================================
# met_optimize_famoptg.R
#
# OptiSparseMET version of optimize_famoptg() from OptiDesign.
# The met_ prefix prevents namespace conflicts when both packages are loaded.
# Function body is identical to optimize_famoptg(). Only the name and all
# internal met_ call references change.
# ==============================================================================

#' Search for a criterion-optimal repeated-check block design
#'
#' @description
#' `met_optimize_famoptg()` is the OptiSparseMET version of
#' `optimize_famoptg()` from the OptiDesign package. The `met_` prefix avoids
#' namespace conflicts when both packages are loaded simultaneously. All
#' arguments, return values, and internal logic are identical to
#' `optimize_famoptg()` in OptiDesign.
#'
#' `met_optimize_famoptg()` wraps [met_prep_famoptg()] and
#' [met_evaluate_famoptg_efficiency()] in a Random Restart (RS) optimisation
#' loop that searches for the design with the best optimality criterion
#' value among `n_restarts` independent randomisations.
#'
#' **Why Random Restart only?** [met_prep_famoptg()] enforces the p-rep
#' constraint - that no replicated treatment appears twice in the same block
#' - by construction at every call. Permutation-based methods such as
#' Simulated Annealing or a Genetic Algorithm would require block-aware swap
#' logic to preserve this constraint, making them substantially more complex
#' for modest criterion improvement. RS generates fully valid designs at
#' every restart with no risk of constraint violation.
#'
#' **Constraint preservation**: every candidate is produced by
#' [met_prep_famoptg()], so the following structural guarantees hold by
#' construction for all candidates, regardless of `n_restarts`:
#' \itemize{
#'   \item Check treatments appear in every block.
#'   \item P-rep treatments appear in exactly `p_rep_reps[i]` blocks each,
#'     always in distinct blocks - never twice in the same block.
#'   \item Unreplicated treatments appear exactly once.
#'   \item Family/cluster adjacency is minimised within blocks.
#'   \item Optional genetic dispersion is applied.
#' }
#'
#' **Integrity guarantee**: every candidate additionally passes
#' `.check_famoptg_integrity()` before being scored or stored as the best.
#' The running best is updated only when both the criterion score improves
#' and an integrity re-check passes. A final integrity check is performed
#' before returning. If no valid design is found an emergency fallback
#' returns a single freshly constructed valid design.
#'
#' @details
#' ## Optimisation target
#'
#' All restarts minimise an internal score (lower = better internally):
#'
#' | `criterion` | Internal score | Direction reported to user |
#' |---|---|---|
#' | `"A"` | `A_criterion` | Lower is better |
#' | `"D"` | `D_criterion` | Lower is better |
#' | `"both"` | Mean of `A_criterion` and `D_criterion` | Lower is better |
#' | `"CDmean"` | Negated CDmean | Higher CDmean is better |
#'
#' For `"CDmean"` the internal negation is transparent to the user:
#' `best_score` and `score_history` are always reported as positive CDmean
#' values (higher = better).
#'
#' ## Integrity checking
#'
#' Five structural constraints are verified for every candidate:
#' \enumerate{
#'   \item No non-check entry appears more than once in any single block.
#'   \item Each p-rep treatment appears in exactly `p_rep_reps[i]` blocks.
#'   \item Each unreplicated treatment appears exactly once.
#'   \item All check treatments appear in every block.
#'   \item No p-rep treatment occupies the same block twice (core p-rep
#'     constraint).
#' }
#'
#' Checks 2 and 5 together constitute the p-rep guarantee. Failure of any
#' check discards the candidate and counts toward `n_failed`.
#'
#' ## Failure handling
#'
#' Failed restarts (construction error or integrity failure) are counted.
#' If the failure rate exceeds `max_failure_rate`, the function stops with a
#' diagnostic message. Below the threshold a warning is issued. If no valid
#' design is found after all restarts, up to 10 emergency fallback attempts
#' are made. If all fail, the function stops with an informative error.
#'
#' @param check_treatments,check_families,p_rep_treatments,p_rep_reps,p_rep_families,unreplicated_treatments,unreplicated_families,n_blocks,n_rows,n_cols,order,serpentine,seed,attempts,warn_and_correct,fix_rows,cluster_source,GRM,A,id_map,cluster_method,cluster_seed,cluster_attempts,n_pcs_use,check_placement,check_opt_attempts,use_dispersion,dispersion_source,dispersion_radius,dispersion_iters,dispersion_seed,K,line_id_map
#'   Passed directly to [met_prep_famoptg()]. See that function for full
#'   documentation.
#'
#' @param treatment_effect,prediction_type,check_as_fixed,residual_structure,rho_row,rho_col,varcomp,spatial_engine,dense_max_n,eff_trace_samples,eff_full_max
#'   Passed directly to [met_evaluate_famoptg_efficiency()]. See that function
#'   for full documentation. Note that `varcomp` here uses `sigma_b2` (block
#'   variance) instead of `sigma_rep2` and `sigma_ib2`.
#'
#' @param criterion Character. Optimality criterion to drive the search:
#' \describe{
#'   \item{`"A"`}{Minimise `A_criterion` (mean pairwise contrast variance
#'     for fixed effects; mean PEV for random effects). Valid for both
#'     `treatment_effect = "fixed"` and `"random"`. Recommended default.}
#'   \item{`"D"`}{Minimise `D_criterion` (geometric mean of contrast
#'     covariance eigenvalues). Fixed effects only. Falls back to
#'     `A_criterion` with a warning for random effects.}
#'   \item{`"both"`}{Minimise mean of `A_criterion` and `D_criterion`.
#'     Fixed effects only.}
#'   \item{`"CDmean"`}{Maximise CDmean (mean coefficient of determination
#'     for GEBV prediction; Rincent et al. 2012). Requires
#'     `treatment_effect = "random"` and
#'     `prediction_type %in% c("IID", "GBLUP", "PBLUP")`. Best suited for
#'     optimising designs for genomic selection.}
#' }
#'
#' @param n_restarts Positive integer. Number of independent random restarts.
#'   Each restart calls [met_prep_famoptg()] with a new seed and evaluates the
#'   resulting design. Default 50.
#'
#' @param max_failure_rate Numeric in \eqn{[0, 1]}. Maximum tolerated
#'   fraction of restarts that may fail construction or integrity checking
#'   before the function stops with a diagnostic error. Default `0.5`.
#'   Below the threshold a warning is issued. Increase for very constrained
#'   field geometries where some failure is expected.
#'
#' @param verbose_opt Logical. If `TRUE`, prints per-restart progress
#'   messages showing current score, running best, and failure status.
#'   Default `TRUE`.
#'
#' @return The return value of [met_prep_famoptg()] for the best valid design
#'   found, augmented with two additional components:
#' \describe{
#'   \item{`efficiency`}{The full output of [met_evaluate_famoptg_efficiency()]
#'     for the best design, including all applicable criterion values:
#'     `A_criterion`, `D_criterion`, `A_efficiency`, `D_efficiency`,
#'     `CDmean`, and `CD_per_line`.}
#'   \item{`optimization`}{Named list of optimisation metadata:
#'     \describe{
#'       \item{`method`}{Character. Always `"RS"`.}
#'       \item{`criterion`}{Character. As supplied.}
#'       \item{`best_score`}{Numeric. Best criterion value found, in natural
#'         direction: lower is better for `"A"`, `"D"`, `"both"`; higher is
#'         better for `"CDmean"` (positive value, negation is internal only).}
#'       \item{`score_history`}{Numeric vector of length `n_restarts`.
#'         Criterion value for each restart (`NA` for failed restarts). For
#'         `criterion = "CDmean"`, values are positive CDmean (higher =
#'         better). For all other criteria, lower = better.}
#'       \item{`master_seed`}{Integer. The master random seed used.}
#'       \item{`n_restarts`}{Integer. As supplied.}
#'       \item{`n_failed`}{Integer. Number of restarts that failed
#'         construction or integrity checking.}
#'     }
#'   }
#' }
#'
#' If no valid optimised design is found after all restarts, a warning is
#' issued and the function returns a single freshly constructed valid design
#' via emergency fallback. If even that fails after 10 attempts, the function
#' stops with a diagnostic error pointing to the likely cause.
#'
#' @references
#' Rincent, R., Laloe, D., Nicolas, S., Altmann, T., Brunel, D., Revilla, P.,
#' ... & Moreau, L. (2012). Maximizing the reliability of genomic selection by
#' optimizing the calibration set of reference individuals.
#' *Genetics*, 192(2), 715-728.
#'
#' Jones, B., Allen-Moyer, K., & Goos, P. (2021). A-optimal versus D-optimal
#' design of screening experiments. *Journal of Quality Technology*, 53(4),
#' 369-382.
#'
#' @seealso [met_prep_famoptg()] for single-design construction without
#'   optimisation.
#'   [met_evaluate_famoptg_efficiency()] for standalone criterion evaluation.
#'   [met_optimize_alpha_rc()] for the equivalent optimizer for
#'   [met_alpha_rc_stream()] designs (also supports SA and GA methods).
#'
#' @examples
#' ## Optimise a p-rep design by A-criterion (fixed effects)
#' result <- met_optimize_famoptg(
#'   check_treatments        = c("CHK1", "CHK2"),
#'   check_families          = c("CHECK", "CHECK"),
#'   p_rep_treatments        = paste0("P", 1:20),
#'   p_rep_reps              = rep(2L, 20),
#'   p_rep_families          = rep(paste0("F", 1:4), 5),
#'   unreplicated_treatments = paste0("U", 1:60),
#'   unreplicated_families   = rep(paste0("F", 1:4), 15),
#'   n_blocks           = 5,
#'   n_rows             = 15,
#'   n_cols             = 20,
#'   treatment_effect   = "fixed",
#'   residual_structure = "AR1xAR1",
#'   rho_row            = 0.10,
#'   rho_col            = 0.10,
#'   criterion          = "A",
#'   n_restarts         = 50
#' )
#'
#' result$efficiency$A_criterion      # best A-criterion found
#' result$optimization$best_score     # same value, lower is better
#' result$optimization$score_history  # per-restart scores
#' result$optimization$n_failed       # failed restarts
#'
#' \dontrun{
#' ## Optimise an augmented design by CDmean (GBLUP)
#' result_cdmean <- met_optimize_famoptg(
#'   check_treatments        = c("CHK1", "CHK2"),
#'   check_families          = c("CHECK", "CHECK"),
#'   p_rep_treatments        = character(0),
#'   p_rep_reps              = integer(0),
#'   p_rep_families          = character(0),
#'   unreplicated_treatments = paste0("E", 1:120),
#'   unreplicated_families   = rep(paste0("F", 1:6), 20),
#'   n_blocks           = 6,
#'   n_rows             = 13,
#'   n_cols             = 10,
#'   treatment_effect   = "random",
#'   prediction_type    = "GBLUP",
#'   K                  = my_kinship_matrix,
#'   varcomp            = list(
#'     sigma_g2 = 0.4, sigma_e2 = 0.6,
#'     sigma_b2 = 0.1, sigma_r2 = 0.02, sigma_c2 = 0.02
#'   ),
#'   criterion          = "CDmean",
#'   n_restarts         = 50
#' )
#'
#' result_cdmean$efficiency$CDmean        # mean GEBV reliability [0,1]
#' result_cdmean$optimization$best_score  # positive CDmean, higher = better
#' result_cdmean$optimization$n_failed    # failed restarts
#' }
#'
#' @importFrom stats runif
#' @export
met_optimize_famoptg <- function(

    # -- Construction arguments (forwarded to met_prep_famoptg) -------------------
    check_treatments,
    check_families,
    p_rep_treatments,
    p_rep_reps,
    p_rep_families,
    unreplicated_treatments,
    unreplicated_families,
    n_blocks,
    n_rows,
    n_cols,
    order              = "column",
    serpentine         = FALSE,
    seed               = NULL,
    attempts           = 1000,
    warn_and_correct   = TRUE,
    fix_rows           = TRUE,
    cluster_source     = c("Family", "GRM", "A"),
    GRM                = NULL,
    A                  = NULL,
    id_map             = NULL,
    cluster_method     = c("kmeans", "hclust"),
    cluster_seed       = 1,
    cluster_attempts   = 25,
    n_pcs_use          = Inf,
    check_placement    = c("random", "systematic", "optimal"),
    check_opt_attempts = 200,
    use_dispersion     = FALSE,
    dispersion_source  = c("K", "A", "GRM"),
    dispersion_radius  = 1,
    dispersion_iters   = 2000,
    dispersion_seed    = 1,
    K                  = NULL,
    line_id_map        = NULL,

    # -- Evaluation arguments (forwarded to met_evaluate_famoptg_efficiency) ------
    treatment_effect   = c("random", "fixed"),
    prediction_type    = c("IID", "GBLUP", "PBLUP", "none"),
    check_as_fixed     = TRUE,
    residual_structure = c("IID", "AR1", "AR1xAR1"),
    rho_row            = 0,
    rho_col            = 0,
    varcomp            = list(
      sigma_e2 = 1,
      sigma_g2 = 1,
      sigma_b2 = 1,
      sigma_r2 = 1,
      sigma_c2 = 1
    ),
    spatial_engine     = c("auto", "sparse", "dense"),
    dense_max_n        = 5000,
    eff_trace_samples  = 80,
    eff_full_max       = 400,

    # -- Optimizer arguments ---------------------------------------------------
    criterion           = c("A", "D", "both", "CDmean"),
    n_restarts          = 50L,
    max_failure_rate    = 0.5,   # stop if > this fraction of restarts fail

    # Verbosity
    verbose_opt = TRUE
) {

  # -- Argument matching -------------------------------------------------------
  criterion          <- match.arg(criterion)
  treatment_effect   <- match.arg(treatment_effect)
  prediction_type    <- match.arg(prediction_type)
  cluster_source     <- match.arg(cluster_source)
  cluster_method     <- match.arg(cluster_method)
  residual_structure <- match.arg(residual_structure)
  spatial_engine     <- match.arg(spatial_engine)
  check_placement    <- match.arg(check_placement)
  dispersion_source  <- match.arg(dispersion_source)

  # -- Pre-flight: criterion / model compatibility ----------------------------
  if (treatment_effect == "random" && prediction_type == "none") {
    stop(paste0(
      "Cannot optimise: treatment_effect = 'random' and prediction_type = 'none' ",
      "produces no estimable efficiency criterion.\n",
      "Set treatment_effect = 'fixed', or choose prediction_type in ",
      "c('IID', 'GBLUP', 'PBLUP')."
    ))
  }
  if (criterion %in% c("D", "both") && treatment_effect == "random") {
    warning(paste0(
      "D-criterion is not computed for random treatment effects. ",
      "criterion will fall back to 'A'."
    ))
  }
  if (criterion == "CDmean" && treatment_effect == "fixed") {
    stop(paste0(
      "criterion = 'CDmean' requires treatment_effect = 'random'. ",
      "CDmean is only defined for random treatment / genomic prediction models."
    ))
  }
  if (criterion == "CDmean" && prediction_type == "none") {
    stop("criterion = 'CDmean' requires prediction_type in c('IID', 'GBLUP', 'PBLUP').")
  }
  if (n_restarts < 1L) stop("n_restarts must be >= 1.")

  # -- Master seed -------------------------------------------------------------
  master_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else seed
  set.seed(master_seed)

  # -- Shared construction arguments ------------------------------------------
  build_args <- list(
    check_treatments        = check_treatments,
    check_families          = check_families,
    p_rep_treatments        = p_rep_treatments,
    p_rep_reps              = p_rep_reps,
    p_rep_families          = p_rep_families,
    unreplicated_treatments = unreplicated_treatments,
    unreplicated_families   = unreplicated_families,
    n_blocks                = n_blocks,
    n_rows                  = n_rows,
    n_cols                  = n_cols,
    order                   = order,
    serpentine              = serpentine,
    attempts                = attempts,
    warn_and_correct        = warn_and_correct,
    fix_rows                = fix_rows,
    cluster_source          = cluster_source,
    GRM                     = GRM,
    A                       = A,
    id_map                  = id_map,
    cluster_method          = cluster_method,
    cluster_seed            = cluster_seed,
    cluster_attempts        = cluster_attempts,
    n_pcs_use               = n_pcs_use,
    check_placement         = check_placement,
    check_opt_attempts      = check_opt_attempts,
    use_dispersion          = use_dispersion,
    dispersion_source       = dispersion_source,
    dispersion_radius       = dispersion_radius,
    dispersion_iters        = dispersion_iters,
    dispersion_seed         = dispersion_seed,
    K                       = K,
    line_id_map             = line_id_map,
    verbose                 = FALSE
  )

  # -- Shared evaluation arguments ---------------------------------------------
  eval_args <- list(
    n_rows             = n_rows,
    n_cols             = n_cols,
    check_treatments   = check_treatments,
    treatment_effect   = treatment_effect,
    prediction_type    = prediction_type,
    check_as_fixed     = check_as_fixed,
    residual_structure = residual_structure,
    rho_row            = rho_row,
    rho_col            = rho_col,
    varcomp            = varcomp,
    K                  = K,
    line_id_map        = line_id_map,
    spatial_engine     = spatial_engine,
    dense_max_n        = dense_max_n,
    eff_trace_samples  = eff_trace_samples,
    eff_full_max       = eff_full_max
  )

  # -- Criterion extraction (lower = better internally) -----------------------
  # CDmean is negated internally so the same "keep minimum" logic applies.
  # best_score and score_history are always reported in natural direction.
  extract_score <- function(eff) {
    if (is.null(eff)) return(NA_real_)
    a_val  <- if (!is.null(eff$A_criterion)) eff$A_criterion else NA_real_
    d_val  <- if (!is.null(eff$D_criterion)) eff$D_criterion else NA_real_
    cd_val <- if (!is.null(eff$CDmean))      eff$CDmean      else NA_real_

    if (criterion == "A") return(a_val)

    if (criterion == "D") {
      if (is.na(d_val)) {
        warning("D_criterion is NA; falling back to A_criterion.")
        return(a_val)
      }
      return(d_val)
    }

    if (criterion == "CDmean") {
      if (is.na(cd_val)) {
        warning("CDmean is NA; falling back to A_criterion.")
        return(a_val)
      }
      return(-cd_val)   # negate: higher CDmean -> lower internal score
    }

    # "both": mean of available lower-is-better criteria
    vals <- c(a_val, d_val); vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_real_)
    mean(vals)
  }

  # ============================================================
  # INTEGRITY CHECK
  # ============================================================
  # Validates five structural constraints:
  #   1. No entry appears more than once within any block
  #   2. P-rep treatments appear in exactly p_rep_reps[i] blocks each
  #   3. Unreplicated treatments appear exactly once
  #   4. Checks appear in every block
  #   5. P-rep treatments are not duplicated within a single block
  # ============================================================
  .check_famoptg_integrity <- function(design) {
    fb <- design$field_book
    is_check_plot <- fb$Treatment %in% check_treatments

    non_check_fb <- fb[!is_check_plot, , drop = FALSE]

    # 1. No non-check entry appears more than once within any block
    if (nrow(non_check_fb) > 0) {
      counts <- table(non_check_fb$Block, non_check_fb$Treatment)
      if (any(counts > 1)) {
        warning(paste0(
          "[met_optimize_famoptg] Integrity check failed: ",
          "non-check treatment appears more than once in a single block. ",
          "Design discarded."
        ))
        return(FALSE)
      }
    }

    # 2 & 5. P-rep treatments: correct total count and no within-block duplicates
    if (length(p_rep_treatments) > 0) {
      for (i in seq_along(p_rep_treatments)) {
        trt_i     <- p_rep_treatments[i]
        expected  <- p_rep_reps[i]
        actual    <- sum(fb$Treatment == trt_i)
        if (actual != expected) {
          warning(paste0(
            "[met_optimize_famoptg] Integrity check failed: p-rep treatment '",
            trt_i, "' appears ", actual, " times but expected ", expected, ". ",
            "Design discarded."
          ))
          return(FALSE)
        }
        # Check it appears in distinct blocks only (constraint 5)
        trt_blocks <- fb$Block[fb$Treatment == trt_i]
        if (length(trt_blocks) != length(unique(trt_blocks))) {
          warning(paste0(
            "[met_optimize_famoptg] Integrity check failed: p-rep treatment '",
            trt_i, "' appears in the same block more than once. ",
            "Design discarded."
          ))
          return(FALSE)
        }
      }
    }

    # 3. Unreplicated treatments appear exactly once
    if (length(unreplicated_treatments) > 0) {
      for (trt_u in unreplicated_treatments) {
        actual <- sum(fb$Treatment == trt_u)
        if (actual != 1L) {
          warning(paste0(
            "[met_optimize_famoptg] Integrity check failed: unreplicated treatment '",
            trt_u, "' appears ", actual, " times (expected 1). ",
            "Design discarded."
          ))
          return(FALSE)
        }
      }
    }

    # 4. All checks appear in every block
    check_fb <- fb[is_check_plot, , drop = FALSE]
    for (chk in check_treatments) {
      blocks_with_chk <- unique(check_fb$Block[check_fb$Treatment == chk])
      missing_blocks  <- setdiff(seq_len(n_blocks), blocks_with_chk)
      if (length(missing_blocks) > 0) {
        warning(paste0(
          "[met_optimize_famoptg] Integrity check failed: check treatment '",
          chk, "' is missing from block(s): ",
          paste(missing_blocks, collapse = ", "), ". Design discarded."
        ))
        return(FALSE)
      }
    }

    TRUE
  }

  # ============================================================
  # EMERGENCY FALLBACK
  # ============================================================
  .emergency_fallback <- function() {
    warning(paste0(
      "[met_optimize_famoptg] No valid optimised design found. ",
      "Attempting emergency fallback: generating a single valid design."
    ))
    for (attempt in seq_len(10L)) {
      s        <- sample.int(.Machine$integer.max, 1)
      fallback <- run_pipeline(s)
      if (!is.null(fallback) && !is.na(fallback$score)) {
        warning(paste0(
          "[met_optimize_famoptg] Emergency fallback succeeded on attempt ", attempt,
          ". Returning a non-optimised valid design."
        ))
        return(fallback)
      }
    }
    stop(paste0(
      "[met_optimize_famoptg] Emergency fallback failed after 10 attempts. ",
      "Design parameters may be infeasible. ",
      "Check n_blocks, n_rows, n_cols, p_rep_reps, and block size constraints."
    ))
  }

  # ============================================================
  # PIPELINE: construct + integrity check + evaluate
  # ============================================================
  run_pipeline <- function(s) {
    args      <- build_args
    args$seed <- s

    design <- tryCatch(do.call(met_prep_famoptg, args), error = function(e) NULL)
    if (is.null(design)) return(NULL)

    if (!.check_famoptg_integrity(design)) return(NULL)

    ea            <- eval_args
    ea$field_book <- design$field_book
    eff   <- tryCatch(do.call(met_evaluate_famoptg_efficiency, ea), error = function(e) NULL)
    score <- extract_score(eff)

    list(design = design, efficiency = eff, score = score, seed = s)
  }

  # ============================================================
  # SAFE BEST UPDATE
  # ============================================================
  # Updates running best only when score improves AND integrity re-passes.
  .update_best <- function(best, best_score, cand) {
    if (is.null(cand) || is.na(cand$score))
      return(list(best = best, best_score = best_score))
    if (cand$score >= best_score)
      return(list(best = best, best_score = best_score))
    if (!.check_famoptg_integrity(cand$design)) {
      warning(paste0(
        "[met_optimize_famoptg] Candidate passed initial integrity check but ",
        "failed re-check at best-update stage. Candidate discarded."
      ))
      return(list(best = best, best_score = best_score))
    }
    list(best = cand, best_score = cand$score)
  }

  # ============================================================
  # RANDOM RESTART
  # ============================================================
  if (verbose_opt) message(sprintf(
    "\n[met_optimize_famoptg] criterion = '%s' | n_restarts = %d | master seed = %d\n",
    criterion, n_restarts, master_seed
  ))

  best       <- NULL
  best_score <- Inf
  scores     <- rep(NA_real_, n_restarts)
  n_failed   <- 0L

  for (i in seq_len(n_restarts)) {
    s    <- sample.int(.Machine$integer.max, 1)
    cand <- run_pipeline(s)
    sc   <- if (is.null(cand)) NA_real_ else cand$score
    scores[i] <- sc

    if (is.null(cand) || is.na(sc)) {
      n_failed <- n_failed + 1L
    } else {
      upd        <- .update_best(best, best_score, cand)
      best       <- upd$best
      best_score <- upd$best_score
    }

    if (verbose_opt) message(sprintf(
      "  [RS] %d/%d | score = %s | best = %s%s",
      i, n_restarts,
      ifelse(is.na(sc), "FAILED", sprintf("%.6f", sc)),
      ifelse(is.infinite(best_score), "none yet", sprintf("%.6f", best_score)),
      if (is.null(cand)) " [construction/integrity failed]" else ""
    ))
  }

  # Report failure rate
  if (n_failed > 0L) {
    rate <- n_failed / n_restarts
    msg  <- sprintf(
      "[met_optimize_famoptg] %d/%d restarts failed construction or integrity check (%.0f%%).",
      n_failed, n_restarts, 100 * rate
    )
    if (rate > max_failure_rate) {
      stop(paste0(
        msg, "\nFailure rate exceeds threshold (",
        round(100 * max_failure_rate), "%). ",
        "Check n_blocks, n_rows, n_cols, p_rep_reps."
      ))
    } else {
      warning(msg)
    }
  }

  # ============================================================
  # FINAL INTEGRITY CHECK
  # ============================================================
  if (!is.null(best) && !is.null(best$design)) {
    if (!.check_famoptg_integrity(best$design)) {
      warning(paste0(
        "[met_optimize_famoptg] Stored best design failed final integrity check. ",
        "Triggering emergency fallback. Please report this as a bug."
      ))
      best <- NULL
    }
  }

  # ============================================================
  # FALLBACK
  # ============================================================
  if (is.null(best) || is.null(best$design)) {
    best <- .emergency_fallback()
  }

  # ============================================================
  # ASSEMBLE OUTPUT
  # ============================================================
  out            <- best$design
  out$efficiency <- best$efficiency

  # CDmean: convert internal negated score back to positive for reporting
  reported_best_score <- if (criterion == "CDmean" && !is.na(best$score))
    -best$score else best$score

  out$optimization <- list(
    method      = "RS",
    criterion   = criterion,
    best_score  = reported_best_score,
    score_history = if (criterion == "CDmean") {
      ifelse(is.na(scores), NA_real_, -scores)
    } else {
      scores
    },
    master_seed = master_seed,
    n_restarts  = n_restarts,
    n_failed    = n_failed
  )

  if (verbose_opt) {
    score_display <- ifelse(is.na(reported_best_score), NaN, reported_best_score)
    direction     <- if (criterion == "CDmean") "higher is better" else "lower is better"
    message(sprintf(
      "\n[met_optimize_famoptg] Done. Best %s = %.6f (%s) | valid restarts = %d/%d\n",
      criterion, score_display, direction,
      n_restarts - n_failed, n_restarts
    ))
  }

  out
}
