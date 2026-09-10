#' Create a validated sparse-MET design object
#'
#' @description
#' Packages the allocation, realised integer replication, fieldbooks, covariance
#' assumptions, and optimisation provenance into one versioned object. The
#' constructor runs cross-component QA so a design cannot silently disagree
#' with its plantable fieldbooks.
#'
#' @param allocation_matrix Named treatment-by-environment 0/1 matrix.
#' @param reps Optional named integer plot-count matrix. Defaults to one plot per
#'   allocated cell.
#' @param fieldbooks Optional named list of environment fieldbooks, each with a
#'   `Treatment` column. Treatments not in the allocation (for example checks)
#'   are allowed, but allocated treatments must reconcile exactly with `reps`.
#' @param G,Sigma_E Optional genetic and environment covariance matrices.
#' @param tpe_weights Optional target-population weights.
#' @param sigma_g2,sigma_e2 Variance assumptions.
#' @param diagnostics,provenance Named metadata lists.
#' @param validate Logical; run [validate_sparse_met_design()] immediately.
#' @return An object of class `sparse_met_design`.
#' @export
sparse_met_design <- function(
    allocation_matrix, reps = NULL, fieldbooks = NULL,
    G = NULL, Sigma_E = NULL, tpe_weights = NULL,
    sigma_g2 = 1, sigma_e2 = 1,
    diagnostics = list(), provenance = list(), validate = TRUE) {
  M <- as.matrix(allocation_matrix)
  if (is.null(reps)) reps <- M
  out <- structure(list(
    schema_version = "1.0.0",
    allocation_matrix = M,
    reps = as.matrix(reps),
    fieldbooks = fieldbooks,
    G = G, Sigma_E = Sigma_E,
    tpe_weights = if (is.null(colnames(M))) tpe_weights else
      .normalise_tpe_weights(tpe_weights, colnames(M)),
    sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
    diagnostics = diagnostics, provenance = provenance),
    class = "sparse_met_design")
  if (isTRUE(validate)) validate_sparse_met_design(out, error = TRUE)
  out
}

#' Validate a sparse-MET design object
#'
#' @param x A `sparse_met_design` object.
#' @param error If `TRUE`, stop on failure; otherwise return an audit list.
#' @return Invisibly `TRUE` when valid, or an audit list when `error = FALSE`.
#' @export
validate_sparse_met_design <- function(x, error = TRUE) {
  issues <- character(0)
  add <- function(msg) issues <<- c(issues, msg)
  if (!inherits(x, "sparse_met_design")) add("Object does not inherit `sparse_met_design`.")
  M <- tryCatch(as.matrix(x$allocation_matrix), error = function(e) NULL)
  if (is.null(M) || !is.numeric(M) || anyNA(M) || !all(M %in% c(0, 1)) ||
      is.null(rownames(M)) || is.null(colnames(M)) ||
      anyDuplicated(rownames(M)) || anyDuplicated(colnames(M))) {
    add("Allocation must be a uniquely named finite 0/1 matrix.")
  } else {
    if (!nrow(M) || !ncol(M) || any(colSums(M) == 0L))
      add("Allocation must be non-empty and every environment must contain a treatment.")
  }
  R <- tryCatch(as.matrix(x$reps), error = function(e) NULL)
  if (!is.null(M) &&
      (is.null(R) || !identical(dim(R), dim(M)) || !is.numeric(R) ||
       any(!is.finite(R)) || any(R < 0) || any(abs(R - round(R)) > 1e-8) ||
       any(R[M > 0] < 1) || any(R[M == 0] != 0)))
    add("Replication must be an integer matrix positive exactly in allocated cells.")
  if (!is.null(M) && !is.null(R) &&
      (!identical(rownames(R), rownames(M)) || !identical(colnames(R), colnames(M))))
    add("Replication dimnames must exactly match the allocation.")

  if (!is.null(x$fieldbooks) && !is.null(M) && !is.null(R)) {
    if (!is.list(x$fieldbooks) || is.null(names(x$fieldbooks)) ||
        anyDuplicated(names(x$fieldbooks)) ||
        !all(colnames(M) %in% names(x$fieldbooks))) {
      add("Fieldbooks must be a uniquely named list covering every environment.")
    } else for (e in colnames(M)) {
      fb <- x$fieldbooks[[e]]
      if (!is.data.frame(fb) || !"Treatment" %in% names(fb)) {
        add(paste0("Fieldbook `", e, "` lacks a Treatment column.")); next
      }
      counts <- table(as.character(fb$Treatment))
      observed <- numeric(nrow(M)); names(observed) <- rownames(M)
      ids <- intersect(names(counts), rownames(M))
      observed[ids] <- as.numeric(counts[ids])
      if (!identical(as.integer(observed), as.integer(R[, e])))
        add(paste0("Fieldbook `", e, "` does not reconcile with replication."))
    }
  }

  check_cov <- function(A, ids, nm) {
    if (is.null(A)) return()
    A <- as.matrix(A)
    if (!is.numeric(A) || nrow(A) != ncol(A) || any(!is.finite(A)) ||
        is.null(rownames(A)) || is.null(colnames(A)) ||
        !all(ids %in% rownames(A)) || !all(ids %in% colnames(A))) {
      add(paste0(nm, " must be a finite named square matrix covering its IDs.")); return()
    }
    A <- A[ids, ids, drop = FALSE]
    if (!isTRUE(all.equal(A, t(A), tolerance = 1e-8)))
      add(paste0(nm, " must be symmetric."))
    else if (min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) <
             -1e-8 * max(abs(A), 1)) add(paste0(nm, " must be positive semidefinite."))
  }
  if (!is.null(M)) {
    check_cov(x$G, rownames(M), "G")
    check_cov(x$Sigma_E, colnames(M), "Sigma_E")
  }
  if (!is.list(x$diagnostics) || !is.list(x$provenance))
    add("Diagnostics and provenance must be lists.")
  valid <- !length(issues)
  if (!valid && isTRUE(error))
    stop("Invalid sparse MET design:\n- ", paste(issues, collapse = "\n- "),
         call. = FALSE)
  if (isTRUE(error)) invisible(TRUE) else
    list(valid = valid, issues = issues,
         n_errors = length(issues), schema_version = x$schema_version)
}

#' @export
print.sparse_met_design <- function(x, ...) {
  M <- x$allocation_matrix
  cat("Sparse MET design (schema ", x$schema_version, ")\n", sep = "")
  cat("  Treatments:", nrow(M), " Environments:", ncol(M),
      " Allocated cells:", sum(M), " Plots:", sum(x$reps), "\n")
  cat("  TPE weights:", paste(names(x$tpe_weights),
                               format(x$tpe_weights, digits = 3),
                               sep = "=", collapse = ", "), "\n")
  invisible(x)
}
