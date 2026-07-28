#' Make a relationship matrix well-conditioned and invertible
#'
#' @description
#' Genomic relationship matrices are frequently singular or ill-conditioned
#' (more markers than individuals, clones/near-duplicates, redundant entries),
#' which breaks the mixed-model equations. `tune_relationship_matrix()` returns a
#' tuned, invertible version. When the \pkg{ASRgenomics} package is available it
#' uses its `G.tuneup()` (blending toward an identity/`A` matrix and/or bending
#' negative eigenvalues, optionally targeting a reciprocal condition number) --
#' the standard, expert-trusted procedure. Otherwise it applies a simple
#' diagonal blend so the pipeline still runs.
#'
#' @param G A relationship matrix (e.g. from [build_relationship_matrix()]).
#' @param method How to tune. `"auto"` (default) uses `ASRgenomics::G.tuneup()`
#'   if the package is installed, otherwise the native `"bend"`. `"asrgenomics"`
#'   forces ASRgenomics. `"bend"` is a native nearest-positive-definite
#'   correction (eigenvalue flooring; needs no external package). `"blend"` is a
#'   native diagonal blend.
#' @param blend,bend,rcn Passed to `ASRgenomics::G.tuneup()` when it is used.
#' @param blend_weight Weight of the identity blend for the native `"blend"`
#'   method. Default 0.02.
#' @param eps Eigenvalue floor (as a fraction of the largest eigenvalue) for the
#'   native `"bend"` method. Default 1e-6.
#' @param ... Further arguments passed to `ASRgenomics::G.tuneup()`.
#' @return The tuned relationship matrix (with the original dimnames).
#' @references
#' Gezan, S. A., Borgognone, M. G., et al. ASRgenomics. VSN International.
#' Higham, N. J. (2002). Computing the nearest correlation matrix.
#' @seealso [build_relationship_matrix()], [kinship_pca()].
#' @examples
#' set.seed(1)
#' M <- matrix(sample(0:2, 30 * 200, TRUE), 30, 200)
#' rownames(M) <- paste0("L", 1:30)
#' G <- build_relationship_matrix(markers = M, type = "genomic")
#' Gb <- tune_relationship_matrix(G)                       # auto
#' Gb2 <- tune_relationship_matrix(G, method = "bend")     # no external package
#' @export
tune_relationship_matrix <- function(G,
                                     method = c("auto", "asrgenomics",
                                                "bend", "blend"),
                                     blend = TRUE, bend = FALSE, rcn = TRUE,
                                     blend_weight = 0.02, eps = 1e-6, ...) {
  method <- match.arg(method)
  G <- as.matrix(G)
  if (method == "auto")
    method <- if (requireNamespace("ASRgenomics", quietly = TRUE))
      "asrgenomics" else "bend"

  if (method == "asrgenomics") {
    if (!requireNamespace("ASRgenomics", quietly = TRUE))
      stop("Package 'ASRgenomics' is not installed; use method = 'bend' or 'blend'.")
    Gb <- ASRgenomics::G.tuneup(G = G, blend = blend, bend = bend, rcn = rcn, ...)$Gb
    if (is.null(dimnames(Gb))) dimnames(Gb) <- dimnames(G)
    return(Gb)
  }

  if (method == "bend") {
    e <- eigen((G + t(G)) / 2, symmetric = TRUE)
    floor_val <- eps * max(e$values)
    vals <- pmax(e$values, floor_val)
    Gb <- e$vectors %*% (vals * t(e$vectors))     # V diag(vals) V'
    dimnames(Gb) <- dimnames(G)
    return((Gb + t(Gb)) / 2)
  }

  ## native diagonal blend
  n <- nrow(G)
  Gb <- (1 - blend_weight) * G + blend_weight * diag(mean(diag(G)), n)
  dimnames(Gb) <- dimnames(G)
  Gb
}


#' Principal components of a kinship/relationship matrix
#'
#' @description
#' Population-structure principal components from a relationship matrix, useful
#' for structure-aware allocation grouping ([derive_allocation_groups()]) and
#' diagnostics. Uses `ASRgenomics::kinship.pca()` when available, otherwise a
#' symmetric eigen-decomposition fallback.
#'
#' @param K A kinship/relationship matrix (ideally a tuned, PD one).
#' @param ncp Number of principal components to return. Default: all.
#' @return A list with `values` (eigenvalues), `scores` (individual PC scores),
#'   and `prop_var` (proportion of variance per PC). When \pkg{ASRgenomics} is
#'   used its native object is returned instead.
#' @seealso [tune_relationship_matrix()], [derive_allocation_groups()].
#' @examples
#' set.seed(1)
#' M <- matrix(sample(0:2, 30 * 200, TRUE), 30, 200); rownames(M) <- paste0("L", 1:30)
#' G <- tune_relationship_matrix(build_relationship_matrix(markers = M, type = "genomic"))
#' pcs <- kinship_pca(G, ncp = 5)
#' @export
kinship_pca <- function(K, ncp = NULL) {
  K <- as.matrix(K)
  if (requireNamespace("ASRgenomics", quietly = TRUE)) {
    return(ASRgenomics::kinship.pca(K = K, ncp = if (is.null(ncp)) nrow(K) else ncp))
  }
  e <- eigen((K + t(K)) / 2, symmetric = TRUE)
  vals <- pmax(e$values, 0)
  k <- if (is.null(ncp)) length(vals) else min(ncp, length(vals))
  scores <- e$vectors[, seq_len(k), drop = FALSE] %*% diag(sqrt(vals[seq_len(k)]),
                                                           k, k)
  rownames(scores) <- rownames(K); colnames(scores) <- paste0("PC", seq_len(k))
  list(values = vals[seq_len(k)], scores = scores,
       prop_var = vals[seq_len(k)] / sum(vals))
}
