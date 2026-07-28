#' Build a relationship matrix from markers, pedigree, or both
#'
#' @description
#' The prediction and evaluation tools in this package take a single relationship
#' matrix. `build_relationship_matrix()` constructs it from raw inputs:
#' \describe{
#'   \item{`"genomic"`}{A VanRaden (2008) genomic relationship matrix from a
#'     marker matrix coded 0/1/2.}
#'   \item{`"pedigree"`}{The numerator relationship matrix \eqn{A} from a
#'     three-column pedigree, via Henderson's recursion (captures inbreeding).}
#'   \item{`"hmatrix"`}{A single-step \eqn{H} matrix (Legarra et al. 2009;
#'     Aguilar et al. 2010) that blends the genomic matrix on the genotyped block
#'     into the pedigree matrix -- the standard way to use markers *and* pedigree
#'     together and to include non-genotyped lines (predicted through pedigree).}
#' }
#'
#' @param markers Numeric marker matrix (genotypes in rows, markers in columns),
#'   coded 0/1/2. Required for `"genomic"` and `"hmatrix"`.
#' @param pedigree Data frame with columns `id`, `sire`, `dam` (unknown parents
#'   as `0` or `NA`). Required for `"pedigree"` and `"hmatrix"`.
#' @param type `"genomic"`, `"pedigree"`, or `"hmatrix"`.
#' @param omega For `"hmatrix"`: weight on the pedigree block within the genomic
#'   block (small, for invertibility). Default 0.05.
#' @param tau For `"hmatrix"`: overall weight on the genomic information.
#' @param method Genomic estimator: `"vanraden"` (built-in, default) or
#'   `"AGHmatrix"` to use `AGHmatrix::Gmatrix()` when the package is installed
#'   (falls back to the built-in with a warning otherwise).
#' @param relationship `"additive"` (default) or `"dominance"` (a Vitezica-2013
#'   dominance genomic relationship, intended as an explicit option for
#'   hybrid/SCA effects). Dominance is never added automatically and should not
#'   be used for ordinary inbred-line additive analyses. AGHmatrix uses
#'   `"VanRaden"` vs `"Vitezica"` accordingly; the built-in dominance needs
#'   integer 0/1/2 genotypes and ploidy 2.
#' @param ploidy Ploidy passed to the VanRaden estimator / AGHmatrix. Default 2.
#' @param tuneup If `TRUE`, pass the result through [tune_relationship_matrix()]
#'   (ASRgenomics `G.tuneup()`) to guarantee a well-conditioned, invertible
#'   matrix. Default `FALSE`.
#' @param blend,bend,rcn Passed to [tune_relationship_matrix()] when
#'   `tuneup = TRUE`.
#'   Defaults are `TRUE`, `FALSE`, and `TRUE`, respectively.
#'
#' @return A named, symmetric relationship matrix.
#'
#' @references
#' VanRaden, P. M. (2008). Efficient methods to compute genomic predictions.
#' *Journal of Dairy Science*, 91, 4414-4423.
#' Legarra, A., Aguilar, I., & Misztal, I. (2009). A relationship matrix
#' including full pedigree and genomic information. *Journal of Dairy Science*,
#' 92, 4656-4663.
#'
#' @seealso [met_information()], [select_individuals()], [tune_common_treatments()].
#' @export
build_relationship_matrix <- function(markers = NULL, pedigree = NULL,
                                      type = c("genomic", "pedigree", "hmatrix"),
                                      omega = 0.05, tau = 1,
                                      method = c("vanraden", "AGHmatrix"),
                                      relationship = c("additive", "dominance"),
                                      ploidy = 2, tuneup = FALSE,
                                      blend = TRUE, bend = FALSE, rcn = TRUE) {
  type <- match.arg(type)
  method <- match.arg(method)
  relationship <- match.arg(relationship)
  if (!is.numeric(ploidy) || length(ploidy) != 1L ||
      !is.finite(ploidy) || ploidy <= 0)
    stop("`ploidy` must be one finite positive number.")
  if (!is.numeric(omega) || length(omega) != 1L || !is.finite(omega) ||
      omega < 0 || omega > 1 ||
      !is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
      tau < 0 || tau > 1)
    stop("`omega` and `tau` must be finite values in [0, 1].")
  if (relationship == "dominance" && type != "genomic")
    stop("`relationship = 'dominance'` is only available for type = 'genomic'.")

  g_mat <- function(M) {
    M <- as.matrix(M)
    if (!is.numeric(M) || nrow(M) < 2L || ncol(M) < 1L ||
        any(!is.finite(M[!is.na(M)])))
      stop("`markers` must be a numeric genotype-by-marker matrix.")
    if (any(M < 0 | M > ploidy, na.rm = TRUE))
      stop("Marker dosages must lie between 0 and `ploidy`.")
    if (is.null(rownames(M)) || anyDuplicated(rownames(M)))
      stop("`markers` must have unique genotype row names.")
    p <- colMeans(M, na.rm = TRUE) / ploidy
    polymorphic <- is.finite(p) & p > 0 & p < 1
    if (!any(polymorphic))
      stop("No polymorphic markers remain after removing missing/monomorphic columns.")
    M <- M[, polymorphic, drop = FALSE]
    p <- p[polymorphic]
    q <- 1 - p

    use_agh <- method == "AGHmatrix"
    if (use_agh && relationship == "dominance" && anyNA(M)) {
      warning("Dominance markers contain missing values; using the built-in ",
              "centred-missing estimator instead of AGHmatrix.", call. = FALSE)
      use_agh <- FALSE
    }
    if (use_agh && relationship == "additive" && anyNA(M)) {
      for (j in seq_len(ncol(M)))
        M[is.na(M[, j]), j] <- ploidy * p[j]
    }
    if (use_agh) {
      if (!requireNamespace("AGHmatrix", quietly = TRUE)) {
        warning("Package 'AGHmatrix' not installed; using the built-in ",
                "estimator.", call. = FALSE)
      } else {
        agh <- if (relationship == "dominance") "Vitezica" else "VanRaden"
        G <- AGHmatrix::Gmatrix(SNPmatrix = M, method = agh, ploidy = ploidy)
        if (is.null(rownames(G))) dimnames(G) <- list(rownames(M), rownames(M))
        return(G)
      }
    }
    if (relationship == "dominance") {
      if (ploidy != 2)
        stop("Built-in dominance supports ploidy 2; install AGHmatrix for others.")
      if (any(!is.na(M) & !(M %in% 0:2)))
        stop("Built-in dominance needs integer 0/1/2 genotypes.")
      # Vitezica (2013) dominance coding: aa -> -2p^2, Aa -> 2pq, AA -> -2q^2
      H <- vapply(seq_len(ncol(M)), function(j)
        ifelse(M[, j] == 2, -2 * q[j]^2,
        ifelse(M[, j] == 1,  2 * p[j] * q[j],
        ifelse(M[, j] == 0, -2 * p[j]^2, 0))), numeric(nrow(M)))
      H[is.na(H)] <- 0
      denom <- sum((2 * p * q)^2); if (denom <= 0) denom <- 1
      G <- tcrossprod(H) / denom
    } else {
      # Mean-impute missing dosages before centring. These cells then contribute
      # zero to the centred marker matrix, the standard VanRaden treatment.
      for (j in seq_len(ncol(M)))
        M[is.na(M[, j]), j] <- ploidy * p[j]
      Z <- sweep(M, 2, ploidy * p)
      denom <- ploidy * sum(p * (1 - p))
      G <- tcrossprod(Z) / denom
    }
    if (is.null(rownames(G))) dimnames(G) <- list(rownames(M), rownames(M))
    G
  }

  a_mat <- function(ped) {
    ped <- as.data.frame(ped)
    if (!all(c("id", "sire", "dam") %in% names(ped)))
      stop("`pedigree` must have columns id, sire, dam.")
    id <- as.character(ped$id)
    sire <- as.character(ped$sire); dam <- as.character(ped$dam)
    if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id))
      stop("Pedigree `id` values must be non-missing and unique.")
    sire[is.na(sire)] <- "0"; dam[is.na(dam)] <- "0"
    sire[!nzchar(sire)] <- "0"; dam[!nzchar(dam)] <- "0"
    missing_parents <- setdiff(unique(c(sire, dam)), c("0", id))
    if (length(missing_parents))
      stop("Pedigree parents absent from `id`: ",
           paste(missing_parents, collapse = ", "), ".")
    # order so parents precede offspring
    ord <- .pedigree_order(id, sire, dam)
    id <- id[ord]; sire <- sire[ord]; dam <- dam[ord]
    pos <- stats::setNames(seq_along(id), id)
    n <- length(id); A <- matrix(0, n, n, dimnames = list(id, id))
    for (i in seq_len(n)) {
      s <- pos[sire[i]]; d <- pos[dam[i]]
      s <- if (is.na(s)) 0L else s; d <- if (is.na(d)) 0L else d
      if (i > 1L) for (j in seq_len(i - 1L)) {
        asj <- if (s > 0L) A[s, j] else 0
        adj <- if (d > 0L) A[d, j] else 0
        A[i, j] <- A[j, i] <- 0.5 * (asj + adj)
      }
      asd <- if (s > 0L && d > 0L) A[s, d] else 0
      A[i, i] <- 1 + 0.5 * asd
    }
    A
  }

  finish <- function(M) {
    M <- as.matrix(M)
    if (!is.numeric(M) || nrow(M) != ncol(M) || any(!is.finite(M)) ||
        is.null(rownames(M)) || is.null(colnames(M)) ||
        anyDuplicated(rownames(M)) || anyDuplicated(colnames(M)) ||
        !setequal(rownames(M), colnames(M)))
      stop("Constructed relationship matrix is not finite, named, and square.")
    M <- M[rownames(M), rownames(M), drop = FALSE]
    M <- (M + t(M)) / 2
    ee <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    if (min(ee) < -1e-8 * max(abs(ee), 1))
      stop("Constructed relationship matrix is materially indefinite.")
    if (isTRUE(tuneup))
      tune_relationship_matrix(M, blend = blend, bend = bend, rcn = rcn)
    else M
  }

  if (type == "genomic") {
    if (is.null(markers)) stop("`markers` is required for type = 'genomic'.")
    return(finish(g_mat(markers)))
  }
  if (type == "pedigree") {
    if (is.null(pedigree)) stop("`pedigree` is required for type = 'pedigree'.")
    return(finish(a_mat(pedigree)))
  }

  ## single-step H
  if (is.null(markers) || is.null(pedigree))
    stop("`markers` and `pedigree` are both required for type = 'hmatrix'.")
  G <- g_mat(markers)
  A <- a_mat(pedigree)
  gen <- intersect(rownames(G), rownames(A))
  if (length(gen) < 1L) stop("No genotyped individuals are present in the pedigree.")
  ung <- setdiff(rownames(A), gen)
  A22 <- A[gen, gen, drop = FALSE]
  Gg  <- G[gen, gen, drop = FALSE]

  # Tune G to the scale of A22 (Vitezica et al. 2011): G* = a + b G matching the
  # mean diagonal and mean off-diagonal of A22.
  md <- function(M) mean(diag(M))
  mo <- function(M) if (nrow(M) < 2L) 0 else mean(M[upper.tri(M)])
  bA <- md(A22) - mo(A22); bG <- md(Gg) - mo(Gg)
  b  <- if (abs(bG) > 1e-8) bA / bG else 1
  a  <- mo(A22) - b * mo(Gg)
  Gs <- a + b * Gg
  H22 <- tau * ((1 - omega) * Gs + omega * A22) + (1 - tau) * A22

  if (length(ung) == 0L) { dimnames(H22) <- list(gen, gen); return(finish(H22)) }

  A11 <- A[ung, ung, drop = FALSE]
  A12 <- A[ung, gen, drop = FALSE]; A21 <- t(A12)
  A22i <- tryCatch(solve(A22), error = function(e) .pinv_sym_dense(A22))
  T21 <- A22i %*% A21
  H11 <- A11 + t(T21) %*% (H22 - A22) %*% T21
  H12 <- t(T21) %*% H22
  H <- matrix(0, nrow(A), ncol(A), dimnames = dimnames(A))
  H[ung, ung] <- H11; H[gen, gen] <- H22
  H[ung, gen] <- H12; H[gen, ung] <- t(H12)
  finish((H + t(H)) / 2)
}


# Topological order of a pedigree so every parent precedes its offspring.
.pedigree_order <- function(id, sire, dam) {
  n <- length(id); pos <- stats::setNames(seq_along(id), id)
  depth <- rep(NA_integer_, n)
  get_depth <- function(i, stack = integer(0)) {
    if (!is.na(depth[i])) return(depth[i])
    if (i %in% stack) stop("Pedigree contains a cycle.")
    ps <- c(pos[sire[i]], pos[dam[i]]); ps <- ps[!is.na(ps)]
    d <- if (!length(ps)) 0L else 1L + max(vapply(ps, get_depth, integer(1), stack = c(stack, i)))
    depth[i] <<- d; d
  }
  for (i in seq_len(n)) get_depth(i)
  order(depth, seq_len(n))
}
