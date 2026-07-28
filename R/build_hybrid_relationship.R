#' Genomic relationship among hybrids/testcrosses from parental G
#'
#' @description
#' Hybrids (single crosses, testcrosses) are usually not genotyped directly --
#' only their parents are. And tools like \pkg{AGHmatrix} require integer 0/1/2
#' marker dosages, which a hybrid's parent-mean genotype (0, 0.5, 1, ...) is not.
#' `build_hybrid_relationship()` sidesteps both problems: build the **parental**
#' genomic relationship matrix from integer marker data (e.g.
#' `build_relationship_matrix(method = "AGHmatrix")`), then derive the additive
#' genomic relationship among the hybrids from the cross design.
#'
#' @details
#' For a single cross \eqn{H = a \times b}, the additive genomic value is the sum
#' of the parental gamete contributions, so the additive relationship between two
#' hybrids \eqn{a\times b} and \eqn{c\times d} is
#' \deqn{G_H[ab, cd] = \tfrac{1}{4}\,(G[a,c] + G[a,d] + G[b,c] + G[b,d]).}
#' This equals VanRaden's G computed on the parent-mean hybrid dosages **when
#' both are centred on the parental (base) allele frequencies** -- the
#' appropriate centring for hybrids (Technow et al. 2014). It is therefore not an
#' approximation: it derives the hybrid GRM on the parental frequency scale while
#' avoiding non-integer dosages in the estimator, rather than recomputing
#' allele frequencies from the specific set of hybrids. Applies to testcrosses
#' (one common tester/female) and full factorials alike.
#'
#' @param G_parents Parental genomic relationship matrix (row/column names are
#'   parent ids), e.g. from [build_relationship_matrix()].
#' @param crosses A data frame or matrix with two columns naming the two parents
#'   of each hybrid (both must be rows of `G_parents`).
#' @param hybrid_ids Optional names for the hybrids (length `nrow(crosses)`);
#'   defaults to `crosses` row names, else `parent1_x_parent2`.
#' @return The hybrid-by-hybrid additive genomic relationship matrix (dimnames =
#'   `hybrid_ids`). Pass it through [tune_relationship_matrix()] if you need a
#'   guaranteed-invertible version.
#' @references
#' Technow, F. et al. (2014). Genome properties and prospects of genomic
#' prediction of hybrid performance. *Genetics*, 197, 1343-1355.
#' @seealso [build_relationship_matrix()], [tune_relationship_matrix()].
#' @examples
#' set.seed(1)
#' Mpar <- matrix(sample(0:2, 6 * 300, TRUE), 6, 300,
#'                dimnames = list(c("A","B","C","D","T1","T2"), NULL))
#' Gp <- build_relationship_matrix(markers = Mpar, type = "genomic")
#' crosses <- data.frame(p1 = c("A","B","C","D"), p2 = c("T1","T1","T2","T2"))
#' Ghyb <- build_hybrid_relationship(Gp, crosses)
#' dim(Ghyb)
#' @export
build_hybrid_relationship <- function(G_parents, crosses, hybrid_ids = NULL) {
  Gp <- as.matrix(G_parents)
  if (!is.numeric(Gp) || nrow(Gp) != ncol(Gp) ||
      any(!is.finite(Gp)) || is.null(rownames(Gp)) ||
      is.null(colnames(Gp)) || anyDuplicated(rownames(Gp)) ||
      anyDuplicated(colnames(Gp)) ||
      !setequal(rownames(Gp), colnames(Gp)))
    stop("`G_parents` must be a finite, named square relationship matrix.")
  Gp <- Gp[rownames(Gp), rownames(Gp), drop = FALSE]
  if (!isTRUE(all.equal(Gp, t(Gp), tolerance = 1e-8)))
    stop("`G_parents` must be symmetric.")
  ee <- eigen((Gp + t(Gp)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ee) < -1e-8 * max(abs(ee), 1))
    stop("`G_parents` must be positive semidefinite.")
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  if (nrow(crosses) < 1L || ncol(crosses) < 2L)
    stop("`crosses` must have at least one row and two parent columns.")
  p1 <- as.character(crosses[[1]]); p2 <- as.character(crosses[[2]])
  if (anyNA(p1) || anyNA(p2) || any(!nzchar(p1)) || any(!nzchar(p2)))
    stop("Cross parent IDs must be non-missing.")
  miss <- setdiff(unique(c(p1, p2)), rownames(Gp))
  if (length(miss))
    stop("Parents not found in `G_parents`: ", paste(miss, collapse = ", "), ".")

  G <- 0.25 * (Gp[p1, p1, drop = FALSE] + Gp[p1, p2, drop = FALSE] +
               Gp[p2, p1, drop = FALSE] + Gp[p2, p2, drop = FALSE])
  automatic_rownames <- is.null(rownames(crosses)) ||
    identical(rownames(crosses), as.character(seq_len(nrow(crosses))))
  ids <- if (!is.null(hybrid_ids)) as.character(hybrid_ids)
         else if (!automatic_rownames) rownames(crosses)
         else paste(p1, p2, sep = "_x_")
  if (length(ids) != nrow(crosses) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids))
    stop("Hybrid IDs must be non-missing and unique.")
  dimnames(G) <- list(ids, ids)
  (G + t(G)) / 2
}


#' Combine relationship matrices (e.g. additive + dominance) for total merit
#'
#' @description
#' When more than one genetic effect matters -- classically **additive** and
#' **dominance** for hybrids -- the relevant relationship for predicting *total*
#' genetic merit is a variance-weighted combination of the component matrices.
#' `combine_relationship_matrices()` aligns the matrices by name and returns
#' \eqn{\sum_k w_k K_k / \sum_k w_k}, where the weights are the variance
#' proportions of each component (e.g. `c(additive = sigma2_A, dominance =
#' sigma2_D)`). The result is a single relationship matrix that the design/
#' evaluation functions can use so the design accounts for dominance, not just
#' additive effects.
#'
#' @param matrices A named list of relationship matrices over the same
#'   individuals (e.g. `list(additive = G_A, dominance = G_D)`).
#' @param weights One finite non-negative weight per matrix (variance components
#'   or proportions), normalised to sum to 1. Named weights are aligned to a
#'   named `matrices` list. Defaults to equal weights.
#' @return The combined relationship matrix (dimnames from the first matrix).
#' @seealso [build_relationship_matrix()], [build_hybrid_relationship()],
#'   [tune_relationship_matrix()].
#' @examples
#' set.seed(1)
#' M <- matrix(sample(0:2, 20 * 300, TRUE), 20, 300); rownames(M) <- paste0("H", 1:20)
#' G_A <- build_relationship_matrix(markers = M, type = "genomic")
#' G_D <- build_relationship_matrix(markers = M, type = "genomic",
#'                                  relationship = "dominance")
#' G_tot <- combine_relationship_matrices(list(additive = G_A, dominance = G_D),
#'                                        weights = c(2, 1))   # sigma2_A : sigma2_D
#' @export
combine_relationship_matrices <- function(matrices, weights = NULL) {
  if (!is.list(matrices) || !length(matrices)) stop("`matrices` is empty.")
  mats <- lapply(matrices, function(m) {
    m <- as.matrix(m)
    if (!is.numeric(m) || nrow(m) != ncol(m) || any(!is.finite(m)) ||
        is.null(rownames(m)) || is.null(colnames(m)) ||
        anyDuplicated(rownames(m)) || anyDuplicated(colnames(m)) ||
        !setequal(rownames(m), colnames(m)))
      stop("Every relationship matrix must be finite, named, and square.")
    m <- m[rownames(m), rownames(m), drop = FALSE]
    if (!isTRUE(all.equal(m, t(m), tolerance = 1e-8)))
      stop("Every relationship matrix must be symmetric.")
    ee <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    if (min(ee) < -1e-8 * max(abs(ee), 1))
      stop("Every relationship matrix must be positive semidefinite.")
    m
  })
  ids <- rownames(mats[[1]])
  mats <- lapply(mats, function(m) {
    if (!setequal(ids, rownames(m)))
      stop("All matrices must cover exactly the same individuals.")
    m[ids, ids, drop = FALSE]
  })
  if (is.null(weights)) weights <- rep(1, length(mats))
  if (!is.null(names(weights)) && any(names(weights) != "")) {
    matrix_names <- names(matrices)
    if (is.null(matrix_names) || any(matrix_names == "") ||
        anyDuplicated(matrix_names) || anyDuplicated(names(weights)) ||
        !setequal(matrix_names, names(weights)))
      stop("Named `weights` require unique names matching the named `matrices` list.")
    weights <- weights[matrix_names]
  }
  if (!is.numeric(weights) || length(weights) != length(mats) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) == 0)
    stop("`weights` must be finite, non-negative, have one value per matrix, and not all be zero.")
  weights <- weights / sum(weights)
  Reduce(`+`, Map(function(m, w) w * m, mats, weights))
}
