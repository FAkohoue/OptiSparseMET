#' Select a representative subset of environments from a TPE (decision 1)
#'
#' @description
#' `select_environments()` chooses `n` environments to represent a target
#' population of environments (TPE), given an environment relationship matrix
#' `D` (e.g. genetic correlations across environments, or a covariance derived
#' from enviromic data). It is a corrected, modified version of the
#' environment-selection idea in Colmant et al. (2026, decision 1).
#'
#' @details
#' Methods:
#' \describe{
#'   \item{`"representative"` (default)}{Maximises coverage of the whole TPE:
#'     the selected set \eqn{S} maximises \eqn{\sum_e \max_{s \in S} D_{es}}
#'     (a k-medoid / facility-location criterion). This targets
#'     *representativeness* directly.}
#'   \item{`"optcontrib"`}{The modified optimal-contribution / spread criterion
#'     \eqn{-q'Dq}: greedily minimises the total similarity among selected
#'     environments, maximising diversity. Note (per the package critique) that
#'     spread is not the same as representativeness -- a maximally diverse set can
#'     be a set of mutually atypical outliers.}
#'   \item{`"kmeans"`, `"hclust"`}{Clustering baselines (as in the paper):
#'     partition environments into `n` groups and return the medoid of each.}
#'   \item{`"random"`}{Uniform random subset (control).}
#' }
#'
#' @param D Square environment relationship/similarity matrix (higher = more
#'   similar), with row/column names.
#' @param n Number of environments to select.
#' @param method One of `"representative"`, `"optcontrib"`, `"kmeans"`,
#'   `"hclust"`, `"random"`.
#' @param seed Optional RNG seed (used by `"kmeans"` and `"random"`).
#'
#' @return A list with `selected` (environment names), `method`, and `score`
#'   (coverage for `"representative"`, total within-set similarity for
#'   `"optcontrib"`; `NA` otherwise).
#'
#' @seealso [select_individuals()], [met_information()].
#' @export
select_environments <- function(D, n,
                                method = c("representative", "optcontrib",
                                           "kmeans", "hclust", "random"),
                                seed = NULL) {
  method <- match.arg(method)
  D <- as.matrix(D)
  E <- nrow(D)
  if (is.null(rownames(D))) rownames(D) <- colnames(D) <- paste0("E", seq_len(E))
  envs <- rownames(D)
  if (n < 1L || n > E) stop("`n` must be between 1 and nrow(D).")
  if (!is.null(seed)) set.seed(seed)

  coverage <- function(S) sum(apply(D[, S, drop = FALSE], 1, max))

  if (method == "representative") {
    S <- integer(0)
    for (step in seq_len(n)) {
      cand <- setdiff(seq_len(E), S)
      gains <- vapply(cand, function(j) coverage(c(S, j)), numeric(1))
      S <- c(S, cand[which.max(gains)])
    }
    score <- coverage(S)

  } else if (method == "optcontrib") {
    # Greedy spread: add the environment with least similarity to those chosen.
    start <- which(rowSums(D) == min(rowSums(D)))[1]
    S <- start
    while (length(S) < n) {
      cand <- setdiff(seq_len(E), S)
      addsim <- vapply(cand, function(j) sum(D[j, S]), numeric(1))
      S <- c(S, cand[which.min(addsim)])
    }
    score <- sum(D[S, S])

  } else if (method %in% c("kmeans", "hclust")) {
    if (method == "kmeans") {
      km  <- stats::kmeans(D, centers = n, nstart = 10)
      cl  <- km$cluster
    } else {
      dst <- stats::as.dist(max(D) - D)
      cl  <- stats::cutree(stats::hclust(dst), k = n)
    }
    S <- integer(0)
    for (g in sort(unique(cl))) {
      members <- which(cl == g)
      # medoid: member with highest mean similarity to its cluster
      medoid <- members[which.max(rowMeans(D[members, members, drop = FALSE]))]
      S <- c(S, medoid)
    }
    score <- NA_real_

  } else { # random
    S <- sample(seq_len(E), n)
    score <- NA_real_
  }

  list(selected = envs[S], method = method, score = score)
}
