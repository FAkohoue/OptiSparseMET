#' Optimise entry-across-environment allocation for G x E prediction (decision 3)
#'
#' @description
#' `optimize_allocation_gxe()` refines a sparse allocation to maximise the
#' across-TPE prediction reliability (minimise mean PEV) under the coupled G x E
#' information matrix of [met_information()]. It is the corrected, reproducible
#' form of decision 3 in Colmant et al. (2026): instead of the unbounded spread
#' proxy \eqn{-q'Dq} (with \eqn{D = G \otimes E}) or the underspecified
#' connectedness objective functions, it optimises the quantity that actually
#' governs prediction accuracy.
#'
#' @details
#' Swaps exchange a treatment between two environments (treatment `a` from `e1`
#' to `e2`, treatment `b` from `e2` to `e1`), so every treatment keeps its
#' replication count and every environment keeps its size -- the equireplicate
#' margins are preserved. A swap is accepted when it lowers the across-TPE mean
#' PEV. Each evaluation solves the (\eqn{JE \times JE}) information matrix, so
#' keep `J * E` within `max_dim`.
#'
#' @inheritParams met_information
#' @param iter Number of swap iterations.
#' @param seed Optional RNG seed.
#'
#' @return A list with the refined `allocation_matrix`, `mean_PEV_before`,
#'   `mean_PEV_after`, and `CDmean_after`.
#'
#' @seealso [allocate_sparse_met()] for the initial equireplicate allocation;
#'   [met_information()] for the objective.
#' @export
optimize_allocation_gxe <- function(allocation_matrix, G, Sigma_E = NULL,
                                    sigma_g2 = 1, sigma_e2 = 1,
                                    env_efficiency = NULL,
                                    iter = 200L, seed = NULL, max_dim = 2000L) {
  M <- allocation_matrix
  if (is.null(rownames(M)) || is.null(colnames(M)))
    stop("`allocation_matrix` must have row and column names.")
  E <- ncol(M)
  if (E < 2L) stop("Need at least two environments.")

  score <- function(mat) {
    met_information(mat, G = G, Sigma_E = Sigma_E,
                    sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                    env_efficiency = env_efficiency,
                    target = "across_tpe", max_dim = max_dim)$mean_PEV
  }

  before <- score(M)
  cur <- before
  if (!is.null(seed)) set.seed(seed)

  for (i in seq_len(iter)) {
    es <- sample.int(E, 2L); e1 <- es[1L]; e2 <- es[2L]
    A <- which(M[, e1] == 1 & M[, e2] == 0)
    B <- which(M[, e2] == 1 & M[, e1] == 0)
    if (!length(A) || !length(B)) next
    a <- A[sample.int(length(A), 1L)]
    b <- B[sample.int(length(B), 1L)]
    M2 <- M
    M2[a, e1] <- 0; M2[a, e2] <- 1
    M2[b, e2] <- 0; M2[b, e1] <- 1
    new <- score(M2)
    if (new < cur) { M <- M2; cur <- new }
  }

  info <- met_information(M, G = G, Sigma_E = Sigma_E,
                          sigma_g2 = sigma_g2, sigma_e2 = sigma_e2,
                          env_efficiency = env_efficiency,
                          target = "across_tpe", max_dim = max_dim)
  list(allocation_matrix = M,
       mean_PEV_before   = before,
       mean_PEV_after    = cur,
       CDmean_after      = info$CDmean)
}
