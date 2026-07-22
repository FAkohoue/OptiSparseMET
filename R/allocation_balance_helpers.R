# ==============================================================================
# allocation_balance_helpers.R
#
# Internal helpers that refine an equireplicate sparse allocation toward a
# near-balanced (regular-graph) design without disturbing its margins.
#
# A "swap" exchanges a sparse treatment between two environments: treatment a
# moves from environment e1 to e2, and treatment b moves from e2 to e1. Because
# each environment loses one line and gains one, and each treatment stays in the
# same number of environments, swaps preserve BOTH equal replication (row sums)
# and equal environment size (column sums). Only the pairwise concurrence /
# connectivity structure changes.
#
# A strict balanced incomplete block design (constant pairwise concurrence
# lambda) generally cannot exist when the number of treatments greatly exceeds
# the number of environments (Fisher's inequality b >= v, and integrality of
# lambda = r(k-1)/(v-1)). These helpers therefore target the achievable
# relaxation: minimise the variance of the off-diagonal concurrences.
#   - kind = "env_pair":  off-diagonals of crossprod(S)   (n_env x n_env) --
#                         the number of lines shared by each pair of
#                         environments (cross-environment connectivity).
#   - kind = "line_pair": off-diagonals of tcrossprod(S)  (n_sparse x n_sparse)
#                         -- the number of environments shared by each pair of
#                         treatments (pairwise concurrence).
# ==============================================================================

# Off-diagonal environment-pair shared-line counts (upper triangle).
.env_pair_offdiag <- function(S) {
  M <- crossprod(S)                     # n_env x n_env
  M[upper.tri(M)]
}

# Off-diagonal line-pair concurrence counts (upper triangle).
.line_pair_offdiag <- function(S) {
  C <- tcrossprod(S)                    # n_sparse x n_sparse
  C[upper.tri(C)]
}

# Objective for the swap search: variance of the relevant off-diagonals
# (0 => perfectly balanced). Uses population variance so a single-pair design
# does not return NA.
.balance_objective <- function(S, kind) {
  off <- if (kind == "env_pair") .env_pair_offdiag(S) else .line_pair_offdiag(S)
  n <- length(off)
  if (n < 2L) return(0)
  m <- mean(off)
  sum((off - m)^2) / n
}

# Summary metrics reported back to the caller (before/after the refinement).
# The line-pair metrics require an n_sparse x n_sparse crossproduct, so they are
# skipped (returned NA) for large designs to avoid a costly allocation.
.balance_metrics <- function(S, line_pair_max_n = 2000L) {
  ep <- .env_pair_offdiag(S)
  do_lp <- nrow(S) <= line_pair_max_n
  lp <- if (do_lp) .line_pair_offdiag(S) else NULL
  list(
    env_pair_var     = if (length(ep) > 1L) stats::var(ep) else 0,
    env_pair_range   = if (length(ep)) range(ep) else c(NA, NA),
    line_pair_var    = if (do_lp && length(lp) > 1L) stats::var(lp) else NA_real_,
    line_pair_values = if (do_lp) sort(unique(lp)) else NA_integer_
  )
}

# Swap-improvement search (greedy: accept any non-worsening swap).
# S is an integer 0/1 matrix of dimension n_sparse x n_env. Returns the refined
# matrix with identical row and column sums.
.balance_allocation <- function(S, kind = c("env_pair", "line_pair"),
                                iter = 2000L, max_line_pair_n = 800L) {
  kind     <- match.arg(kind)
  n_sparse <- nrow(S)
  n_env    <- ncol(S)
  if (n_sparse < 2L || n_env < 2L) return(S)

  if (kind == "line_pair" && n_sparse > max_line_pair_n) {
    message("line_pair balancing skipped: n_sparse (", n_sparse,
            ") exceeds max_line_pair_n (", max_line_pair_n,
            "); use balance = 'env_pair' for large designs.")
    return(S)
  }

  cur <- .balance_objective(S, kind)
  if (cur == 0) return(S)

  for (i in seq_len(iter)) {
    es <- sample.int(n_env, 2L)
    e1 <- es[1L]; e2 <- es[2L]
    A <- which(S[, e1] == 1L & S[, e2] == 0L)   # candidates to move e1 -> e2
    B <- which(S[, e2] == 1L & S[, e1] == 0L)   # candidates to move e2 -> e1
    if (length(A) == 0L || length(B) == 0L) next
    a <- if (length(A) == 1L) A else sample(A, 1L)
    b <- if (length(B) == 1L) B else sample(B, 1L)

    S[a, e1] <- 0L; S[a, e2] <- 1L
    S[b, e2] <- 0L; S[b, e1] <- 1L
    new <- .balance_objective(S, kind)

    if (new <= cur) {
      cur <- new
    } else {
      S[a, e1] <- 1L; S[a, e2] <- 0L            # revert
      S[b, e2] <- 1L; S[b, e1] <- 0L
    }
    if (cur == 0) break
  }
  S
}
