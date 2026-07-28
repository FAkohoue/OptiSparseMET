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
# and the requested environment loads (column sums, which may be unequal). Only
# the pairwise concurrence /
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
                                iter = 2000L, max_line_pair_n = 800L,
                                groups = NULL, seed_budget = NULL,
                                environment_cost = NULL) {
  kind     <- match.arg(kind)
  n_sparse <- nrow(S)
  n_env    <- ncol(S)
  if (n_sparse < 2L || n_env < 2L) return(S)
  if (!is.null(groups)) {
    if (is.null(names(groups)) || !all(rownames(S) %in% names(groups)))
      stop("`groups` must be named and cover every row of `S`.")
    groups <- as.character(groups[rownames(S)])
  }
  seed_constrained <- !is.null(seed_budget) || !is.null(environment_cost)
  current_seed <- NULL
  if (seed_constrained) {
    if (is.null(seed_budget) || is.null(environment_cost))
      stop("Supply both `seed_budget` and `environment_cost`.")
    if (is.null(names(seed_budget)) ||
        !all(rownames(S) %in% names(seed_budget)))
      stop("`seed_budget` must be named and cover every row of `S`.")
    if (is.null(names(environment_cost)) ||
        !all(colnames(S) %in% names(environment_cost)))
      stop("`environment_cost` must be named and cover every column of `S`.")
    seed_budget <- seed_budget[rownames(S)]
    environment_cost <- environment_cost[colnames(S)]
    current_seed <- as.numeric(S %*% environment_cost)
  }

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
    if (!is.null(groups)) B <- B[groups[B] == groups[a]]
    if (length(B) == 0L) next
    b <- if (length(B) == 1L) B else sample(B, 1L)
    if (seed_constrained) {
      new_a <- current_seed[a] - environment_cost[e1] + environment_cost[e2]
      new_b <- current_seed[b] - environment_cost[e2] + environment_cost[e1]
      if (new_a > seed_budget[a] + 1e-8 ||
          new_b > seed_budget[b] + 1e-8) next
    }

    S[a, e1] <- 0L; S[a, e2] <- 1L
    S[b, e2] <- 0L; S[b, e1] <- 1L
    new <- .balance_objective(S, kind)

    if (new <= cur) {
      cur <- new
      if (seed_constrained) {
        current_seed[a] <- new_a
        current_seed[b] <- new_b
      }
    } else {
      S[a, e1] <- 1L; S[a, e2] <- 0L            # revert
      S[b, e2] <- 1L; S[b, e1] <- 0L
    }
    if (cur == 0) break
  }
  S
}

# Correlation-adaptive environment-pair refinement. Connectivity is the number
# of distinct treatments shared by two environments. Replication within an
# environment is deliberately absent from this calculation.
.adaptive_pair_refinement <- function(
    S, pair_targets, common_offset = 0L,
    aggregate = c("maximin", "cvar", "mean"), cvar_alpha = 0.25,
    iter = 2000L, groups = NULL, seed_budget = NULL,
    environment_cost = NULL) {
  aggregate <- match.arg(aggregate)
  n_sparse <- nrow(S)
  n_env <- ncol(S)
  n_pairs <- n_env * (n_env - 1L) / 2L
  if (!is.numeric(pair_targets) || length(pair_targets) != n_pairs ||
      any(!is.finite(pair_targets)) || any(pair_targets <= 0))
    stop("`pair_targets` must contain one finite positive target per ",
         "environment pair.")
  if (!is.numeric(common_offset) || length(common_offset) != 1L ||
      !is.finite(common_offset) || common_offset < 0)
    stop("`common_offset` must be one finite non-negative value.")
  if (!is.numeric(cvar_alpha) || length(cvar_alpha) != 1L ||
      !is.finite(cvar_alpha) || cvar_alpha <= 0 || cvar_alpha > 1)
    stop("`cvar_alpha` must be in (0, 1].")

  pair_fraction <- function(M) {
    shared <- .env_pair_offdiag(M) + common_offset
    pmin(1, shared / pair_targets)
  }
  lower_tail_mean <- function(x, alpha) {
    x <- sort(as.numeric(x))
    mass <- rep(1 / length(x), length(x))
    take <- pmin(
      mass,
      pmax(0, alpha - c(0, head(cumsum(mass), -1L)))
    )
    sum(x * take) / alpha
  }
  score <- function(M) {
    fraction <- pair_fraction(M)
    primary <- switch(
      aggregate,
      maximin = min(fraction),
      cvar = lower_tail_mean(fraction, cvar_alpha),
      mean = mean(fraction)
    )
    # The secondary mean resolves flat maximin/CVaR steps without permitting a
    # loss in the primary risk criterion.
    c(primary = primary, secondary = mean(fraction))
  }
  is_not_worse <- function(new, old, tol = 1e-12) {
    new[1L] > old[1L] + tol ||
      (abs(new[1L] - old[1L]) <= tol &&
         new[2L] >= old[2L] - tol)
  }

  before_fraction <- pair_fraction(S)
  before_score <- score(S)
  if (n_sparse < 2L || n_env < 2L || iter < 1L)
    return(list(
      allocation = S, before_fraction = before_fraction,
      after_fraction = before_fraction,
      before_score = before_score, after_score = before_score
    ))

  if (!is.null(groups)) {
    if (is.null(names(groups)) || !all(rownames(S) %in% names(groups)))
      stop("`groups` must be named and cover every sparse treatment.")
    groups <- as.character(groups[rownames(S)])
  }
  seed_constrained <- !is.null(seed_budget) || !is.null(environment_cost)
  current_seed <- NULL
  if (seed_constrained) {
    if (is.null(seed_budget) || is.null(environment_cost))
      stop("Supply both `seed_budget` and `environment_cost`.")
    seed_budget <- seed_budget[rownames(S)]
    environment_cost <- environment_cost[colnames(S)]
    current_seed <- as.numeric(S %*% environment_cost)
  }

  current_score <- before_score
  for (i in seq_len(as.integer(iter))) {
    es <- sample.int(n_env, 2L)
    e1 <- es[1L]
    e2 <- es[2L]
    move_12 <- which(S[, e1] == 1L & S[, e2] == 0L)
    move_21 <- which(S[, e2] == 1L & S[, e1] == 0L)
    if (!length(move_12) || !length(move_21)) next
    a <- if (length(move_12) == 1L) move_12 else sample(move_12, 1L)
    if (!is.null(groups))
      move_21 <- move_21[groups[move_21] == groups[a]]
    if (!length(move_21)) next
    b <- if (length(move_21) == 1L) move_21 else sample(move_21, 1L)

    if (seed_constrained) {
      new_a <- current_seed[a] - environment_cost[e1] + environment_cost[e2]
      new_b <- current_seed[b] - environment_cost[e2] + environment_cost[e1]
      if (new_a > seed_budget[a] + 1e-8 ||
          new_b > seed_budget[b] + 1e-8) next
    }

    S[a, e1] <- 0L
    S[a, e2] <- 1L
    S[b, e2] <- 0L
    S[b, e1] <- 1L
    proposed_score <- score(S)
    if (is_not_worse(proposed_score, current_score)) {
      current_score <- proposed_score
      if (seed_constrained) {
        current_seed[a] <- new_a
        current_seed[b] <- new_b
      }
    } else {
      S[a, e1] <- 1L
      S[a, e2] <- 0L
      S[b, e2] <- 1L
      S[b, e1] <- 0L
    }
  }

  list(
    allocation = S,
    before_fraction = before_fraction,
    after_fraction = pair_fraction(S),
    before_score = before_score,
    after_score = current_score
  )
}
