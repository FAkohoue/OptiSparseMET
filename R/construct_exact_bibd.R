#' Construct an exact balanced incomplete block design (when one exists)
#'
#' @description
#' `construct_exact_bibd()` builds a strict balanced incomplete block design
#' (BIBD) -- every pair of treatments co-occurring in the same environment an
#' equal number of times (constant \eqn{\lambda}) -- by delegating to
#' `crossdes::find.BIB()` and verifying the balance property (constant pairwise
#' concurrence lambda and constant replication) directly from the concurrence
#' matrix.
#'
#' A strict BIBD exists only when the necessary conditions hold: the count
#' identity \eqn{J r = I k}, integrality of
#' \eqn{\lambda = r(k-1)/(J-1)}, and Fisher's inequality \eqn{I \ge J}. For
#' sparse multi-environment testing the number of treatments \eqn{J} usually far
#' exceeds the number of environments \eqn{I}, so Fisher's inequality fails and
#' no BIBD exists; in that regime use [allocate_sparse_met()] with
#' `allocation_method = "equireplicate"` and `balance = "env_pair"` /
#' `"line_pair"` for the achievable near-balanced design.
#'
#' @param treatments Character vector of treatment (line) IDs, length \eqn{J}.
#' @param environments Character vector of environment names, length \eqn{I}
#'   (the BIBD blocks).
#' @param k Integer. Number of treatments per environment (block size).
#'
#' @return A list with `allocation_matrix` (a \eqn{J \times I} 0/1 incidence
#'   matrix with treatments as rows and environments as columns), `lambda`, and
#'   `is_bibd` (logical; constant-lambda balance check). Errors if `crossdes` is not
#'   installed or if the parameters admit no BIBD.
#'
#' @seealso [allocate_sparse_met()] for the near-balanced equireplicate
#'   allocation used when a strict BIBD does not exist;
#'   [check_equireplicate_feasibility()] to test `strict_bibd_possible` first.
#'
#' @examples
#' \dontrun{
#' ## Symmetric BIBD: v = b = 7, r = k = 3, lambda = 1 (Fano plane)
#' bibd <- construct_exact_bibd(
#'   treatments   = paste0("L", 1:7),
#'   environments = paste0("E", 1:7),
#'   k            = 3
#' )
#' bibd$is_bibd     # TRUE
#' bibd$lambda      # 1
#' }
#'
#' @export
construct_exact_bibd <- function(treatments, environments, k) {
  if (!requireNamespace("crossdes", quietly = TRUE))
    stop("Package 'crossdes' is required for construct_exact_bibd(). ",
         "Install it, or use allocate_sparse_met(allocation_method = ",
         "\"equireplicate\") for the near-balanced design.")

  treatments   <- unique(as.character(treatments))
  environments <- unique(as.character(environments))
  J <- length(treatments)
  I <- length(environments)
  k <- as.integer(k)

  if (J < 2L || I < 2L) stop("Need at least two treatments and two environments.")
  if (k < 1L || k >= J)
    stop("`k` (treatments per environment) must satisfy 1 <= k < J.")

  # Necessary BIBD conditions.
  if ((I * k) %% J != 0L)
    stop("No BIBD: I*k must be divisible by J (r = I*k/J must be an integer).")
  r <- (I * k) %/% J
  lambda_num <- r * (k - 1L)
  if (lambda_num %% (J - 1L) != 0L)
    stop("No BIBD: lambda = r(k-1)/(J-1) is not an integer for these parameters.")
  lambda <- lambda_num %/% (J - 1L)
  if (I < J)
    stop("No BIBD: Fisher's inequality requires #environments (", I,
         ") >= #treatments (", J, "). Use allocate_sparse_met() instead.")

  # crossdes::find.BIB(trt, b, k): treatments, blocks, block size.
  design <- crossdes::find.BIB(trt = J, b = I, k = k)

  # `design` is an I x k matrix: row = environment (block), entries = treatment
  # indices in that block. Convert to a J x I incidence matrix.
  alloc <- matrix(0L, nrow = J, ncol = I,
                  dimnames = list(treatments, environments))
  for (b in seq_len(I))
    alloc[design[b, ], b] <- 1L

  # Verify the BIBD balance property directly from the concurrence matrix
  # (constant off-diagonal lambda, constant replication r), rather than relying
  # on a specific crossdes verifier name.
  concurrence <- tcrossprod(alloc)                 # J x J
  offdiag     <- concurrence[upper.tri(concurrence)]
  is_bibd     <- all(rowSums(alloc) == r) && length(unique(offdiag)) == 1L &&
                 offdiag[1] == lambda

  list(
    allocation_matrix = alloc,
    lambda            = lambda,
    r                 = r,
    is_bibd           = isTRUE(is_bibd)
  )
}
