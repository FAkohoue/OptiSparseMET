#' Multi-trait index weights from desired genetic gains
#'
#' @description
#' Economic weights are hard to elicit and subjective. A more intuitive
#' alternative is to state *desired genetic gains* directly -- "raise yield by
#' 1.5 units while lowering disease score by 1 unit" -- and derive the index
#' weights that deliver gains in those proportions. `desired_gain_weights()`
#' returns those weights, ready to use as the multi-trait objective in
#' [optimize_design()] / [design_objective()] / [expected_genetic_gain()] (pass
#' them as `trait_weights` with `trait_gencov = gen_cov`).
#'
#' @details
#' Two sources are supported:
#' \describe{
#'   \item{`"pesek_baker"` (default, built-in)}{The classical desired-gain index
#'     of Pesek & Baker (1969): with genetic (co)variance matrix \eqn{G} and
#'     signed desired gains \eqn{d}, the weights are \eqn{b = G^{-1} d}. Because
#'     the expected genetic gains from index selection are proportional to
#'     \eqn{G b}, this makes the realised gains proportional to \eqn{d} exactly.
#'     Deterministic and dependency-free.}
#'   \item{`"dgqgsi"`}{Extracts the optimised index-weight vector from a result
#'     produced by the \pkg{DGQGSI} package (Akohoue 2026;
#'     `run_desired_gain_index_Joukhadar2024()` / `run_dgsi_qgsi_pipeline()`),
#'     which optimises desired-gain indices (and a quadratic genomic variant)
#'     from breeding data. Pass the result as `dgqgsi_result`; if the weight
#'     vector is stored under a non-standard name, give it with `weights_field`.}
#' }
#' Desired gains are supplied as positive magnitudes; traits named in
#' `lower_is_better` have their sign flipped internally (the same convention as
#' \pkg{DGQGSI}), so a smaller value counts as a gain.
#'
#' @param desired_gain Named (recommended) numeric vector of desired gain
#'   magnitudes, one per trait.
#' @param gen_cov Trait genetic (co)variance matrix; if it has dimnames it is
#'   aligned to `names(desired_gain)`.
#' @param lower_is_better Optional trait names (or indices) for which a decrease
#'   is the improvement; their desired gain is negated.
#' @param method `"pesek_baker"` (default) or `"dgqgsi"`.
#' @param dgqgsi_result For `method = "dgqgsi"`: a result object/list from the
#'   \pkg{DGQGSI} package.
#' @param weights_field Optional name of the element in `dgqgsi_result` holding
#'   the numeric weight vector, if auto-detection fails.
#' @return A list with `weights` (the index weights \eqn{b}), `sigma_index`
#'   (\eqn{\sqrt{b'Gb}}, the index genetic SD), `desired_gain` (signed),
#'   `gen_cov` (aligned), and `method`.
#' @references
#' Pesek, J., & Baker, R. J. (1969). Desired improvement in relation to selection
#' indices. *Canadian Journal of Plant Science*, 49, 803-804.
#'
#' Akohoue, F. (2026). DGQGSI: Desired-Gain Selection Index and Quadratic Genomic
#' Selection Index tools for breeding programs. R package.
#' @seealso [expected_genetic_gain()], [design_objective()], [optimize_design()].
#' @examples
#' G <- matrix(c(1, .3, -.2, .3, .8, .1, -.2, .1, .5), 3,
#'             dimnames = list(c("YLD", "DIS", "MOI"), c("YLD", "DIS", "MOI")))
#' dg <- c(YLD = 1.5, DIS = 1.0, MOI = 1.0)
#' w <- desired_gain_weights(dg, gen_cov = G, lower_is_better = "DIS")
#' w$weights
#' @export
desired_gain_weights <- function(desired_gain, gen_cov,
                                 lower_is_better = NULL,
                                 method = c("pesek_baker", "dgqgsi"),
                                 dgqgsi_result = NULL, weights_field = NULL) {
  method <- match.arg(method)
  d <- as.numeric(desired_gain)
  traits <- names(desired_gain)
  G <- as.matrix(gen_cov)

  if (!is.null(traits) && !is.null(colnames(G))) {
    if (!all(traits %in% colnames(G)))
      stop("`gen_cov` is missing traits: ",
           paste(setdiff(traits, colnames(G)), collapse = ", "), ".")
    G <- G[traits, traits, drop = FALSE]
  }
  if (nrow(G) != length(d) || ncol(G) != length(d))
    stop("`gen_cov` must be a length(desired_gain) square matrix.")

  d_signed <- d
  if (!is.null(lower_is_better)) {
    idx <- if (is.character(lower_is_better)) match(lower_is_better, traits)
           else as.integer(lower_is_better)
    idx <- idx[!is.na(idx)]
    if (length(idx)) d_signed[idx] <- -d_signed[idx]
  }

  b <- if (method == "dgqgsi")
    .extract_dgqgsi_weights(dgqgsi_result, weights_field, length(d))
  else
    # Pesek-Baker: b = G^{-1} d. Guard against a singular trait covariance
    # (e.g. perfectly correlated traits) with a symmetric pseudo-inverse.
    tryCatch(as.numeric(solve(G, d_signed)),
             error = function(e) as.numeric(.pinv_sym_dense(G) %*% d_signed))

  if (!is.null(traits)) names(b) <- traits
  sigma_index <- sqrt(as.numeric(t(b) %*% G %*% b))

  list(weights = b, sigma_index = sigma_index,
       desired_gain = stats::setNames(d_signed, traits),
       gen_cov = G, method = method)
}


# Pull the index-weight vector out of a DGQGSI result object.
.extract_dgqgsi_weights <- function(res, weights_field, n_traits) {
  if (is.null(res))
    stop("`method = 'dgqgsi'` needs `dgqgsi_result` (a DGQGSI output object).")
  res <- if (is.list(res) && !is.null(res$dg_result)) res$dg_result else res

  pick <- function(x, nm) if (is.list(x) && nm %in% names(x)) x[[nm]] else NULL

  cand <- weights_field
  if (is.null(cand))
    cand <- c("weights", "b", "index_weights", "b_weights",
              "index_coefficients", "index_vector", "coefficients")
  for (nm in cand) {
    v <- pick(res, nm)
    if (is.numeric(v) && length(v) == n_traits) return(as.numeric(v))
  }
  stop("Could not find a length-", n_traits, " numeric weight vector in the ",
       "DGQGSI result. Available names: ",
       paste(names(res), collapse = ", "),
       ". Pass the correct one via `weights_field`.")
}
