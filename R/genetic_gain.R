#' Standardised selection intensity
#'
#' @description
#' The selection intensity \eqn{i} for truncation selection of the top
#' proportion \eqn{p} of a normally distributed population, \eqn{i = \phi(z_p)/p}
#' with \eqn{z_p = \Phi^{-1}(1-p)}. It is the multiplier that converts prediction
#' accuracy into expected genetic gain in the breeder's equation.
#'
#' @param prop Selected proportion(s) in (0, 1).
#' @return Selection intensity, one value per `prop`.
#' @examples
#' selection_intensity(c(0.01, 0.05, 0.1, 0.2))
#' @export
selection_intensity <- function(prop) {
  if (any(prop <= 0 | prop >= 1))
    stop("`prop` must be in (0, 1).")
  stats::dnorm(stats::qnorm(1 - prop)) / prop
}


#' Expected genetic gain of a design
#'
#' @description
#' Turns a design's prediction reliability into expected genetic gain per cycle
#' via the breeder's equation \eqn{\Delta G = i\,r\,\sigma}, where \eqn{i} is the
#' [selection_intensity()], \eqn{r=\sqrt{\text{reliability}}} is the accuracy of
#' the (across-TPE) breeding-value prediction, and \eqn{\sigma} is the genetic
#' standard deviation of the selection target. Reliability is taken from a
#' [met_information()] result (its `CDmean`, the mean coefficient of
#' determination) or supplied directly.
#'
#' @details
#' **Single trait.** \eqn{\sigma = \sqrt{\sigma_g^2}} and
#' \eqn{\Delta G = i\sqrt{r^2}\,\sigma_g}. Note that for one trait \eqn{\Delta G}
#' is a monotone transform of reliability, so maximising gain and maximising
#' CDmean are equivalent objectives; the gain scale becomes distinct only in the
#' multi-trait case or when comparing targets with different \eqn{\sigma_g}.
#'
#' **Multi-trait index.** With index weights \eqn{w} and a genetic (co)variance
#' matrix among traits \eqn{\Sigma_T}, the aggregate genotype is \eqn{H = w'g}
#' with standard deviation \eqn{\sigma_I = \sqrt{w'\Sigma_T w}}, and
#' \eqn{\Delta H = i\,r\,\sigma_I}. Rather than subjective *economic* weights, the
#' recommended way to obtain \eqn{w} is from breeder-defined **desired gains**
#' via [desired_gain_weights()] (classical Pesek-Baker, or the \pkg{DGQGSI}
#' package). This uses the design's reliability for the index and assumes the
#' traits share the design (measured on the same plots) -- a first-order
#' approximation, not a full multi-trait mixed-model PEV.
#'
#' @param reliability Reliability (squared accuracy) in \[0, 1], or `NULL` to
#'   take it from `info`.
#' @param info Optional [met_information()] result; its `CDmean` is used as the
#'   reliability when `reliability` is `NULL`.
#' @param sigma_g Genetic standard deviation of the (single) selection target.
#'   Default 1.
#' @param prop Selected proportion. Default 0.1.
#' @param trait_weights Optional index weights for a multi-trait target,
#'   typically from [desired_gain_weights()] (desired gains, not economic
#'   weights).
#' @param trait_gencov Optional trait genetic (co)variance matrix
#'   (`length(trait_weights)` square) for the index; if `NULL` the identity is
#'   used.
#' @return A list with `gain` (\eqn{\Delta G} or \eqn{\Delta H}), `intensity`,
#'   `accuracy` (\eqn{r}), and `sigma` (the target SD used).
#' @seealso [met_information()], [design_objective()].
#' @examples
#' expected_genetic_gain(reliability = 0.5, sigma_g = 1, prop = 0.1)$gain
#' expected_genetic_gain(reliability = 0.5, prop = 0.1,
#'                       trait_weights = c(2, -1),
#'                       trait_gencov = matrix(c(1, .3, .3, .5), 2))$gain
#' @export
expected_genetic_gain <- function(reliability = NULL, info = NULL,
                                  sigma_g = 1, prop = 0.1,
                                  trait_weights = NULL, trait_gencov = NULL) {
  if (is.null(reliability)) {
    if (is.null(info) || is.null(info$CDmean))
      stop("Supply `reliability` or a met_information() `info` with CDmean.")
    reliability <- info$CDmean
  }
  reliability <- max(0, min(1, reliability))
  i  <- selection_intensity(prop)
  r  <- sqrt(reliability)

  if (!is.null(trait_weights)) {
    w  <- as.numeric(trait_weights)
    ST <- if (is.null(trait_gencov)) diag(length(w)) else as.matrix(trait_gencov)
    if (nrow(ST) != length(w) || ncol(ST) != length(w))
      stop("`trait_gencov` must be length(trait_weights) square.")
    sigma <- sqrt(as.numeric(t(w) %*% ST %*% w))
  } else {
    sigma <- sqrt(sigma_g^2)
  }

  list(gain = i * r * sigma, intensity = i, accuracy = r, sigma = sigma)
}
