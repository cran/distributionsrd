#' The Exponential distribution
#'
#' Raw moments for the exponential distribution.
#' @param truncation lower truncation parameter, defaults to 0.
#' @param r rth raw moment of the distribution, defaults to 1.
#' @param rate rate of the distribution with default values of 1.
#' @param lower.tail logical; if TRUE (default), moments are \eqn{E[x^r|X \le y]}, otherwise, \eqn{E[x^r|X > y]}
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) = \frac{1}{s}e^{-\frac{\omega}{s}} , \qquad F_X(x) = 1-e^{-\frac{\omega}{s}}}
#'
#'  The y-bounded r-th raw moment of the distribution equals:
#'
#'  \deqn{s^{\sigma_s - 1} \Gamma\left(\sigma_s +1, \frac{y}{s} \right)}
#'
#'  where \eqn{\Gamma(,)} denotes the upper incomplete gamma function.
#'
#' @return Returns the truncated rth raw moment of the distribution.
#'
#' @examples
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' pexp(2, rate = 1)
#' mexp(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rexp(1e5, rate = 1)
#' mean(x)
#' mexp(r = 1, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mexp(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name exp

NULL

#' @rdname exp
#' @export
#'

mexp <- function(r = 0, truncation = 0, rate = 1, lower.tail = TRUE) {

  # Recover code from Weibull (exponential is weibull with shape=1)
  m <- mweibull(r = r, truncation = truncation, shape = 1, scale = 1 / rate, lower.tail = lower.tail)

  return(m)
}
