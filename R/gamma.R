#' The Gamma distribution
#'
#' Raw moments for the Gamma distribution.
#' @param truncation lower truncation parameter, defaults to 0.
#' @param r rth raw moment of the distribution, defaults to 1.
#' @param shape,rate,scale shape, rate and scale of the distribution with default values of 2 and 1 respectively.
#' @param lower.tail logical; if TRUE (default), moments are \eqn{E[x^r|X \le y]}, otherwise, \eqn{E[x^r|X > y]}
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) = \frac{1}{s^k\Gamma(k)}\omega^{k-1}e^{-\frac{\omega}{s}},\qquad F_X(x) = \frac{1}{\Gamma(k)}\gamma(k,\frac{\omega}{s})},
#'
#'  where \eqn{\Gamma(x)} stands for the upper incomplete gamma function function, while \eqn{\gamma(s,x)} stands for the lower incomplete Gamma function with upper bound \eqn{x}.
#'
#'  The y-bounded r-th raw moment of the distribution equals:
#'
#'  \deqn{\mu^r_y =   \frac{s^{r}}{\Gamma(k)} \Gamma\left(r + k , \frac{y}{s} \right)  }
#'
#' @return Provides the truncated rth raw moment of the distribution.
#'
#'  ## The zeroth truncated moment is equivalent to the probability function
#'  pgamma(2,shape=2,rate=1)
#'  mgamma(truncation=2)
#'
#'  ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#'  #for large enough samples.
#'  x = rgamma(1e5,shape=2,rate=1)
#'  mean(x)
#'  mgamma(r=1,lower.tail=FALSE)
#'
#'  sum(x[x>quantile(x,0.1)])/length(x)
#'  mgamma(r=1,truncation=quantile(x,0.1),lower.tail=FALSE)
#'
#' @name gamma

NULL

#' @rdname gamma
#' @export
#'

mgamma <- function(r = 0, truncation = 0, shape = 2, rate = 1, scale = 1 / rate, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  m <- (scale^(r)) / gamma(shape / 1) * pgamma(((truncation / scale)^1), ((r + shape) / 1), lower.tail = FALSE) * gamma(((r + shape) / 1))
  notrunc <- (scale^(r)) / gamma(shape / 1) * pgamma(((0 / scale)^1), ((r + shape) / 1), lower.tail = FALSE) * gamma(((r + shape) / 1))

  if (lower.tail) {
    m <- notrunc - m
  }

  return(m)
}
