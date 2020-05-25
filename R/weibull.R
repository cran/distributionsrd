#' The Weibull distribution
#'
#' Raw moments for the Weibull distribution.
#' @param truncation lower truncation parameter, defaults to 0.
#' @param r rth raw moment of the distribution, defaults to 1.
#' @param shape,scale shape and scale of the distribution with default values of 2 and 1 respectively.
#' @param lower.tail logical; if TRUE (default), moments are \eqn{E[x^r|X \le y]}, otherwise, \eqn{E[x^r|X > y]}
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) = \frac{shape}{scale}(\frac{\omega}{scale})^{shape-1}e^{-(\frac{\omega}{scale})^shape} , \qquad F_X(x) = 1-e^{-(\frac{\omega}{scale})^shape}}
#'
#'  The y-bounded r-th raw moment of the distribution equals:
#'
#'  \deqn{\mu^r_y =   scale^{r} \Gamma(\frac{r}{shape} +1, (\frac{y}{scale})^shape ) }
#'
#'  where \eqn{\Gamma(,)} denotes the upper incomplete gamma function.
#'
#' @return returns the truncated rth raw moment of the distribution.
#'
#' @examples
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' pweibull(2, shape = 2, scale = 1)
#' mweibull(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rweibull(1e5, shape = 2, scale = 1)
#' mean(x)
#' mweibull(r = 1, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mweibull(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name weibull

NULL

#' @rdname weibull
#' @export
#'

mweibull <- function(r = 0, truncation = 0, shape = 2, scale = 1, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  m <- (scale^(r)) * pgamma(((truncation / scale)^shape), ((r) / shape + 1), lower.tail = FALSE) * gamma(((r) / shape + 1))
  notrunc <- (scale^(r)) * pgamma(((0 / scale)^shape), ((r) / shape + 1), lower.tail = FALSE) * gamma(((r) / shape + 1))

  if (lower.tail) {
    m <- notrunc - m
  }

  return(m)
}

#' Weibull coefficients of power-law transformed Weibull
#'
#' Coefficients of a power-law transformed Weibull distribution
#' @param shape,scale shape and scale of the distribution with default values of 2 and 1 respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable y is Weibull distributed with mean meanlog and standard deviation sdlog, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Weibull distributed with scale \eqn{ ( \frac{scale}{a})^{\frac{1}{b}} } and shape \eqn{b*shape}.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  pweibull(3,shape=2,scale=1)
#'  coeff = weibull_plt(shape=2,scale=1,a=5,b=7)$coefficients
#'  pweibull(5*3^7,shape=coeff[["shape"]],scale=coeff[["scale"]])
#'
#'  pweibull(5*0.8^7,shape=2,scale=1)
#'  coeff = weibull_plt(shape=2,scale=1,a=5,b=7,inv=TRUE)$coefficients
#'  pweibull(0.8,shape=coeff[["shape"]],scale=coeff[["scale"]])
#'
#'  ## Comparing the first moments and sample means of power-law transformed variables for large enough samples
#'  x = rweibull(1e5,shape=2,scale=1)
#'  coeff = weibull_plt(shape=2,scale=1,a=2,b=0.5)$coefficients
#'  y = rweibull(1e5,shape=coeff[["shape"]],scale=coeff[["scale"]])
#'  mean(2*x^0.5)
#'  mean(y)
#'  mweibull(r=1,shape=coeff[["shape"]],scale=coeff[["scale"]],lower.tail=FALSE)
#'
#' @export

weibull_plt <- function(scale = 1, shape = 2, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  newshape <- b * shape
  newscale <- (scale / a)^(1 / b)

  return_list <- list(coefficients = c(scale = newscale, shape = newshape))
  return(return_list)
}
