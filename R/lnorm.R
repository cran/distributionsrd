#' The Lognormal distribution
#'
#' Raw moments for the Lognormal distribution.
#' @param truncation lower truncation parameter, defaults to 0.
#' @param r rth raw moment of the distribution, defaults to 1.
#' @param meanlog,sdlog mean and standard deviation of the distribution on the log scale with default values of 0 and 1 respectively.
#' @param lower.tail logical; if TRUE (default), moments are \eqn{E[x^r|X \le y]}, otherwise, \eqn{E[x^r|X > y]}
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) = \frac{1}{{x Var \sqrt {2\pi } }}e^{- (lnx - \mu  )^2/ 2Var^2} , \qquad F_X(x) = \Phi(\frac{lnx- \mu}{Var})}
#'
#'  The y-bounded r-th raw moment of the Lognormal distribution equals:
#'
#'  \deqn{\mu^r_y =  e^{\frac{r (rVar^2 + 2\mu)}{2}}[1-\Phi(\frac{lny - (rVar^2 + \mu)}{Var})]  }
#'
#' @return Provides the y-bounded, rth raw moment of the distribution.
#'
#' @examples
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' plnorm(2, meanlog = -0.5, sdlog = 0.5)
#' mlnorm(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rlnorm(1e5, meanlog = -0.5, sdlog = 0.5)
#' mean(x)
#' mlnorm(r = 1, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mlnorm(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name lnorm

NULL

#' @rdname lnorm
#' @export

mlnorm <- function(r = 0, truncation = 0, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  m <- exp(1 / 2 * (sdlog * r)^2 + meanlog * r) * (pnorm((log(truncation) - (sdlog^2) * r - meanlog) / sdlog))
  notrunc <- exp(1 / 2 * (sdlog * r)^2 + meanlog * r) * (pnorm((log(0) - meanlog) / sdlog - sdlog * r, lower.tail = FALSE))

  if (!lower.tail) {
    m <- notrunc - m
  }

  return(m)
}

#' Log Normal coefficients of power-law transformed log normal
#'
#' Coefficients of a power-law transformed log normal distribution
#' @param meanlog,sdlog mean and standard deviation of the log normal distributed variable, defaults to 0 and 1 respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable y is log normally distributed with mean meanlog and standard deviation sdlog, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is log normally distributed with mean \eqn{ \frac{meanlog - ln(a)}{b} } and standard deviation \eqn{\frac{sdlog}{b}}.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  plnorm(3,meanlog=-0.5,sdlog=0.5)
#'  coeff = lnorm_plt(meanlog=-0.5,sdlog=0.5,a=5,b=7)$coefficients
#'  plnorm(5*3^7,meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]])
#'
#'  plnorm(5*0.8^7,meanlog=-0.5,sdlog=0.5)
#'  coeff = lnorm_plt(meanlog=-0.5,sdlog=0.5,a=5,b=7,inv=TRUE)$coefficients
#'  plnorm(0.8,meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]])
#'
#'  ## Comparing the first moments and sample means of power-law transformed variables for large enough samples
#'  x = rlnorm(1e5,meanlog=-0.5,sdlog=0.5)
#'  coeff = lnorm_plt(meanlog=-0.5,sdlog=0.5,a=2,b=0.5)$coefficients
#'  y = rlnorm(1e5,meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]])
#'  mean(2*x^0.5)
#'  mean(y)
#'  mlnorm(r=1,meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]],lower.tail=FALSE)
#'
#' @export

lnorm_plt <- function(meanlog = 0, sdlog = 1, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  newmeanlog <- (meanlog - log(a)) / b
  newsdlog <- sdlog / b

  return_list <- list(coefficients = c(meanlog = newmeanlog, sdlog = newsdlog))
  return(return_list)
}
