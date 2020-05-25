#' The Right-Pareto Lognormal distribution
#'
#' Density, distribution function, quantile function and random generation for the Right-Pareto Lognormal distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param shape2,meanlog,sdlog Shape, mean and variance of the Right-Pareto Lognormal distribution respectively.
#' @param r rth raw moment of the Pareto distribution
#' @param truncation lower truncation parameter, defaults to xmin
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{(E[x^r|X \le y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#'
#' @details Probability and Cumulative Distribution Function as provided by \insertCite{reed2004double}{distributionsrd}:
#'
#'  \eqn{f(x) = shape2 \omega^{-shape2-1}e^{shape2 meanlog + \frac{shape2^2sdlog^2}{2}}\Phi(\frac{lnx - meanlog - shape2 sdlog^2}{sdlog}), \newline F_X(x) = \Phi(\frac{lnx - meanlog }{sdlog}) -  \omega^{-shape2}e^{shape2 meanlog + \frac{shape2^2sdlog^2}{2}}\Phi(\frac{lnx - meanlog - shape2 sdlog^2}{sdlog})}
#'
#'  The y-bounded r-th raw moment of the Right-Pareto Lognormal distribution equals:
#'
#'  \eqn{meanlog^{r}_{y} = -shape2e^{shape2 meanlog + \frac{shape2^2sdlog^2}{2}}\frac{y^{\sigma_s - shape2-1}}{\sigma_s - shape2 - 1}\Phi(\frac{lny - meanlog - shape2 sdlog^2}{sdlog}) \newline \qquad - \frac{shape2}{r-shape2} e^{\frac{ r^2sdlog^2 + 2meanlog r }{2}}\Phi^c(\frac{lny -  rsdlog^2 + meanlog}{sdlog}),  \qquad shape2>r}
#'
#' @return drightparetolognormal gives the density, prightparetolognormal gives the distribution function, qrightparetolognormal gives the quantile function, mrightparetolognormal gives the rth moment of the distribution and rrightparetolognormal generates random deviates.
#'
#'  The length of the result is determined by n for rrightparetolognormal, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#'
#' ## Right-Pareto Lognormal density
#' plot(x = seq(0, 5, length.out = 100), y = drightparetolognormal(x = seq(0, 5, length.out = 100)))
#' plot(x = seq(0, 5, length.out = 100), y = drightparetolognormal(x = seq(0, 5, length.out = 100),
#' shape2 = 1))
#'
#' ## Right-Pareto Lognormal relates to the Lognormal if the shape parameter goes to infinity
#' prightparetolognormal(q = 6, shape2 = 1e20, meanlog = -0.5, sdlog = 0.5)
#' plnorm(q = 6, meanlog = -0.5, sdlog = 0.5)
#'
#' ## Demonstration of log functionality for probability and quantile function
#' qrightparetolognormal(prightparetolognormal(2, log.p = TRUE), log.p = TRUE)
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' prightparetolognormal(2)
#' mrightparetolognormal(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rrightparetolognormal(1e5, shape2 = 3)
#'
#' mean(x)
#' mrightparetolognormal(r = 1, shape2 = 3, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mrightparetolognormal(r = 1, shape2 = 3, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name rightparetolognormal

NULL

#' @rdname rightparetolognormal
#' @export

drightparetolognormal <- function(x, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, log = FALSE) {
  d <- shape2 * x^(-shape2 - 1) * exp((shape2 * meanlog + (shape2^2 * sdlog^2) / 2)) * pnorm(((log(x) - meanlog - shape2 * sdlog^2) / sdlog))
  if (log) {
    d <- log(d)
  }

  return(d)
}

#' @rdname rightparetolognormal
#' @export

prightparetolognormal <- function(q, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE, log.p = FALSE) {
  P1 <- pnorm(((log(q) - meanlog) / sdlog))
  P2 <- q^(-shape2) * exp(shape2 * meanlog + (shape2^2 * sdlog^2) / 2) * pnorm(((log(q) - meanlog - shape2 * sdlog^2) / sdlog))

  P2[pnorm(((log(q) - meanlog - shape2 * sdlog^2) / sdlog)) == 0] <- 0
  p <- P1 - P2

  p <- if (lower.tail) {
    p
  } else {
    1 - p
  }
  p <- if (log.p) {
    log(p)
  } else {
    p
  }

  return(p)
}

#' @rdname rightparetolognormal
#' @export

qrightparetolognormal <- function(p, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE, log.p = FALSE) {
  p <- if (log.p) {
    exp(p)
  } else {
    p
  }
  p <- if (lower.tail) {
    p
  } else {
    1 - p
  }


  q <- suppressWarnings(uniroot(function(q, shape2, meanlog, sdlog) {
    prightparetolognormal(q, shape2, meanlog, sdlog) - p
  }, interval = c(1e-6, 1), extendInt = c("yes"), shape2, meanlog, sdlog)$root)

  return(q)
}
qrightparetolognormal <- Vectorize(qrightparetolognormal, vectorize.args = "p")

#' @rdname rightparetolognormal
#' @export

mrightparetolognormal <- function(r = 0, truncation = 0, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  if (shape2 > r) {
    I1 <- -((shape2) / (1)) * exp(shape2 * meanlog + 1 / 2 * (shape2^2 * sdlog^2)) * (truncation^(r - shape2)) / (r - shape2) * pnorm(((log(truncation) - meanlog - shape2 * sdlog^2) / sdlog))
    I2 <- -(((shape2) / (1)) * (exp(r * meanlog + (r * sdlog)^2 / 2)) / (r - shape2) * (1 - pnorm((log(truncation) - (meanlog + shape2 * sdlog^2)) / sdlog - sdlog * (r - shape2))))
    I1[is.na(I1)] <- 0
    moment <- I1 + I2

    notrunc_I1 <- -((shape2) / (1)) * exp(shape2 * meanlog + 1 / 2 * (shape2^2 * sdlog^2)) * (0^(r - shape2)) / (r - shape2) * pnorm(((log(0) - meanlog - shape2 * sdlog^2) / sdlog))
    notrunc_I2 <- -(((shape2) / (1)) * (exp(r * meanlog + (r * sdlog)^2 / 2)) / (r - shape2) * (1 - pnorm((log(0) - (meanlog + shape2 * sdlog^2)) / sdlog - sdlog * (r - shape2))))
    notrunc_I1[is.na(notrunc_I1)] <- 0
    notrunc_moment <- notrunc_I1 + notrunc_I2

    if (lower.tail) {
      moment <- notrunc_moment - moment
    }
  } else {
    moment <- rep(NA, length(truncation))
  }
  return(moment)
}

#' @rdname rightparetolognormal
#' @export

rrightparetolognormal <- function(n, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE) {

  # r = qrightparetolognormal(runif(n),shape2=shape2,meanlog=meanlog,sdlog=sdlog)
  # See equation 1.2.16 in master thesis Chuan Chuan Zhang
  r <- exp(meanlog + sdlog * rnorm(n) + rexp(n) / shape2)
  return(r)
}

#' Right-Pareto Lognormal coefficients of power-law transformed Right-Pareto Lognormal
#'
#' Coefficients of a power-law transformed Right-Pareto Lognormal distribution
#' @param shape2,meanlog,sdlog Shapes, mean and variance of the Right-Pareto Lognormal distribution respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable y is Right-Pareto Lognormal distributed with mean meanlog and standard deviation sdlog, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Right-Pareto Lognormal distributed with \eqn{ \frac{meanlog-log(a)}{b}, \frac{sdlog}{b}, shape2*b }.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  prightparetolognormal(3, shape2 = 3, meanlog = -0.5, sdlog = 0.5)
#'  coeff = rightparetolognormal_plt( shape2 = 3, meanlog = -0.5, sdlog = 0.5,a=5,b=7)$coefficients
#'  prightparetolognormal(5*3^7,shape2=coeff[["shape2"]],meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]])
#'
#'  prightparetolognormal(5*0.9^7,shape2 = 3, meanlog = -0.5, sdlog = 0.5)
#'  coeff = rightparetolognormal_plt(shape2 = 3, meanlog = -0.5, sdlog = 0.5,a=5,b=7, inv=TRUE)$coefficients
#'  prightparetolognormal(0.9,shape2=coeff[["shape2"]],meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]])
#'
#'
#' @export

rightparetolognormal_plt <- function(shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  # newmeanlog = (b*meanlog + log(a))
  # newshape2 = shape2/b

  newmeanlog <- (meanlog - log(a)) / b
  newshape2 <- shape2 * b
  newsdlog <- sdlog / b

  return_list <- list(coefficients = c(shape2 = newshape2, meanlog = newmeanlog, sdlog = newsdlog))
  return(return_list)
}

#' Right-Pareto Lognormal MLE
#'
#' Maximum likelihood estimation of the parameters of the Right-Pareto Lognormal distribution.
#' @param x data vector
#' @param lower,upper Upper and lower bounds for the estimation procedure on the parameters c(shape2,sdlog), defaults to c(1e-10,1e-10) and c(Inf,Inf) respectively.
#' @param start named vector with starting values, default to c(shape2=2,sdlog=sd(log(x)))
#'
#' @return Returns a named list containing a
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' \item{convergence}{logical indicator of convergence}
#' \item{n}{Length of the fitted data vector}
#' \item{np}{Nr. of coefficients}
#' }
#'
#' @examples
#'
#' x <- rrightparetolognormal(1e3)
#'
#' ## Pareto fit with xmin set to the minium of the sample
#' rightparetolognormal.mle(x = x)
#' @export

rightparetolognormal.mle <- function(x, lower = c(1e-10, 1e-10), upper = c(Inf, Inf), start = NULL) {
  y <- log(x)
  n <- length(y)

  logl <- function(x, shape2, meanlog, sdlog, n) {
    n * log(shape2) + sum(log(dnorm((y - meanlog) / sdlog))) +
      sum(log((pnorm((shape2 * sdlog - (y - meanlog) / sdlog), lower.tail = FALSE) / dnorm((shape2 * sdlog - (y - meanlog) / sdlog)))))
  }

  mle <- function(par, y, n) {
    shape2 <- par[1]
    meanlog <- mean(y) - 1 / shape2
    sdlog <- par[2]

    nlogl <- -logl(x = y, shape2 = shape2, meanlog = meanlog, sdlog = sdlog, n = n)
    return(nlogl)
  }

  if (is.null(start)) {
    start <- c(shape2 = 2, sdlog = sd(log(x)))
  }

  start <- start[c("shape2", "sdlog")]

  optim.out <- suppressWarnings(nlminb(as.numeric(start), mle, y = y, n = n, lower = lower, upper = upper, control = list(maxit = 1e5)))

  par <- optim.out$par

  shape2 <- par[1]
  meanlog <- mean(y) - 1 / shape2
  sdlog <- par[2]

  convergence <- ifelse(optim.out$convergence == 0, TRUE, FALSE)

  return_list <- list(coefficients = c(meanlog = meanlog, sdlog = sdlog, shape2 = shape2), convergence = convergence, n = n, np = 3)
  return(return_list)
}
