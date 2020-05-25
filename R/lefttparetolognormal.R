#' The Left-Pareto Lognormal distribution
#'
#' Density, distribution function, quantile function and random generation for the Left-Pareto Lognormal distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param shape1,meanlog,sdlog Shape, mean and variance of the Left-Pareto Lognormal distribution respectively.
#' @param r rth raw moment of the Pareto distribution
#' @param truncation lower truncation parameter, defaults to xmin
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{(E[x^r|X \le y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#'
#' @details Probability and Cumulative Distribution Function as provided by \insertCite{reed2004double}{distributionsrd}:
#'
#'  \eqn{f(x) = shape1\omega^{shape1-1}e^{-shape1 meanlog + \frac{shape1^2sdlog^2}{2}}\Phi^c(\frac{ln\omega - meanlog +shape1 sdlog^2}{sdlog}), \newline F_X(x) = \Phi(\frac{ln\omega - meanlog }{sdlog}) -  		\omega^{shape1}e^{-shape1 meanlog + \frac{shape1^2sdlog^2}{2}}\Phi^c(\frac{ln\omega - meanlog +shape1 sdlog^2}{sdlog})}
#'
#'  The y-bounded r-th raw moment of the Let-Pareto Lognormal distribution equals:
#'
#'  \eqn{meanlog^{r}_{y} = -shape1e^{-shape1 meanlog + \frac{shape1^2sdlog^2}{2}}\frac{y^{\sigma_s + shape1-1}}{\sigma_s + shape1 - 1}\Phi^c(\frac{lny - meanlog + shape1 sdlog^2}{sdlog}) \newline \qquad + \frac{shape1}{r+shape1} e^{\frac{ r^2sdlog^2 + 2meanlog r }{2}}\Phi^c(\frac{lny -  rsdlog^2 + meanlog}{sdlog})}
#'
#' @return dleftparetolognormal gives the density, pleftparetolognormal gives the distribution function, qleftparetolognormal gives the quantile function, mleftparetolognormal gives the rth moment of the distribution and rleftparetolognormal generates random deviates.
#'
#'  The length of the result is determined by n for rleftparetolognormal, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references
#' \insertAllCited{}
#'
#'  ## Left-Pareto Lognormal density
#'  plot(x=seq(0,5,length.out=100),y=dleftparetolognormal(x=seq(0,5,length.out=100)))
#'  plot(x=seq(0,5,length.out=100),y=dleftparetolognormal(x=seq(0,5,length.out=100),shape1=1))
#'
#'  ## Left-Pareto Lognormal relates to the Lognormal if the shape parameter goes to infinity
#'  pleftparetolognormal(q=6,shape1=1e20,meanlog=-0.5,sdlog=0.5)
#'  plnorm(q=6,meanlog=-0.5,sdlog=0.5)
#'
#'  ## Demonstration of log functionality for probability and quantile function
#'  qleftparetolognormal(pleftparetolognormal(2,log.p=TRUE),log.p=TRUE)
#'
#'  ## The zeroth truncated moment is equivalent to the probability function
#'  pleftparetolognormal(2)
#'  mleftparetolognormal(truncation=2)
#'
#'  ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#'  #for large enough samples.
#'  x = rleftparetolognormal(1e5)
#'
#'  mean(x)
#'  mleftparetolognormal(r=1,lower.tail=FALSE)
#'
#'  sum(x[x>quantile(x,0.1)])/length(x)
#'  mleftparetolognormal(r=1,truncation=quantile(x,0.1),lower.tail=FALSE)
#'
#' @name leftparetolognormal

NULL

#' @rdname leftparetolognormal
#' @export

dleftparetolognormal <- function(x, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, log = FALSE) {
  d <- shape1 * x^(shape1 - 1) * exp((-shape1 * meanlog + ((-shape1)^2 * sdlog^2) / 2)) * pnorm(((log(x) - meanlog + shape1 * sdlog^2) / sdlog), lower.tail = FALSE)
  if (log) {
    d <- log(d)
  }

  return(d)
}

#' @rdname leftparetolognormal
#' @export

pleftparetolognormal <- function(q, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE, log.p = FALSE) {
  P1 <- pnorm(((log(q) - meanlog) / sdlog))
  P2 <- q^(shape1) * exp(-shape1 * meanlog + ((-shape1)^2 * sdlog^2) / 2) * pnorm(((log(q) - meanlog + shape1 * sdlog^2) / sdlog), lower.tail = FALSE)

  P2[pnorm(((log(q) - meanlog + shape1 * sdlog^2) / sdlog), lower.tail = FALSE) == 0] <- 0

  p <- P1 + P2

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

#' @rdname leftparetolognormal
#' @export

qleftparetolognormal <- function(p, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE, log.p = FALSE) {
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


  q <- suppressWarnings(uniroot(function(q, shape1, meanlog, sdlog) {
    pleftparetolognormal(q, shape1, meanlog, sdlog) - p
  }, interval = c(1e-6, 1), extendInt = c("yes"), shape1, meanlog, sdlog)$root)
  return(q)
}
qleftparetolognormal <- Vectorize(qleftparetolognormal, vectorize.args = "p")

#' @rdname leftparetolognormal
#' @export

mleftparetolognormal <- function(r = 0, truncation = 0, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  I1 <- -shape1 * exp(-shape1 * meanlog + 1 / 2 * (shape1^2 * sdlog^2)) * (truncation^(r + shape1)) / (r + shape1) * (1 - pnorm(((log(truncation) - meanlog + shape1 * sdlog^2) / sdlog)))
  I2 <- (shape1 * (exp((r) * meanlog + ((r) * sdlog)^2 / 2)) / (r + shape1) * (1 - pnorm((log(truncation) - (meanlog - shape1 * sdlog^2)) / sdlog - sdlog * ((r) + shape1))))
  I1[is.na(I1)] <- 0
  moment <- I1 + I2

  notrunc_I1 <- -shape1 * exp(-shape1 * meanlog + 1 / 2 * (shape1^2 * sdlog^2)) * (0^(r + shape1)) / (r + shape1) * (1 - pnorm(((log(0) - meanlog + shape1 * sdlog^2) / sdlog)))
  notrunc_I2 <- (shape1 * (exp((r) * meanlog + ((r) * sdlog)^2 / 2)) / (r + shape1) * (1 - pnorm((log(0) - (meanlog - shape1 * sdlog^2)) / sdlog - sdlog * ((r) + shape1))))
  notrunc_I1[is.na(notrunc_I1)] <- 0
  notrunc_moment <- notrunc_I1 + notrunc_I2

  if (lower.tail) {
    moment <- notrunc_moment - moment
  }

  return(moment)
}


#' @rdname leftparetolognormal
#' @export

rleftparetolognormal <- function(n, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5) {

  # r = qleftparetolognormal(runif(n),shape1=shape1,meanlog=meanlog,sdlog=sdlog)
  # See equation 1.2.16 in master thesis Chuan Chuan Zhang
  r <- exp(meanlog + sdlog * rnorm(n) - rexp(n) / shape1)
  return(r)
}

#' Left-Pareto Lognormal coefficients of power-law transformed Left-Pareto Lognormal
#'
#' Coefficients of a power-law transformed Left-Pareto Lognormal distribution
#' @param shape1,meanlog,sdlog Shapes, mean and variance of the Left-Pareto Lognormal distribution respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable y is Left-Pareto Lognormal distributed with mean meanlog and standard deviation sdlog, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Left-Pareto Lognormal distributed with \eqn{ shape1*b, \frac{meanlog-log(a)}{b}, \frac{sdlog}{b} }.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#' @examples
#'
#' ## Comparing probabilites of power-law transformed transformed variables
#' pleftparetolognormal(3, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5)
#' coeff <- leftparetolognormal_plt(shape1 = 1.5, meanlog = -0.5, sdlog = 0.5,
#' a = 5, b = 7)$coefficients
#' pleftparetolognormal(5 * 3^7, shape1 = coeff[["shape1"]], meanlog = coeff[["meanlog"]],
#' sdlog = coeff[["sdlog"]])
#'
#' pleftparetolognormal(5 * 0.9^7, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5)
#' coeff <- leftparetolognormal_plt(shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, a = 5, b = 7,
#' inv = TRUE)$coefficients
#' pleftparetolognormal(0.9, shape1 = coeff[["shape1"]], meanlog = coeff[["meanlog"]],
#'  sdlog = coeff[["sdlog"]])
#' @export

leftparetolognormal_plt <- function(shape1 = 1.5, meanlog = -0.5, sdlog = 0.5, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  # newmeanlog = (b*meanlog + log(a))
  # newshape1 = shape1/b

  newmeanlog <- (meanlog - log(a)) / b
  newshape1 <- shape1 * b
  newsdlog <- sdlog / b

  return_list <- list(coefficients = c(shape1 = newshape1, meanlog = newmeanlog, sdlog = newsdlog))
  return(return_list)
}

#' Left-Pareto Lognormal MLE
#'
#' Maximum likelihood estimation of the parameters of the Left-Pareto Lognormal distribution.
#' @param x data vector
#' @param lower,upper Upper and lower bounds for the estimation procedure on the parameters c(shape1,sdlog), defaults to c(1e-10,1e-10) and c(Inf,Inf) respectively.
#' @param start named vector with starting values, default to c(shape1=2,sdlog=sd(log(x)))
#'
#' @return Returns a named list containing a
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' \item{convergence}{logical indicator of convergence}
#' \item{n}{Length of the fitted data vector}
#' \item{np}{Nr. of coefficients}
#' }
#'
#' x = rleftparetolognormal(1e3)
#'
#' ## Pareto fit with xmin set to the minium of the sample
#' leftparetolognormal.mle(x=x)
#'
#' @export

leftparetolognormal.mle <- function(x, lower = c(1e-10, 1e-10), upper = c(Inf, Inf), start = NULL) {
  y <- log(x)
  n <- length(y)

  logl <- function(x, shape1, meanlog, sdlog, n) {
    n * log(shape1) + sum(log(dnorm((y - meanlog) / sdlog))) +
      sum(log(pnorm((shape1 * sdlog + (y - meanlog) / sdlog), lower.tail = FALSE) / dnorm((shape1 * sdlog + (y - meanlog) / sdlog))))
  }

  mle <- function(par, y, n) {
    shape1 <- par[1]
    meanlog <- mean(y) + 1 / shape1
    sdlog <- par[2]

    nlogl <- -logl(x = y, shape1 = shape1, meanlog = meanlog, sdlog = sdlog, n = n)
    return(nlogl)
  }

  if (is.null(start)) {
    start <- c(shape1 = 2, sdlog = sd(log(x)))
  }

  start <- start[c("shape1", "sdlog")]

  optim.out <- suppressWarnings(nlminb(as.numeric(start), mle, y = y, n = n, lower = lower, upper = upper, control = list(maxit = 1e5)))

  par <- optim.out$par

  shape1 <- par[1]
  meanlog <- mean(y) + 1 / shape1
  sdlog <- par[2]

  convergence <- ifelse(optim.out$convergence == 0, TRUE, FALSE)


  return_list <- list(coefficients = c(meanlog = meanlog, sdlog = sdlog, shape1 = shape1), convergence = convergence, n = n, np = 3)
  return(return_list)
}
