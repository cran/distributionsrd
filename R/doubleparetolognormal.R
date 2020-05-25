#' The Double-Pareto Lognormal distribution
#'
#' Density, distribution function, quantile function and random generation for the Double-Pareto Lognormal distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param shape1,shape2,meanlog,sdlog Shapes, mean and variance of the Double-Pareto Lognormal distribution respectively, defaults to shape1=1.5, shape2=1.5, meanlog=-0.5, sdlog=0.5.
#' @param r rth raw moment of the Pareto distribution
#' @param truncation lower truncation parameter, defaults to xmin
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{(E[x^r|X \le y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#'
#' @details Probability and Cumulative Distribution Function as provided by \insertCite{reed2004double}{distributionsrd}:
#'
#'  \eqn{f(x) = \frac{shape2shape1}{shape2 + shape1}[ x^{-shape2-1}e^{shape2 meanlog + \frac{shape2^2sdlog^2}{2}}\Phi(\frac{lnx - meanlog - shape2 sdlog^2}{sdlog}) +  x^{shape1-1}e^{-shape1 meanlog + \frac{shape1^2sdlog^2}{2}}\Phi^c(\frac{lnx - meanlog +shape1 sdlog^2}{sdlog}) ], \newline F_X(x) = \Phi(\frac{lnx - meanlog }{sdlog}) - \frac{1}{shape2 + shape1}[  shape1 x^{-shape2}e^{shape2 meanlog + \frac{shape2^2sdlog^2}{2}}\Phi(\frac{lnx - meanlog - shape2 sdlog^2}{sdlog}) - shape2 x^{shape1}e^{-shape1 meanlog + \frac{shape1^2sdlog^2}{2}}\Phi^c(\frac{lnx - meanlog +shape1 sdlog^2}{sdlog})  ]}
#'
#'  The y-bounded r-th raw moment of the Double-Pareto Lognormal distribution equals:
#'
#'  \eqn{meanlog^{r}_{y} = -\frac{shape2shape1}{shape2 + shape1}e^{shape2 meanlog + \frac{shape2^2sdlog^2}{2}}\frac{y^{r- shape2}}{r - shape2}\Phi(\frac{lny - meanlog - shape2 sdlog^2}{sdlog}) \newline \qquad- \frac{shape2shape1}{shape2 + shape1} \frac{1}{r-shape2} e^{\frac{ r^2sdlog^2 + 2meanlog r }{2}}\Phi^c(\frac{lny -  rsdlog^2 - meanlog}{sdlog}) \newline  \qquad -\frac{shape2shape1}{shape2 + shape1}e^{-shape1 meanlog + \frac{shape1^2sdlog^2}{2}}\frac{y^{r + shape1}}{r + shape1}\Phi^c(\frac{lny - meanlog + shape1 sdlog^2}{sdlog}) \newline  \qquad + \frac{shape2shape1}{shape2 + shape1} \frac{1}{r+shape1} e^{\frac{ r^2sdlog^2 + 2meanlog r }{2}}\Phi^c(\frac{lny -  rsdlog^2 - meanlog}{sdlog}),  \qquad shape2>r}
#'
#' @return ddoubleparetolognormal returns the density, pdoubleparetolognormal the distribution function, qdoubleparetolognormal the quantile function, mdoubleparetolognormal the rth moment of the distribution and rdoubleparetolognormal generates random deviates.
#'
#'  The length of the result is determined by n for rdoubleparetolognormal, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#'
#' ## Double-Pareto Lognormal density
#' plot(x = seq(0, 5, length.out = 100), y = ddoubleparetolognormal(x = seq(0, 5, length.out = 100)))
#' plot(x = seq(0, 5, length.out = 100), y = ddoubleparetolognormal(x = seq(0, 5, length.out = 100),
#'  shape2 = 1))
#'
#' ## Double-Pareto Lognormal relates to the right-pareto Lognormal distribution if
#' #shape1 goes to infinity
#' pdoubleparetolognormal(q = 6, shape1 = 1e20, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5)
#' prightparetolognormal(q = 6, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5)
#'
#' ## Double-Pareto Lognormal relates to the left-pareto Lognormal distribution if
#' # shape2 goes to infinity
#' pdoubleparetolognormal(q = 6, shape1 = 1.5, shape2 = 1e20, meanlog = -0.5, sdlog = 0.5)
#' pleftparetolognormal(q = 6, shape1 = 1.5, meanlog = -0.5, sdlog = 0.5)
#'
#' ## Double-Pareto Lognormal relates to the Lognormal if both shape parameters go to infinity
#' pdoubleparetolognormal(q = 6, shape1 = 1e20, shape2 = 1e20, meanlog = -0.5, sdlog = 0.5)
#' plnorm(q = 6, meanlog = -0.5, sdlog = 0.5)
#'
#' ## Demonstration of log functionality for probability and quantile function
#' qdoubleparetolognormal(pdoubleparetolognormal(2, log.p = TRUE), log.p = TRUE)
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' pdoubleparetolognormal(2)
#' mdoubleparetolognormal(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rdoubleparetolognormal(1e5, shape2 = 3)
#'
#' mean(x)
#' mdoubleparetolognormal(r = 1, shape2 = 3, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mdoubleparetolognormal(r = 1, shape2 = 3, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name doubleparetolognormal

NULL

#' @rdname doubleparetolognormal
#' @export

ddoubleparetolognormal <- function(x, shape1 = 1.5, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, log = FALSE) {
  d <- (shape2 * shape1) / (shape2 + shape1) * (
    x^(-shape2 - 1) * exp((shape2 * meanlog + (shape2^2 * sdlog^2) / 2)) * pnorm(((log(x) - meanlog - shape2 * sdlog^2) / sdlog)) +
      pnorm(((log(x) - meanlog + shape1 * sdlog^2) / sdlog), lower.tail = FALSE) * x^(shape1 - 1) * exp((-shape1 * meanlog + (shape1^2 * sdlog^2) / 2)))

  if (log) {
    d <- log(shape2) + log(shape1) - log(shape2 + shape1) + log(x^(-shape2 - 1) * exp((shape2 * meanlog + (shape2^2 * sdlog^2) / 2)) * pnorm(((log(x) - meanlog - shape2 * sdlog^2) / sdlog)) +
      pnorm(((log(x) - meanlog + shape1 * sdlog^2) / sdlog), lower.tail = FALSE) * x^(shape1 - 1) * exp((-shape1 * meanlog + (shape1^2 * sdlog^2) / 2)))
  }
  return(d)
}

#' @rdname doubleparetolognormal
#' @export

pdoubleparetolognormal <- function(q, shape1 = 1.5, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE, log.p = FALSE) {
  P1 <- pnorm(((log(q) - meanlog) / sdlog))
  P2 <- shape1 * q^(-shape2) * exp(shape2 * meanlog + (shape2^2 * sdlog^2) / 2) * pnorm(((log(q) - meanlog - shape2 * sdlog^2) / sdlog))
  P2[pnorm(((log(q) - meanlog - shape2 * sdlog^2) / sdlog)) == 0] <- 0
  P3 <- shape2 * q^(shape1) * exp(-shape1 * meanlog + ((-shape1)^2 * sdlog^2) / 2)
  P4 <- pnorm(((log(q) - meanlog + shape1 * sdlog^2) / sdlog), lower.tail = FALSE)
  P3P4 <- P3 * P4
  P3P4[P4 == 0] <- 0

  p <- P1 - (1) / (shape2 + shape1) * (P2 - P3P4)

  # p = pnorm(((log(q)-meanlog)/sdlog)) - (1)/(shape2 + shape1) * (
  #   shape1*q^(-shape2) * exp(shape2*meanlog + (shape2^2 * sdlog^2)/2) * pnorm( ((log(q) - meanlog - shape2*sdlog^2)/sdlog) ) -
  #     shape2*q^(shape1) * exp(-shape1*meanlog + ((-shape1)^2 * sdlog^2)/2) * pnorm( ((log(q) - meanlog + shape1*sdlog^2)/sdlog), lower.tail=FALSE ) )

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

#' @rdname doubleparetolognormal
#' @export

qdoubleparetolognormal <- function(p, shape1 = 1.5, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE, log.p = FALSE) {
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

  q <- suppressWarnings(uniroot(function(q, shape1, shape2, meanlog, sdlog, p) {
    pdoubleparetolognormal(q, shape1, shape2, meanlog, sdlog) - p
  }, interval = c(1e-6, 1), extendInt = c("yes"), shape1, shape2, meanlog, sdlog, p)$root)

  return(q)
}
qdoubleparetolognormal <- Vectorize(qdoubleparetolognormal, vectorize.args = "p")

#' @rdname doubleparetolognormal
#' @export

mdoubleparetolognormal <- function(r = 0, truncation = 0, shape1 = 1.5, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  if (shape2 > r) {
    I1 <- -((shape2 * shape1) / (shape2 + shape1)) * exp(shape2 * meanlog + 1 / 2 * (shape2^2 * sdlog^2)) * (truncation^(r - shape2)) / (r - shape2) * pnorm(((log(truncation) - meanlog - shape2 * sdlog^2) / sdlog))
    I1[is.na(I1)] <- 0
    I2 <- -(((shape2 * shape1) / (shape2 + shape1)) * (exp(r * meanlog + (r * sdlog)^2 / 2)) / (r - shape2) * (1 - pnorm((log(truncation) - (meanlog + shape2 * sdlog^2)) / sdlog - sdlog * (r - shape2))))
    I3 <- -((shape2 * shape1) / (shape2 + shape1)) * exp(-shape1 * meanlog + 1 / 2 * (shape1^2 * sdlog^2)) * (truncation^(r + shape1)) / (r + shape1) * (1 - pnorm(((log(truncation) - meanlog + shape1 * sdlog^2) / sdlog)))
    I3[(1 - pnorm(((log(truncation) - meanlog + shape1 * sdlog^2) / sdlog))) == 0] <- 0
    I4 <- (((shape2 * shape1) / (shape2 + shape1)) * (exp(r * meanlog + (r * sdlog)^2 / 2)) / (r + shape1) * (1 - pnorm((log(truncation) - (meanlog - shape1 * sdlog^2)) / sdlog - sdlog * (r + shape1))))
    I4[is.na(I4)] <- 0
    moment <- I1 + I2 + I3 + I4

    notrunc_I1 <- -((shape2 * shape1) / (shape2 + shape1)) * exp(shape2 * meanlog + 1 / 2 * (shape2^2 * sdlog^2)) * (0^(r - shape2)) / (r - shape2) * pnorm(((log(0) - meanlog - shape2 * sdlog^2) / sdlog))
    notrunc_I1[is.na(notrunc_I1)] <- 0
    notrunc_I2 <- -(((shape2 * shape1) / (shape2 + shape1)) * (exp(r * meanlog + (r * sdlog)^2 / 2)) / (r - shape2) * (1 - pnorm((log(0) - (meanlog + shape2 * sdlog^2)) / sdlog - sdlog * (r - shape2))))
    notrunc_I3 <- -((shape2 * shape1) / (shape2 + shape1)) * exp(-shape1 * meanlog + 1 / 2 * (shape1^2 * sdlog^2)) * (0^(r + shape1)) / (r + shape1) * (1 - pnorm(((log(0) - meanlog + shape1 * sdlog^2) / sdlog)))
    notrunc_I4 <- (((shape2 * shape1) / (shape2 + shape1)) * (exp(r * meanlog + (r * sdlog)^2 / 2)) / (r + shape1) * (1 - pnorm((log(0) - (meanlog - shape1 * sdlog^2)) / sdlog - sdlog * (r + shape1))))
    notrunc_I4[is.na(notrunc_I4)] <- 0
    notrunc_moment <- notrunc_I1 + notrunc_I2 + notrunc_I3 + notrunc_I4

    if (lower.tail) {
      moment <- notrunc_moment - moment
    }
  } else {
    moment <- rep(NA, length(truncation))
  }

  return(moment)
}

#' @rdname doubleparetolognormal
#' @export

rdoubleparetolognormal <- function(n, shape1 = 1.5, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5) {

  # r = qdoubleparetolognormal(runif(n),shape1=shape1,shape2=shape2,meanlog=meanlog,sdlog=sdlog)
  # See equation 1.2.16 in master thesis Chuan Chuan Zhang
  r <- exp(meanlog + sdlog * rnorm(n) + rexp(n) / shape2 - rexp(n) / shape1)
  return(r)
}

#' Double-Pareto Lognormal coefficients of power-law transformed Double-Pareto Lognormal
#'
#' Coefficients of a power-law transformed Double-Pareto Lognormal distribution
#' @param shape1,shape2,meanlog,sdlog Shapes, mean and variance of the Double-Pareto Lognormal distribution respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable y is Double-Pareto Lognormal distributed with mean meanlog and standard deviation sdlog, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Double-Pareto Lognormal distributed with \eqn{ shape1*b, \frac{meanlog-log(a)}{b}, \frac{sdlog}{b}, shape2*b }.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  pdoubleparetolognormal(3,shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5)
#'  coeff = doubleparetolognormal_plt(shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5,a=5,b=7)$coefficients
#'  pdoubleparetolognormal(5*3^7,shape1=coeff[["shape1"]],shape2=coeff[["shape2"]],meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]])
#'
#'  pdoubleparetolognormal(5*0.9^7,shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5)
#'  coeff = doubleparetolognormal_plt(shape1 = 1.5, shape2 = 3, meanlog = -0.5, sdlog = 0.5,a=5,b=7, inv=TRUE)$coefficients
#'  pdoubleparetolognormal(0.9,shape1=coeff[["shape1"]],shape2=coeff[["shape2"]],meanlog=coeff[["meanlog"]],sdlog=coeff[["sdlog"]])
#'
#' @export

doubleparetolognormal_plt <- function(shape1 = 1.5, shape2 = 1.5, meanlog = -0.5, sdlog = 0.5, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  newmeanlog <- (meanlog - log(a)) / b
  newshape1 <- shape1 * b
  newshape2 <- shape2 * b
  newsdlog <- sdlog / b

  # newmeanlog = (b*meanlog + log(a))
  # newshape1 = shape1/b
  # newshape2 = shape2/b

  return_list <- list(coefficients = c(shape1 = newshape1, shape2 = newshape2, meanlog = newmeanlog, sdlog = newsdlog))
  return(return_list)
}

#' Double-Pareto Lognormal MLE
#'
#' Maximum likelihood estimation of the parameters of the Double-Pareto Lognormal distribution.
#' @param x data vector
#' @param lower,upper Upper and lower bounds for the estimation procedure on the parameters c(shape2,shape1,sdlog), defaults to c(1e-10,1e-10,1e-10) and c(Inf,Inf,Inf) respectively.
#' @param start named vector with starting values, default to c(shape2=2,shape1=2,sdlog=sd(log(x)))
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
#' x <- rdoubleparetolognormal(1e3)
#'
#' ## Pareto fit with xmin set to the minium of the sample
#' doubleparetolognormal.mle(x = x)
#' @export

doubleparetolognormal.mle <- function(x, lower = c(1e-10, 1e-10, 1e-10), upper = c(Inf, Inf, Inf), start = NULL) {
  y <- log(x)
  n <- length(y)
  if (is.null(start)) {
    start <- c(shape2 = 2, shape1 = 2, sdlog = sd(log(x)))
  }
  start <- start[c("shape2", "shape1", "sdlog")]

  logl <- function(x, shape2, shape1, meanlog, sdlog, n) {
    n * log(shape2) + n * log(shape1) - n * log(shape2 + shape1) + sum(log(dnorm((y - meanlog) / sdlog))) +
      sum(log((pnorm((shape2 * sdlog - (y - meanlog) / sdlog), lower.tail = FALSE) / dnorm((shape2 * sdlog - (y - meanlog) / sdlog)) +
        pnorm((shape1 * sdlog + (y - meanlog) / sdlog), lower.tail = FALSE) / dnorm((shape1 * sdlog + (y - meanlog) / sdlog)))))
  }

  fnobj <- function(par, y, n) {
    shape2 <- par[1]
    shape1 <- par[2]
    meanlog <- mean(y) - 1 / shape2 + 1 / shape1
    sdlog <- par[3]

    out <- -logl(x = y, shape2 = shape2, shape1 = shape1, meanlog = meanlog, sdlog = sdlog, n = n)

    if (is.infinite(out)) {
      out <- 1e20
    }

    return(out)
  }

  optim.out <- suppressWarnings(nlminb(as.numeric(start), fnobj, y = y, n = n, lower = lower, upper = upper, control = list(maxit = 1e5)))

  par <- optim.out$par

  shape2 <- par[1]
  shape1 <- par[2]
  meanlog <- mean(y) - 1 / shape2 + 1 / shape1
  sdlog <- par[3]

  n <- length(x)

  convergence <- ifelse(optim.out$convergence == 0, TRUE, FALSE)

  return_list <- list(coefficients = c(meanlog = meanlog, sdlog = sdlog, shape1 = shape1, shape2 = shape2), convergence = convergence, n = n, np = 4)
  return(return_list)
}
