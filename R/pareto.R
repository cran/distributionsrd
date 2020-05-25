#' The Pareto distribution
#'
#' Density, distribution function, quantile function, raw moments and random generation for the Pareto distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param xmin,k Scale and shape of the Pareto distribution, defaults to 1 and 2 respectively.
#' @param r rth raw moment of the Pareto distribution
#' @param truncation lower truncation parameter, defaults to xmin
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{ (E[x^r|X \le y] )}, otherwise, \eqn{P[X > x]} \eqn{ (E[x^r|X > y] )}
#' @param na.rm Removes values that fall outside the support of the distribution
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) = \frac{kx_{min}^{k}}{x^{k+1}}, \qquad F_X(x) = 1-(\frac{x_{min} }{x})^{k}}
#'
#'  The y-bounded r-th raw moment of the Pareto distribution equals:
#'
#'  \deqn{ \mu^{r}_{y} = k x_{min}^k \frac{- y^{r-k} }{r-k},  \qquad k>r}
#'
#' @return dpareto returns the density, ppareto the distribution function, qpareto the quantile function, mpareto the rth moment of the distribution and rpareto generates random deviates.
#'
#'  The length of the result is determined by n for rpareto, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#'
#' ## Pareto density
#' plot(x = seq(1, 5, length.out = 100), y = dpareto(x = seq(1, 5, length.out = 100), k = 2, xmin = 1))
#'
#' ## Pareto relates to the exponential distribution available in the stats package
#' ppareto(q = 5, k = 2, xmin = 3)
#' pexp(q = log(5 / 3), rate = 2)
#'
#' ## Demonstration of log functionality for probability and quantile function
#' qpareto(ppareto(2, log.p = TRUE), log.p = TRUE)
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' ppareto(2)
#' mpareto(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rpareto(1e5)
#'
#' mean(x)
#' mpareto(r = 1, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mpareto(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name pareto

NULL

#' @rdname pareto
#' @export

dpareto <- function(x, k = 2, xmin = 1, log = FALSE, na.rm = FALSE) {
  d <- (k) / (xmin^(-k)) * x^(-k - 1)
  if (log) {
    d <- log(d)
  }
  d[x < xmin] <- NA
  if (na.rm) {
    d <- d[x >= xmin]
  }
  return(d)
}

#' @rdname pareto
#' @export

ppareto <- function(q, k = 2, xmin = 1, lower.tail = TRUE, log.p = FALSE, na.rm = FALSE) {
  p <- 1 - (q^(-k)) / (xmin^(-k))

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

  p[(q < xmin)] <- NA
  if (na.rm) {
    p <- p[q >= xmin]
  }
  return(p)
}

#' @rdname pareto
#' @export

qpareto <- function(p, k = 2, xmin = 1, lower.tail = TRUE, log.p = FALSE) {
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

  q <- ((1 - p) * (xmin^(-k)))^(-1 / k)
  q[(q < xmin)] <- NA
  return(q)
}

#' @rdname pareto
#' @export

mpareto <- function(r = 0, truncation = xmin, k = 2, xmin = 1, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  if (any(truncation < xmin)) {
    truncation[truncation < xmin] <- xmin
    warning("Truncation point lies below Pareto scale parameter.")
  }

  truncation <- as.vector(truncation)
  if (k > r) {
    moment <- (k * xmin^k) / (r - k) * (-truncation^(r - k))

    notrunc <- (k * xmin^k) / (r - k) * (-xmin^(r - k))

    if (lower.tail) {
      moment <- notrunc - moment
    }
  } else {
    moment <- rep(NA, length(truncation))
  }



  return(moment)
}

#' @rdname pareto
#' @export

rpareto <- function(n, k = 2, xmin = 1) {
  r <- qpareto(runif(n), k = k, xmin = xmin)
  return(r)
}

#' Pareto coefficients after power-law transformation
#'
#' Coefficients of a power-law transformed Pareto distribution
#' @param xmin,k Scale and shape of the Pareto distribution, defaults to 1 and 2 respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable x is Pareto-distributed with scale xmin and shape k, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Pareto distributed with scale \eqn{ ( \frac{xmin}{a})^{\frac{1}{b}} } and shape \eqn{b*k}.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#' @examples
#'
#' ## Comparing probabilites of power-law transformed transformed variables
#' ppareto(3, k = 2, xmin = 2)
#' coeff <- pareto_plt(xmin = 2, k = 2, a = 5, b = 7)$coefficients
#' ppareto(5 * 3^7, k = coeff[["k"]], xmin = coeff[["xmin"]])
#'
#' ppareto(5 * 0.9^7, k = 2, xmin = 2)
#' coeff <- pareto_plt(xmin = 2, k = 2, a = 5, b = 7, inv = TRUE)$coefficients
#' ppareto(0.9, k = coeff[["k"]], xmin = coeff[["xmin"]])
#'
#' ## Comparing the first moments and sample means of power-law transformed variables for
#' #large enough samples
#' x <- rpareto(1e5, k = 2, xmin = 2)
#' coeff <- pareto_plt(xmin = 2, k = 2, a = 2, b = 0.5)$coefficients
#' y <- rpareto(1e5, k = coeff[["k"]], xmin = coeff[["xmin"]])
#' mean(2 * x^0.5)
#' mean(y)
#' mpareto(r = 1, k = coeff[["k"]], xmin = coeff[["xmin"]], lower.tail = FALSE)
#' @export

pareto_plt <- function(xmin = 1, k = 2, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  newk <- b * k
  newxmin <- (xmin / a)^(1 / b)

  return_list <- list(coefficients = c(xmin = newxmin, k = newk))
  return(return_list)
}


#' Pareto MLE
#'
#' Maximum likelihood estimation of the Pareto shape parameter using the Hill estimator.
#' @param x data vector
#' @param xmin scale parameter of the Pareto distribution, set to min(x) if not provided
#' @param clauset Indicator variable for calculating the scale parameter using the clauset method, overrides provided xmin
#' @param q Percentage of data to search over (starting from the largest values), defaults to 0.
#' @param lower,upper Lower and upper bounds to the estimated shape parameter, defaults to 1e-10 and Inf respectively
#'
#' @details The Hill estimator equals
#'
#'  \deqn{ \hat{k} = \frac{1}{ \frac{1}{n} \sum_{i=1}^{n} log \frac{x_i}{x_{min}}} }
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
#' x <- rpareto(1e3, k = 2, xmin = 2)
#'
#' ## Pareto fit with xmin set to the minium of the sample
#' pareto.mle(x = x)
#'
#' ## Pareto fit with xmin set to its real value
#' pareto.mle(x = x, xmin = 2)
#'
#' ## Pareto fit with xmin determined by the Clauset method
#' pareto.mle(x = x, clauset = TRUE)
#' @export

pareto.mle <- function(x, xmin = NULL, clauset = FALSE, q = 0, lower = 1e-10, upper = Inf) {
  if (clauset) {
    xmin <- clauset.xmin(x = x, q = q)$coefficients[["xmin"]]
  } else if (is.null(xmin)) {
    xmin <- min(x)
  }

  x <- x[x >= xmin]

  n <- length(x)
  # Hill estimator
  k <- (1 / n * sum(log(x / xmin)))^(-1)
  # Restricted Hill estimator... equals shape parameter equal to the bounds if they are hit...
  if (k < lower) {
    k <- lower
  }
  if (k > upper) {
    k <- upper
  }

  return_list <- list(coefficients = c(xmin = xmin, k = k), convergence = TRUE, n = n, np = 2)
  return(return_list)
}

#' Pareto scale determination à la Clauset
#'
#' This method determines the optimal scale parameter of the Pareto distribution using the iterative method \insertCite{clauset2009power}{distributionsrd}that minimizes the Kolmogorov-Smirnov distance.
#' @param x data vector
#' @param q Percentage of data to search over (starting from the largest values)
#'
#' @return Returns a named list containing a
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' \item{KS}{Minimum Kolmogorov-Smirnov distance}
#' \item{n}{Number of observations in the Pareto tail}
#' \item{coeff.evo}{Evolution of the Pareto shape parameter over the iterations}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#'
#' ## Determine cuttof from compostie lognormal-Pareto distribution using Clauset's method
#' dist <- c("lnorm", "pareto")
#' coeff <- c(coeff1.meanlog = -0.5, coeff1.sdlog = 0.5, coeff2.k = 1.5)
#' x <- rcomposite(1e3, dist = dist, coeff = coeff)
#' out <- clauset.xmin(x = x)
#' out$coefficients
#' coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))$coeff2
#'
#' ## Speed up method by considering values above certain quantile only
#' dist <- c("lnorm", "pareto")
#' coeff <- c(coeff1.meanlog = -0.5, coeff1.sdlog = 0.5, coeff2.k = 1.5)
#' x <- rcomposite(1e3, dist = dist, coeff = coeff)
#' out <- clauset.xmin(x = x, q = 0.5)
#' out$coefficients
#' coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))$coeff2
#' @export

clauset.xmin <- function(x, q = 0) {

  # Disregard values lower than minimum value
  x <- x[x >= quantile(x, q, type = 1)]

  x <- sort(x)

  # possible xmins
  xmins <- unique(x)
  xmins <- xmins[-length(xmins)]

  # Compute fit for each possible xmin
  KS.list <- rep(0, length(xmins))
  pareto.list <- rep(0, length(xmins))

  # Loop
  for (i in 1:length(xmins)) {
    xmin <- xmins[i]
    z <- x[x >= xmin]

    # Fit Pareto to this part of the data
    pareto.coefficients <- pareto.mle(x = z, xmin = xmin)$coefficients
    pareto.list[i] <- pareto.coefficients["k"]

    # Get Kolmogor-Smirnoff statistic
    KS.list[i] <- nmad_test(x = z, r = 0, dist = "pareto", coeff = pareto.coefficients, stat = "max")
  }

  i <- min(which(KS.list == min(KS.list)))
  xmin <- xmins[i]
  z <- x[x >= xmin]
  n <- length(z)
  pareto.coefficients <- pareto.mle(x = z, xmin = xmin)$coefficients

  return_list <- list(coefficients = pareto.coefficients, KS = KS.list[i], ntail = n, coeff.evo = cbind(pareto.list, xmins))
  return(return_list)
}
