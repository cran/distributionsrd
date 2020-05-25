#' The Inverse Pareto distribution
#'
#' Density, distribution function, quantile function, raw moments and random generation for the Pareto distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param xmax,k Scale and shape of the Inverse Pareto distribution, defaults to 5 and 1.5 respectively.
#' @param r rth raw moment of the Inverse Pareto distribution
#' @param truncation lower truncation parameter, defaults to xmin
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{(E[x^r|X \le y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#' @param na.rm Removes values that fall outside the support of the distribution
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) =\frac{k x_{max}^{-k}}{x^{-k+1}}, \qquad F_X(x) = (\frac{x_{max} }{x})^{-k}}
#'
#'  The y-bounded r-th raw moment of the Inverse Pareto distribution equals:
#'
#'  \deqn{\mu^r_y =k\omega_{max}^{-k}\frac{\omega_{max}^{r+k}- y^{r+k}}{r+k}  }
#'
#' @return dinvpareto returns the density, pinvpareto the distribution function, qinvpareto the quantile function, minvpareto the rth moment of the distribution and rinvpareto generates random deviates.
#'
#'  The length of the result is determined by n for rinvpareto, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#'
#' ## Inverse invpareto density
#' plot(x = seq(0, 5, length.out = 100), y = dinvpareto(x = seq(0, 5, length.out = 100)))
#'
#' ## Demonstration of log functionality for probability and quantile function
#' qinvpareto(pinvpareto(2, log.p = TRUE), log.p = TRUE)
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' pinvpareto(2)
#' minvpareto(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rinvpareto(1e5)
#'
#' mean(x)
#' minvpareto(r = 1, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' minvpareto(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name invpareto

NULL

#' @rdname invpareto
#' @export

dinvpareto <- function(x, k = 1.5, xmax = 5, log = FALSE, na.rm = FALSE) {
  d <- (k * xmax^(-k)) / (x^(-k + 1))
  if (log) {
    d <- log(d)
  }
  d[(x > xmax)] <- NA
  if (na.rm) {
    d <- d[x <= xmax]
  }
  return(d)
}

#' @rdname invpareto
#' @export

pinvpareto <- function(q, k = 1.5, xmax = 5, lower.tail = TRUE, log.p = FALSE, log = FALSE, na.rm = FALSE) {
  p <- (q^(k)) / (xmax^(k))

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

  p[(q > xmax)] <- NA
  if (na.rm) {
    p <- p[q <= xmax]
  }
  return(p)
}

#' @rdname invpareto
#' @export

qinvpareto <- function(p, k = 1.5, xmax = 5, lower.tail = TRUE, log.p = FALSE) {
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

  q <- ((p) * (xmax^(k)))^(1 / k)
  q[(q > xmax)] <- NA
  return(q)
}

#' @rdname invpareto
#' @export

minvpareto <- function(r = 0, truncation = 0, k = 1.5, xmax = 5, lower.tail = TRUE) {
  if (any(truncation > xmax)) {
    truncation[truncation > xmax] <- xmax
    warning("Truncation point lies above Inverse Pareto scale parameter.")
  }

  truncation <- as.vector(truncation)

  m <- (k * xmax^(-k)) / (r + k) * (xmax^(r + k) - truncation^(r + k))
  notrunc <- (k * xmax^(-k)) / (r + k) * (xmax^(r + k) - 0^(r + k))

  if (lower.tail) {
    m <- notrunc - m
  }

  return(m)
}

#' @rdname invpareto
#' @export

rinvpareto <- function(n, k = 1.5, xmax = 5) {
  r <- qinvpareto(runif(n), k = k, xmax = xmax)
  return(r)
}

#' Inverse Pareto coefficients after power-law transformation
#'
#' Coefficients of a power-law transformed Inverse Pareto distribution
#' @param xmax,k Scale and shape of the Inverse Pareto distribution, defaults to 5 and 1.5 respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable x is Inverse Pareto-distributed with scale xmin and shape k, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Inverse Pareto distributed with scale \eqn{ ( \frac{xmin}{a})^{\frac{1}{b}} } and shape \eqn{b*k}.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  pinvpareto(3,k=2,xmax=5)
#'  coeff = invpareto_plt(xmax=5,k=2,a=5,b=7)$coefficients
#'  pinvpareto(5*3^7,k=coeff[["k"]],xmax=coeff[["xmax"]])
#'
#'  pinvpareto(5*0.9^7,k=2,xmax=5)
#'  coeff = invpareto_plt(xmax=5,k=2,a=5,b=7, inv=TRUE)$coefficients
#'  pinvpareto(0.9,k=coeff[["k"]],xmax=coeff[["xmax"]])
#'
#' @export

invpareto_plt <- function(xmax = 5, k = 1.5, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  newk <- b * k
  newxmax <- (xmax / a)^(1 / b)

  return_list <- list(coefficients = c(xmax = newxmax, k = newk))
  return(return_list)
}

#' Inverse Pareto MLE
#'
#' Maximum likelihood estimation of the Inverse Pareto shape parameter using the Hill estimator.
#' @param x data vector
#' @param xmax scale parameter of the Inverse Pareto distribution, set to max(x) if not provided
#' @param clauset Indicator variable for calculating the scale parameter using the clauset method, overrides provided xmax
#' @param q Percentage of data to search over (starting from the smallest values), dafults to 1.
#'
#' @details The Hill estimator equals
#'
#'  \deqn{\hat{k} =- \frac{1}{\frac{1}{n}\sum_{i=1}^{n}log\frac{x_{max}}{x_i}}}
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
#' x <- rinvpareto(1e3, k = 1.5, xmax = 5)
#'
#' ## Pareto fit with xmin set to the minium of the sample
#' invpareto.mle(x = x)
#'
#' ## Pareto fit with xmin set to its real value
#' invpareto.mle(x = x, xmax = 5)
#'
#' ## Pareto fit with xmin determined by the Clauset method
#' invpareto.mle(x = x, clauset = TRUE)
#' @export

invpareto.mle <- function(x, xmax = NULL, clauset = FALSE, q = 1) {
  if (clauset) {
    xmax <- clauset.xmax(x = x, q = q)$coefficients[["xmax"]]
  } else if (is.null(xmax)) {
    xmax <- max(x)
  }

  x <- x[x <= xmax]

  n <- length(x)
  # Hill estimator
  k <- (1 / n * sum(log(xmax / x)))^(-1)

  return_list <- list(coefficients = c(k = k, xmax = xmax), convergence = TRUE, n = n, np = 2)
  return(return_list)
}


#' Pareto scale determination à la Clauset
#'
#' This method determines the optimal scale parameter of the Inverse Pareto distribution using the iterative method \insertCite{clauset2009power}{distributionsrd} that minimizes the Kolmogorov-Smirnov distance.
#' @param x data vector
#' @param q Percentage of data to search over (starting from the smallest values)
#'
#' @return Returns a named list containing a
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' \item{KS}{Minimum Kolmogorov-Smirnov distance}
#' \item{n}{Number of observations in the Inverse Pareto tail}
#' \item{coeff.evo}{Evolution of the Inverse Pareto shape parameter over the iterations}
#' }
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#'
#' ## Determine cuttof from compostie InvPareto-Lognormal distribution using Clauset's method
#' dist <- c("invpareto", "lnorm")
#' coeff <- c(coeff1.k = 1.5, coeff2.meanlog = 1, coeff2.sdlog = 0.5)
#' x <- rcomposite(1e3, dist = dist, coeff = coeff)
#' out <- clauset.xmax(x = x)
#' out$coefficients
#' coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))$coeff1
#'
#' ## Speed up method by considering values above certain quantile only
#' dist <- c("invpareto", "lnorm")
#' coeff <- c(coeff1.k = 1.5, coeff2.meanlog = 1, coeff2.sdlog = 0.5)
#' x <- rcomposite(1e3, dist = dist, coeff = coeff)
#' out <- clauset.xmax(x = x, q = 0.5)
#' out$coefficients
#' coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))$coeff1
#' @export

clauset.xmax <- function(x, q = 1) {

  # Disregard values higher than maximum value
  x <- x[x <= quantile(x, q, type = 1)]

  x <- sort(x)

  # possible xmaxs
  xmaxs <- unique(x)
  xmaxs <- xmaxs[-1]

  # Compute fit for each possible xmax
  KS.list <- rep(0, length(xmaxs))

  # Loop
  for (i in 1:length(xmaxs)) {
    xmax <- xmaxs[i]
    z <- x[x <= xmax]

    # Fit Pareto to this part of the data
    invpareto.coefficients <- invpareto.mle(x = z, xmax = xmax)$coefficients

    # Get Kolmogor-Smirnoff statistic
    KS.list[i] <- nmad_test(x = z, r = 0, dist = "invpareto", coeff = invpareto.coefficients, stat = "max")
  }

  i <- max(which(KS.list == min(KS.list)))
  xmax <- xmaxs[i]
  z <- x[x <= xmax]
  n <- length(z)
  invpareto.coefficients <- invpareto.mle(x = z, xmax = xmax)$coefficients

  return_list <- list(coefficients = invpareto.coefficients, KS = KS.list[i], ntail = n)
  return(return_list)
}
