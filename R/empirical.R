#' The empirical distribution
#'
#' Density, distribution function, quantile function, and raw moments for the empirical distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param r rth raw moment of the Pareto distribution
#' @param data data vector
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param truncation lower truncation parameter, defaults to NULL.
#' @param lower.tail logical; if TRUE (default), moments are \eqn{E[x^r|X \le y]}, otherwise, \eqn{E[x^r|X > y]}
#'
#' @details The density function is a standard Kernel density estimation for 1e6 equally spaced points. The cumulative Distribution Function:
#'
#'  \deqn{F_n(x) = \frac{1}{n}\sum_{i=1}^{n}I_{x_i \leq x}}
#'
#'  The y-bounded r-th raw moment of the empirical distribution equals:
#'
#'  \deqn{ \mu^{r}_{y} = \frac{1}{n}\sum_{i=1}^{n}I_{x_i \leq x}x^r}
#'
#' @return dempirical returns the density, pempirical the distribution function, qempirical the quantile function, mempirical gives the rth moment of the distribution or a function that allows to evaluate the rth moment of the distribution if truncation is NULL..
#'
#' @examples
#' #'
#' ## Generate random sample to work with
#' x <- rlnorm(1e5, meanlog = -0.5, sdlog = 0.5)
#'
#' ## Empirical density
#' plot(x = seq(0, 5, length.out = 100), y = dempirical(x = seq(0, 5, length.out = 100), data = x))
#'
#' # Compare empirical and parametric quantities
#' dlnorm(0.5, meanlog = -0.5, sdlog = 0.5)
#' dempirical(0.5, data = x)
#'
#' plnorm(0.5, meanlog = -0.5, sdlog = 0.5)
#' pempirical(0.5, data = x)
#'
#' qlnorm(0.5, meanlog = -0.5, sdlog = 0.5)
#' qempirical(0.5, data = x)
#'
#' mlnorm(r = 0, truncation = 0.5, meanlog = -0.5, sdlog = 0.5)
#' mempirical(r = 0, truncation = 0.5, data = x)
#'
#' mlnorm(r = 1, truncation = 0.5, meanlog = -0.5, sdlog = 0.5)
#' mempirical(r = 1, truncation = 0.5, data = x)
#'
#' ## Demonstration of log functionailty for probability and quantile function
#' quantile(x, 0.5, type = 1)
#' qempirical(p = pempirical(q = quantile(x, 0.5, type = 1), data = x, log.p = TRUE),
#' data = x, log.p = TRUE)
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' pempirical(q = quantile(x, 0.5, type = 1), data = x)
#' mempirical(truncation = quantile(x, 0.5, type = 1), data = x)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' mean(x)
#' mempirical(r = 1, data = x, truncation = 0, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mempirical(r = 1, data = x, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' #'
#' @name empirical

NULL

#' @rdname empirical
#' @export

dempirical <- function(x, data, log = FALSE) {
  d_fun <- approxfun(density(data, n = 1e6))
  d <- d_fun(x)
  if (log) {
    d <- log(d)
  }

  return(d)
}

#' @rdname empirical
#' @export

pempirical <- function(q, data, log.p = FALSE, lower.tail = TRUE) {
  data <- sort(data)
  n <- length(data)
  vals <- unique(data)
  p_fun <- approxfun(vals, cumsum(tabulate(match(data, vals))) / n,
    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered"
  )

  p <- p_fun(q)
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

#' @rdname empirical
#' @export

qempirical <- function(p, data, lower.tail = TRUE, log.p = FALSE) {
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

  # Type 1 is pure discontinuous evaluation without interpolation or smoothing
  q <- as.numeric(quantile(x = data, probs = p, type = 1))

  return(q)
}

#' @rdname empirical
#' @export

mempirical <- function(r = 0, data, truncation = NULL, lower.tail = TRUE) {
  x <- sort(data)
  n <- length(x)
  vals <- unique(x)
  truncation <- as.vector(truncation)

  m_fun <- approxfun(vals, (mean(x^r) - cumsum(tabulate(match(x, vals)) * vals^r) / n),
    method = "linear", yleft = mean(x^r), yright = 0, f = 0, ties = "ordered"
  )

  if (lower.tail) {
    m_fun <- approxfun(vals, (cumsum(tabulate(match(x, vals)) * vals^r) / n),
      method = "linear", yleft = 0, yright = mean(x^r), f = 0, ties = "ordered"
    )
  }

  if (is.null(truncation)) {
    return(m_fun)
  } else {
    return(m_fun(truncation))
  }



  # if(lower.tail){m = m_fun(0) - m}


  # if(!value){
  #   return(m_fun)
  # }else{
  #   return(moment(truncation))
  # }
}
