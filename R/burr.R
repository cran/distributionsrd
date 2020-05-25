#' The Burr distribution
#'
#' Density, distribution function, quantile function, raw moments and random generation for the Burr distribution, also known as the Burr Type XII distribution or the Singh-Maddala distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param shape1,shape2,scale Shape1, shape2 and scale of the Burr distribution, defaults to 2, 1 and 0.5.
#' @param r rth raw moment of the distribution
#' @param truncation lower truncation parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{(E[x^r|X \le y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) =\frac{\frac{kc}{scale}(\frac{\omega}{scale})^{shape2-1}}{(1+(\frac{\omega}{scale})^shape2)^{shape1+1}}, \qquad F_X(x) = 1-\frac{1}{(1+(\frac{\omega}{scale})^shape2)^shape1}}
#'
#'  The y-bounded r-th raw moment of the Fréchet distribution equals:
#'
#'  \deqn{scale^{r} shape1 [\boldsymbol{B}(\frac{r}{shape2} +1,shape1-\frac{r}{shape2}) - \boldsymbol{B}(\frac{(\frac{y}{scale})^{shape2}}{1+(\frac{y}{scale})^{shape2}};\frac{r}{shape2} +1,shape1-\frac{r}{shape2})],  \qquad shape2>r, kc>r}
#'
#' @return dburr returns the density, pburr the distribution function, qburr the quantile function, mburr the rth moment of the distribution and rburr generates random deviates.
#'
#'  The length of the result is determined by n for rburr, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#'
#' ## Burr density
#' plot(x = seq(0, 5, length.out = 100), y = dburr(x = seq(0, 5, length.out = 100)))
#' plot(x = seq(0, 5, length.out = 100), y = dburr(x = seq(0, 5, length.out = 100), shape2 = 3))
#'
#' ## Demonstration of log functionality for probability and quantile function
#' qburr(pburr(2, log.p = TRUE), log.p = TRUE)
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' pburr(2)
#' mburr(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a
#' #(truncated) random sample, for large enough samples.
#' x <- rburr(1e5, shape2 = 3)
#'
#' mean(x)
#' mburr(r = 1, shape2 = 3, lower.tail = FALSE)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mburr(r = 1, shape2 = 3, truncation = quantile(x, 0.1), lower.tail = FALSE)
#' @name burr

NULL

#' @rdname burr
#' @export

dburr <- function(x, shape1 = 2, shape2 = 1, scale = 0.5, log = FALSE) {
  d <- (log((shape1 * shape2) / scale) + (shape2 - 1) * log(x / scale)) - (shape1 + 1) * log(1 + (x / scale)^shape2)

  if (!log) d <- exp(d)
  return(d)
}

#' @rdname burr
#' @export

pburr <- function(q, shape1 = 2, shape2 = 1, scale = 0.5, log.p = FALSE, lower.tail = TRUE) {
  p <- 1 - 1 / ((1 + (q / scale)^shape2)^shape1)

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

#' @rdname burr
#' @export

qburr <- function(p, shape1 = 2, shape2 = 1, scale = 0.5, log.p = FALSE, lower.tail = TRUE) {
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

  q <- (((1 - p)^(-1 / shape1) - 1)^(1 / shape2)) * scale

  return(q)
}

#' @rdname burr
#' @export

mburr <- function(r = 0, truncation = 0, shape1 = 2, shape2 = 1, scale = 0.5, lower.tail = TRUE) {
  if ((shape2 > r) & (shape1 * shape2 > r)) {
    m <- (scale^(r) * shape1 * (beta(((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2))) - scale^(r) * shape1 * (pbeta((((truncation / scale)^shape2) / (1 + ((truncation / scale)^shape2))), ((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2)) * (beta(((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2)))))
    m[is.infinite(truncation)] <- (scale^(r) * shape1 * (beta(((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2))) - scale^(r) * shape1 * (pbeta(Inf, ((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2)) * (beta(((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2)))))

    notrunc <- (scale^(r) * shape1 * (beta(((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2))) - scale^(r) * shape1 * (pbeta((((0 / scale)^shape2) / (1 + ((0 / scale)^shape2))), ((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2)) * (beta(((shape2 + (r)) / shape2), ((shape2 * shape1 - (r)) / shape2)))))

    if (lower.tail) {
      m <- notrunc - m
    }
  } else {
    m <- rep(NA, length(truncation))
  }

  return(m)
}

#' @rdname burr
#' @export

rburr <- function(n, shape1 = 2, shape2 = 1, scale = 0.5) {
  r <- qburr(runif(n), shape1 = shape1, shape2 = shape2, scale = scale)
  return(r)
}

#' Burr coefficients after power-law transformation
#'
#' Coefficients of a power-law transformed Burr distribution
#' @param shape1,shape2,scale Shape1, shape2 and scale of the Burr distribution, defaults to 2, 1 and 1 respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable x is Burr distributed with scale shape and shape scale, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Burr distributed with shape1 \eqn{shape1}, shape2 \eqn{b*shape2} and scale \eqn{ ( \frac{scale}{a})^{\frac{1}{b}} }.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  pburr(3,shape1=2,shape2=3,scale=1)
#'  coeff = burr_plt(shape1=2,shape2=3,scale=1,a=5,b=7)$coefficients
#'  pburr(5*3^7,shape1=coeff[["shape1"]],shape2=coeff[["shape2"]],scale=coeff[["scale"]])
#'
#'  pburr(5*0.9^7,shape1=2,shape2=3,scale=1)
#'  coeff = burr_plt(shape1=2,shape2=3,scale=1,a=5,b=7, inv=TRUE)$coefficients
#'  pburr(0.9,shape1=coeff[["shape1"]],shape2=coeff[["shape2"]],scale=coeff[["scale"]])
#'
#'  ## Comparing the first moments and sample means of power-law transformed variables for large enough samples
#'  x = rburr(1e5,shape1=2,shape2=3,scale=1)
#'  coeff = burr_plt(shape1=2,shape2=3,scale=1,a=2,b=0.5)$coefficients
#'  y = rburr(1e5,shape1=coeff[["shape1"]],shape2=coeff[["shape2"]],scale=coeff[["scale"]])
#'  mean(2*x^0.5)
#'  mean(y)
#'  mburr(r=1,shape1=coeff[["shape1"]],shape2=coeff[["shape2"]],scale=coeff[["scale"]],lower.tail=FALSE)
#'
#' @export

burr_plt <- function(shape1 = 2, shape2 = 1, scale = 0.5, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  newc <- b * shape2
  news <- (scale / a)^(1 / b)

  return_list <- list(coefficients = c(shape1 = shape1, shape2 = newc, scale = news))
  return(return_list)
}
