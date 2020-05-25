#' The Fréchet distribution
#'
#' Density, distribution function, quantile function, raw moments and random generation for the Fréchet distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param shape,scale Shape and scale of the Fréchet distribution, defaults to 1.5 and 0.5 respectively.
#' @param r rth raw moment of the distribution
#' @param truncation lower truncation parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{\left(E[x^r|X \le y]\right)}, otherwise, \eqn{P[X > x]} \eqn{\left(E[x^r|X > y]\right)}
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) =\frac{shape}{scale}\left(\frac{\omega}{scale}\right)^{-1-shape} e^{-\left(\frac{\omega}{scale}\right)^{-shape}}, \qquad F_X(x) = e^{-\left(\frac{\omega}{scale}\right)^{-shape}}}
#'
#'  The y-bounded r-th raw moment of the Fréchet distribution equals:
#'
#'  \deqn{ \mu^{r}_{y} = scale^{\sigma_s - 1} \left[1-\Gamma\left(1-\frac{\sigma_s - 1}{shape}, \left(\frac{y}{scale}\right)^{-shape} \right)\right],  \qquad shape>r}
#'
#' @return dfrechet returns the density, pfrechet the distribution function, qfrechet the quantile function, mfrechet the rth moment of the distribution and rfrechet generates random deviates.
#'
#'  The length of the result is determined by n for rfrechet, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#'
#' ## Frechet density
#' plot(x = seq(0, 5, length.out = 100), y = dfrechet(x = seq(0, 5, length.out = 100),
#' shape = 1, scale = 1))
#' plot(x = seq(0, 5, length.out = 100), y = dfrechet(x = seq(0, 5, length.out = 100),
#' shape = 2, scale = 1))
#' plot(x = seq(0, 5, length.out = 100), y = dfrechet(x = seq(0, 5, length.out = 100),
#' shape = 3, scale = 1))
#' plot(x = seq(0, 5, length.out = 100), y = dfrechet(x = seq(0, 5, length.out = 100),
#' shape = 3, scale = 2))
#'
#' ## frechet is also called the inverse weibull distribution, which is available in the stats package
#' pfrechet(q = 5, shape = 2, scale = 1.5)
#' 1 - pweibull(q = 1 / 5, shape = 2, scale = 1 / 1.5)
#'
#' ## Demonstration of log functionality for probability and quantile function
#' qfrechet(pfrechet(2, log.p = TRUE), log.p = TRUE)
#'
#' ## The zeroth truncated moment is equivalent to the probability function
#' pfrechet(2)
#' mfrechet(truncation = 2)
#'
#' ## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' x <- rfrechet(1e5, scale = 1)
#'
#' mean(x)
#' mfrechet(r = 1, lower.tail = FALSE, scale = 1)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mfrechet(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE, scale = 1)
#' @name frechet

NULL

#' @rdname frechet
#' @export

dfrechet <- function(x, shape = 1.5, scale = 0.5, log = FALSE) {
  d <- log(shape / scale) + (-1 - shape) * log(x / scale) + log(exp(-(x / scale)^(-shape)))

  if (!log) d <- exp(d)
  return(d)
}

#' @rdname frechet
#' @export

pfrechet <- function(q, shape = 1.5, scale = 0.5, log.p = FALSE, lower.tail = TRUE) {
  p <- exp(-(q / scale)^(-shape))

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

#' @rdname frechet
#' @export

qfrechet <- function(p, shape = 1.5, scale = 0.5, log.p = FALSE, lower.tail = TRUE) {
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

  q <- scale * (-log(p))^(-1 / shape)

  return(q)
}

#' @rdname frechet
#' @export

mfrechet <- function(r = 0, truncation = 0, shape = 1.5, scale = 0.5, lower.tail = TRUE) {
  truncation <- as.vector(truncation)

  if (shape > r) {
    m <- ((scale^(r)) * (1 - pgamma(((truncation / scale)^(-shape)), (1 - (r) / shape), lower.tail = FALSE)) * gamma((1 - (r) / shape)))

    notrunc <- ((scale^(r)) * (1 - pgamma(((0 / scale)^(-shape)), (1 - (r) / shape), lower.tail = FALSE)) * gamma((1 - (r) / shape)))

    if (lower.tail) {
      m <- notrunc - m
    }
  } else {
    m <- rep(NA, length(truncation))
  }

  return(m)
}

#' @rdname frechet
#' @export

rfrechet <- function(n, shape = 1.5, scale = 0.5) {
  r <- qfrechet(runif(n), shape = shape, scale = scale)
  return(r)
}

#' Fréchet coefficients after power-law transformation
#'
#' Coefficients of a power-law transformed Fréchet distribution
#' @param shape,scale Scale and shape of the Fréchet distribution, defaults to 1.5 and 0.5 respectively.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @details If the random variable x is Fréchet distributed with scale shape and shape scale, then the power-law transformed variable
#'
#'  \deqn{ y = ax^b }
#'
#'  is Fréchet distributed with scale \eqn{ \left( \frac{scale}{a}\right)^{\frac{1}{b}} } and shape \eqn{b*k}.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  pfrechet(3,shape=2,scale=1)
#'  coeff = frechet_plt(shape=2,scale=1,a=5,b=7)$coefficients
#'  pfrechet(5*3^7,shape=coeff[["shape"]],scale=coeff[["scale"]])
#'
#'  pfrechet(5*0.8^7,shape=2,scale=1)
#'  coeff = frechet_plt(shape=2,scale=1,a=5,b=7,inv=TRUE)$coefficients
#'  pfrechet(0.8,shape=coeff[["shape"]],scale=coeff[["scale"]])
#'
#' @export

frechet_plt <- function(shape = 1.5, scale = 0.5, a = 1, b = 1, inv = FALSE) {
  if (!inv) {
    b <- 1 / b
    a <- (1 / a)^b
  }

  newk <- b * shape
  news <- (scale / a)^(1 / b)

  return_list <- list(coefficients = c(shape = newk, scale = news))
  return(return_list)
}


#' Fréchet MLE
#'
#' Maximum likelihood estimation of the coefficients of the Fréchet distribution
#' @param x data vector
#' @param weights numeric vector for weighted MLE, should have the same length as data vector x
#' @param start named vector with starting values, default to c(shape=1.5,scale=0.5)
#' @param lower,upper Lower and upper bounds to the estimated shape parameter, defaults to 1e-10 and Inf respectively
#'
#' @return Returns a named list containing a
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' \item{convergence}{logical indicator of convergence}
#' \item{n}{Length of the fitted data vector}
#' \item{np}{Nr. of coefficients}
#' }
#'
#' x = rfrechet(1e3)
#'
#' ## Pareto fit with xmin set to the minium of the sample
#' frechet.mle(x=x)
#'
#' @export

frechet.mle <- function(x, weights = NULL, start = c(shape = 1.5, scale = 0.5), lower = c(1e-10, 1e-10), upper = c(Inf, Inf)) {
  if (is.null(weights)) {
    fnobj <- function(par, x, names) {
      names(par) <- names

      out <- -sum(dfrechet(x, shape = par["shape"], scale = par["scale"], log = TRUE))

      if (is.infinite(out)) {
        out <- 1e20
      }
      # if(is.nan(out)){out=1e20}

      return(out)
    }

    # out = optim(par = start, fn = fnobj, x=x)
    out <- nlminb(start = start, objective = fnobj, lower = lower, upper = upper, x = x, control = list(maxit = 1e5), names = names(start))

    shape <- out$par[1]
    scale <- out$par[2]
  } else {
    fnobj_w <- function(par, x, weights, names) {
      names(par) <- names

      out <- -sum(weights * dfrechet(x, shape = par["shape"], scale = par["scale"], log = TRUE))
      if (is.infinite(out)) {
        out <- 1e20
      }
      # if(is.nan(out)){out=1e20}

      return(out)
    }
    out <- nlminb(start, fnobj_w, x = x, weights = weights, lower = lower, upper = upper, control = list(maxit = 1e5), names = names(start))
  }

  shape <- out$par["shape"]
  scale <- out$par["scale"]

  convergence <- ifelse(out$convergence == 0, TRUE, FALSE)

  return_list <- list(coefficients = c(shape = shape, scale = scale), n = length(x), np = 2, convergence = convergence)
  return(return_list)
}
