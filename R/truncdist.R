#' Truncated distribution
#'
#' Density, distribution function, quantile function, raw moments and random generation for a truncated distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param r rth raw moment of the distribution
#' @param dist distribution to be truncated, defaults to lnorm
#' @param coeff list of parameters for the truncated distribution, defaults to list(meanlog=0,sdlog=1)
#' @param lowertrunc,uppertrunc lowertrunc- and uppertrunc truncation points, defaults to 0 and Inf respectively
#' @param truncation lowertrunc truncation parameter, defaults to 0.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{(E[x^r|X \le y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#'
#' @details Probability and Cumulative Distribution Function:
#'
#'  \deqn{f(x) = \frac{g(x)}{F(uppertrunc)-F(lowertrunc)}, \qquad F_X(x) = \frac{F(x)-F(lowertrunc)}{F(uppertrunc)-F(lowertrunc)}}
#'
#' @return dtruncdist gives the density, ptruncdist gives the distribution function, qtruncdist gives the quantile function, mtruncdist gives the rth moment of the distribution and rtruncdist generates random deviates.
#'
#'  The length of the result is determined by n for rpareto, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#'
#' ## Truncated lognormal density
#' plot(x = seq(0.5, 3, length.out = 100), y = dtruncdist(x = seq(0.5, 5, length.out = 100),
#' dist = "lnorm", coeff = list(meanlog = 0.5, sdlog = 0.5), lowertrunc = 0.5, uppertrunc = 5))
#' lines(x = seq(0, 6, length.out = 100), y = dlnorm(x = seq(0, 6, length.out = 100),
#' meanlog = 0.5, sdlog = 0.5))
#'
#' # Compare quantities
#' dtruncdist(0.5)
#' dlnorm(0.5)
#' dtruncdist(0.5, lowertrunc = 0.5, uppertrunc = 3)
#'
#' ptruncdist(2)
#' plnorm(2)
#' ptruncdist(2, lowertrunc = 0.5, uppertrunc = 3)
#'
#' qtruncdist(0.25)
#' qlnorm(0.25)
#' qtruncdist(0.25, lowertrunc = 0.5, uppertrunc = 3)
#'
#' mtruncdist(r = 0, truncation = 2)
#' mlnorm(r = 0, truncation = 2, meanlog = 0, sdlog = 1)
#' mtruncdist(r = 0, truncation = 2, lowertrunc = 0.5, uppertrunc = 3)
#'
#' mtruncdist(r = 1, truncation = 2)
#' mlnorm(r = 1, truncation = 2, meanlog = 0, sdlog = 1)
#' mtruncdist(r = 1, truncation = 2, lowertrunc = 0.5, uppertrunc = 3)
#' @name truncdist
#'

NULL

#' @rdname truncdist
#' @export

dtruncdist <- function(x, dist = c("lnormtrunc"), coeff = list(meanlog = 0, sdlog = 1), lowertrunc = 0, uppertrunc = Inf, log = FALSE) {
  # dist <- str_remove(dist, "trunc")
  dist <- gsub("trunc", "", dist)

  x <- x[(x >= lowertrunc & x <= uppertrunc)]

  dnum <- do.call(paste0("d", dist), args = c(x = list(x), coeff, list(log = TRUE)))
  ddenom <- log(do.call(paste0("p", dist), args = c(q = list(uppertrunc), coeff)) - do.call(paste0("p", dist), args = c(q = list(lowertrunc), coeff)))


  d <- dnum - ddenom

  if (!log) {
    d <- exp(d)
  }
  return(d)
}


#' @rdname truncdist
#' @export

ptruncdist <- function(q, dist = c("lnormtrunc"), coeff = list(meanlog = 0, sdlog = 1), lowertrunc = 0, uppertrunc = Inf, log.p = FALSE, lower.tail = TRUE) {
  dist <- gsub("trunc", "", dist)

  q <- q[(q >= lowertrunc & q <= uppertrunc)]

  p <- (do.call(paste0("p", dist), args = c(q = list(q), coeff)) - do.call(paste0("p", dist), args = c(q = list(lowertrunc), coeff)))
  p <- p / (do.call(paste0("p", dist), args = c(q = list(uppertrunc), coeff)) - do.call(paste0("p", dist), args = c(q = list(lowertrunc), coeff)))

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

#' @rdname truncdist
#' @export

qtruncdist <- function(p, dist = c("lnormtrunc"), coeff = list(meanlog = 0, sdlog = 1), lowertrunc = 0, uppertrunc = Inf, lower.tail = TRUE, log.p = FALSE) {
  dist <- gsub("trunc", "", dist)


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

  q <- (do.call(paste0("p", dist), args = c(q = list(uppertrunc), coeff)) - do.call(paste0("p", dist), args = c(q = list(lowertrunc), coeff)))
  int <- (p * q + do.call(paste0("p", dist), args = c(q = list(lowertrunc), coeff)))
  q <- do.call(paste0("q", dist), args = c(p = list(int), coeff))

  return(q)
}

#' @rdname truncdist
#' @export

mtruncdist <- function(r, truncation = 0, dist = c("lnormtrunc"), coeff = list(meanlog = 0, sdlog = 1), lowertrunc = 0, uppertrunc = Inf, lower.tail = TRUE) {
  dist <- gsub("trunc", "", dist)

  truncation <- truncation[(truncation >= lowertrunc & truncation <= uppertrunc)]

  momentnum <- do.call(paste0("m", dist), args = c(r = list(r), coeff, truncation = list(truncation))) - do.call(paste0("m", dist), args = c(r = list(r), coeff, truncation = list(lowertrunc)))
  momentdenom <- (do.call(paste0("p", dist), args = c(q = list(uppertrunc), coeff)) - do.call(paste0("p", dist), args = c(q = list(lowertrunc), coeff)))


  if (momentdenom == 0) {
    moment <- Inf
  } else {
    moment <- momentnum / momentdenom
  }

  notrunc_momentnum <- do.call(paste0("m", dist), args = c(r = list(r), coeff, truncation = list(uppertrunc))) - do.call(paste0("m", dist), args = c(r = list(r), coeff, truncation = list(lowertrunc)))
  notrunc <- notrunc_momentnum / momentdenom

  if (!lower.tail) {
    moment <- notrunc - moment
  }

  return(moment)
}

#' @rdname truncdist
#' @export

rtruncdist <- function(n, dist = c("lnormtrunc"), coeff = list(meanlog = 0, sdlog = 1), lowertrunc = 0, uppertrunc = Inf) {
  u <- runif(n)
  r <- qtruncdist(p = u, dist = dist, coeff = coeff, lowertrunc = lowertrunc, uppertrunc = uppertrunc)

  return(r)
}
