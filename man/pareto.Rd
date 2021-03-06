% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pareto.R
\name{pareto}
\alias{pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{mpareto}
\alias{rpareto}
\title{The Pareto distribution}
\usage{
dpareto(x, k = 2, xmin = 1, log = FALSE, na.rm = FALSE)

ppareto(q, k = 2, xmin = 1, lower.tail = TRUE, log.p = FALSE, na.rm = FALSE)

qpareto(p, k = 2, xmin = 1, lower.tail = TRUE, log.p = FALSE)

mpareto(r = 0, truncation = xmin, k = 2, xmin = 1, lower.tail = TRUE)

rpareto(n, k = 2, xmin = 1)
}
\arguments{
\item{x, q}{vector of quantiles}

\item{xmin, k}{Scale and shape of the Pareto distribution, defaults to 1 and 2 respectively.}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}

\item{na.rm}{Removes values that fall outside the support of the distribution}

\item{lower.tail}{logical; if TRUE (default), probabilities (moments) are \eqn{P[X \le x]} \eqn{ (E[x^r|X \le y] )}, otherwise, \eqn{P[X > x]} \eqn{ (E[x^r|X > y] )}}

\item{p}{vector of probabilities}

\item{r}{rth raw moment of the Pareto distribution}

\item{truncation}{lower truncation parameter, defaults to xmin}

\item{n}{number of observations}
}
\value{
dpareto returns the density, ppareto the distribution function, qpareto the quantile function, mpareto the rth moment of the distribution and rpareto generates random deviates.

 The length of the result is determined by n for rpareto, and is the maximum of the lengths of the numerical arguments for the other functions.
}
\description{
Density, distribution function, quantile function, raw moments and random generation for the Pareto distribution.
}
\details{
Probability and Cumulative Distribution Function:

 \deqn{f(x) = \frac{kx_{min}^{k}}{x^{k+1}}, \qquad F_X(x) = 1-(\frac{x_{min} }{x})^{k}}

 The y-bounded r-th raw moment of the Pareto distribution equals:

 \deqn{ \mu^{r}_{y} = k x_{min}^k \frac{- y^{r-k} }{r-k},  \qquad k>r}
}
\examples{

## Pareto density
plot(x = seq(1, 5, length.out = 100), y = dpareto(x = seq(1, 5, length.out = 100), k = 2, xmin = 1))

## Pareto relates to the exponential distribution available in the stats package
ppareto(q = 5, k = 2, xmin = 3)
pexp(q = log(5 / 3), rate = 2)

## Demonstration of log functionality for probability and quantile function
qpareto(ppareto(2, log.p = TRUE), log.p = TRUE)

## The zeroth truncated moment is equivalent to the probability function
ppareto(2)
mpareto(truncation = 2)

## The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#for large enough samples.
x <- rpareto(1e5)

mean(x)
mpareto(r = 1, lower.tail = FALSE)

sum(x[x > quantile(x, 0.1)]) / length(x)
mpareto(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE)
}
