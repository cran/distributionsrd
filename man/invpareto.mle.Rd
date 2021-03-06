% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invpareto.R
\name{invpareto.mle}
\alias{invpareto.mle}
\title{Inverse Pareto MLE}
\usage{
invpareto.mle(x, xmax = NULL, clauset = FALSE, q = 1)
}
\arguments{
\item{x}{data vector}

\item{xmax}{scale parameter of the Inverse Pareto distribution, set to max(x) if not provided}

\item{clauset}{Indicator variable for calculating the scale parameter using the clauset method, overrides provided xmax}

\item{q}{Percentage of data to search over (starting from the smallest values), dafults to 1.}
}
\value{
Returns a named list containing a
\describe{
\item{coefficients}{Named vector of coefficients}
\item{convergence}{logical indicator of convergence}
\item{n}{Length of the fitted data vector}
\item{np}{Nr. of coefficients}
}
}
\description{
Maximum likelihood estimation of the Inverse Pareto shape parameter using the Hill estimator.
}
\details{
The Hill estimator equals

 \deqn{\hat{k} =- \frac{1}{\frac{1}{n}\sum_{i=1}^{n}log\frac{x_{max}}{x_i}}}
}
\examples{

x <- rinvpareto(1e3, k = 1.5, xmax = 5)

## Pareto fit with xmin set to the minium of the sample
invpareto.mle(x = x)

## Pareto fit with xmin set to its real value
invpareto.mle(x = x, xmax = 5)

## Pareto fit with xmin determined by the Clauset method
invpareto.mle(x = x, clauset = TRUE)
}
