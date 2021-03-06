% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invpareto.R
\name{clauset.xmax}
\alias{clauset.xmax}
\title{Pareto scale determination à la Clauset}
\usage{
clauset.xmax(x, q = 1)
}
\arguments{
\item{x}{data vector}

\item{q}{Percentage of data to search over (starting from the smallest values)}
}
\value{
Returns a named list containing a
\describe{
\item{coefficients}{Named vector of coefficients}
\item{KS}{Minimum Kolmogorov-Smirnov distance}
\item{n}{Number of observations in the Inverse Pareto tail}
\item{coeff.evo}{Evolution of the Inverse Pareto shape parameter over the iterations}
}
}
\description{
This method determines the optimal scale parameter of the Inverse Pareto distribution using the iterative method \insertCite{clauset2009power}{distributionsrd} that minimizes the Kolmogorov-Smirnov distance.
}
\examples{

## Determine cuttof from compostie InvPareto-Lognormal distribution using Clauset's method
dist <- c("invpareto", "lnorm")
coeff <- c(coeff1.k = 1.5, coeff2.meanlog = 1, coeff2.sdlog = 0.5)
x <- rcomposite(1e3, dist = dist, coeff = coeff)
out <- clauset.xmax(x = x)
out$coefficients
coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))$coeff1

## Speed up method by considering values above certain quantile only
dist <- c("invpareto", "lnorm")
coeff <- c(coeff1.k = 1.5, coeff2.meanlog = 1, coeff2.sdlog = 0.5)
x <- rcomposite(1e3, dist = dist, coeff = coeff)
out <- clauset.xmax(x = x, q = 0.5)
out$coefficients
coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))$coeff1
}
\references{
\insertAllCited{}
}
